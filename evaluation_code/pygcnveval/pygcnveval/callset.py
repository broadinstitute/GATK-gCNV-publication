import pyranges as pr
from typing import FrozenSet, List
import vcf
import pandas as pd

from event import EventType, Event
from interval_collection import IntervalCollection
from interval import Interval


class TruthCallset:
    bed_columns = ['Chromosome', 'Start', 'End', 'Name', 'Genotype', 'Samples', 'Frequency', 'NumBins']
    bed_column_types = {"Chromosome": "category", "Start": "int32", "End": "int32", "Name": "object",
                        "Genotype": "object", "Samples": "object", "Frequency": "float64", "NumBins": "int64"}

    def __init__(self, truth_pyrange: pr.PyRanges, sample_set: FrozenSet):
        self.truth_callset_pyrange = truth_pyrange
        self.filtered_out_events = None
        self.sample_set = sample_set

    @classmethod
    def read_in_callset(cls, **kwargs):
        assert "truth_callset_bed" in kwargs
        assert "interval_collection" in kwargs
        truth_callset_bed = kwargs["truth_callset_bed"]
        interval_collection = kwargs["interval_collection"]
        truth_callset_df = pd.read_csv(open(truth_callset_bed, 'r'), sep='\t', comment='#', names=TruthCallset.bed_columns)
        truth_callset_df['Samples'] = truth_callset_df['Samples'].map(lambda s: s.split(',')).map(frozenset)
        truth_callset_df['Genotype'] = truth_callset_df['Genotype'].map(lambda g: EventType.get_event_type_from_svtype(g))
        sample_set = frozenset.union(*truth_callset_df['Samples'].values)
        truth_callset_df['Frequency'] = truth_callset_df['Samples'].map(lambda s: len(s)/len(sample_set))
        truth_callset_df['NumBins'] = 0

        # Annotate callset with count of overlapping bins
        truth_callset_df = truth_callset_df.astype(TruthCallset.bed_column_types)
        truth_callset_pyrange = pr.PyRanges(truth_callset_df)
        truth_callset_with_coverage_pyrange = truth_callset_pyrange.coverage(interval_collection.pyrange,
                                                                             overlap_col="C", fraction_col="F")
        truth_callset_pyrange.NumBins = truth_callset_with_coverage_pyrange.df["C"]
        return cls(truth_callset_pyrange, sample_set)

    def filter_out_uncovered_events(self, interval_collection: IntervalCollection,
                                    overall_reciprocal_overlap: float = 0.3, padding: int = 5000):
        #padded_interval_collection = interval_collection.pyrange.slack(padding)
        truth_callset_with_coverage = self.truth_callset_pyrange.coverage(interval_collection.pyrange)
        indices = truth_callset_with_coverage.df['FractionOverlaps'] > overall_reciprocal_overlap
        self.truth_callset_pyrange = truth_callset_with_coverage[indices].drop(['NumberOverlaps', 'FractionOverlaps'])
        self.filtered_out_events = truth_callset_with_coverage[~indices].drop(['NumberOverlaps', 'FractionOverlaps'])

    def get_event_generator(self):
        """
        Get a generator that yields all events in the callset, i.e. for each site it will output all events containing
        all samples at this site

        :return: per event generator
        """
        for k, df in self.truth_callset_pyrange:
            for index, truth_event in df.iterrows():
                interval = Interval(truth_event['Chromosome'], truth_event['Start'], truth_event['End'])
                event_type = truth_event['Genotype']
                call_attributes = {'NumBins': truth_event['NumBins']}

                for sample_name in truth_event['Samples']:
                    yield Event(interval, sample_name, event_type, call_attributes)

    def get_overlapping_events_for_sample(self, interval: Interval, sample: str):
        intersecting_calls = self.truth_callset_pyrange[interval.chrom, interval.start:interval.end]
        return [Event(Interval(row['Chromosome'], row['Start'], row['End']), sample, row['Genotype'],
                      {'NumBins': row['NumBins']}) for index, row in intersecting_calls.df.iterrows() if sample in row['Samples']]


class GCNVCallset:
    vcf_columns = ["Chromosome", "Start", "End", "Genotype", "NumBins", "QS"]
    vcf_column_types = {"Chromosome": "category", "Start": "int32", "End": "int32",
                        "Genotype": "object", "NumBins": "int32", "QS": "int32"}

    def __init__(self, sample_to_pyrange_map: dict):
        self.sample_to_pyrange_map = sample_to_pyrange_map
        self.sample_set = set(sample_to_pyrange_map.keys())

    @classmethod
    def read_in_callset(cls, **kwargs):
        assert "gcnv_segment_vcfs" in kwargs
        gcnv_segment_vcfs = kwargs["gcnv_segment_vcfs"]
        sample_to_pyrange_map = {}
        for vcf_file in gcnv_segment_vcfs:
            vcf_reader = vcf.Reader(open(vcf_file, 'r'))
            assert len(vcf_reader.samples) == 1
            sample_name = vcf_reader.samples[0]
            events_df = pd.DataFrame(columns=GCNVCallset.vcf_columns)
            for record in vcf_reader:
                event_type = EventType.gcnv_call_to_event_type(int(record.genotype(sample_name)['GT']))
                event = (record.CHROM, int(record.POS), int(record.INFO['END']), event_type,
                         int(record.genotype(sample_name)['NP']), int(record.genotype(sample_name)['QS']))
                events_df.loc[len(events_df)] = event

            events_df = events_df.astype(GCNVCallset.vcf_column_types)
            events_pr = pr.PyRanges(events_df)
            sample_to_pyrange_map[sample_name] = events_pr
        return cls(sample_to_pyrange_map)

    def get_event_generator(self):
        """
        Get a generator that yields all non reference events in the callset, i.e. for each site it will output all
        events containing all samples at this site

        :return: per event generator
        """
        for sample in self.sample_set:
            for k, df in self.sample_to_pyrange_map[sample]:
                for index, gcnv_event in df.iterrows():
                    if gcnv_event['Genotype'] == EventType.REF:
                        continue
                    interval = Interval(gcnv_event['Chromosome'], gcnv_event['Start'], gcnv_event['End'])
                    event_type = gcnv_event['Genotype']
                    call_attributes = {'NumBins': gcnv_event['NumBins']}

                    yield Event(interval, sample, event_type, call_attributes)

    def get_overlapping_events_for_sample(self, interval: Interval, sample: str) -> List[Event]:
        sample_pyrange = self.sample_to_pyrange_map[sample]
        intersecting_calls = sample_pyrange[interval.chrom, interval.start:interval.end]
        return [Event(Interval(row['Chromosome'], row['Start'], row['End']), sample, row['Genotype'], {'NumBins': row['NumBins']})
                for index, row in intersecting_calls.df.iterrows() if row["Genotype"] != EventType.REF]
