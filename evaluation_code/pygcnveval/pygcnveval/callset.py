import pyranges as pr
import numpy as np
from typing import FrozenSet, List, Generator
import vcf
import pandas as pd
import shutil
import gzip

from event import EventType, Event
from interval_collection import IntervalCollection
from interval import Interval


class CallsetMatrixView:
    """
    Matrix representation of a callset, that contains a corresponding callset genotypes in each of the
    (#samples) * (#intervals) positions. The object may also optionally contain a matrix of corresponding genotype
    qualities.
    """

    def __init__(self, sample_list: List, intervals_by_sample_call_matrix: np.ndarray,
                 intervals_by_sample_quality_matrix: np.ndarray):
        self.sample_list = sample_list
        self.intervals_by_sample_call_matrix = intervals_by_sample_call_matrix
        self.intervals_by_sample_quality_matrix = intervals_by_sample_quality_matrix

    def __eq__(self, other):
        if not isinstance(other, CallsetMatrixView):
            return False

        if set(self.sample_list) != set(set(other.sample_list)):
            return False
        for sample in self.sample_list:
            index_self = self.sample_list.index(sample)
            index_other = other.sample_list.index(sample)
            if (self.intervals_by_sample_call_matrix[index_self] != other.intervals_by_sample_call_matrix[index_other]).all() \
                    or (self.intervals_by_sample_quality_matrix[index_self] != other.intervals_by_sample_quality_matrix[index_other]).all():
                return False

        return True


class Callset:
    CALLSET_COLUMNS = ["Chromosome", "Start", "End", "Genotype", "NumBins", "Quality"]
    CALLSET_COLUMN_TYPES = {"Chromosome": "category", "Start": "int32", "End": "int32",
                            "Genotype": "object", "NumBins": "int32", "Quality": "int32"}

    def __init__(self, sample_to_pyrange_map: dict):
        self.sample_to_pyrange_map = sample_to_pyrange_map
        self.sample_set = frozenset(sample_to_pyrange_map.keys())

    def get_event_generator(self) -> Generator[Event, None, None]:
        """
        Get a generator that yields all events in the callset, i.e. for each site it will output all non-reference
        events containing all samples at this site

        :return: per event generator
        """
        for sample in self.sample_set:
            for k, df in self.sample_to_pyrange_map[sample]:
                for index, event in df.iterrows():
                    if event['Genotype'] == EventType.REF:
                        continue
                    interval = Interval(event['Chromosome'], event['Start'], event['End'])
                    event_type = event['Genotype']
                    call_attributes = {'NumBins': event['NumBins'], 'Quality': event['Quality']}

                    yield Event(interval, sample, event_type, call_attributes)

    def get_callset_matrix_view(self, interval_collection: IntervalCollection, sample_list: List[str]) -> CallsetMatrixView:
        """
        Construct a matrix of dimensions (number_of_samples x number_of_intervals), that contains genotypes for
        corresponding positions

        :param interval_collection: interval collection
        :param sample_list: sample list
        :return: CallsetMatrixView containing event by interval matrix
        """
        sample_by_interval_genotype_matrix = np.zeros((len(self.sample_set), len(interval_collection)), dtype=np.uint8)
        intervals_by_sample_quality_matrix = np.zeros((len(self.sample_set), len(interval_collection)), dtype=np.uint16)

        for index, sample in enumerate(sample_list):
            overlapping_intervals_genotypes = self.sample_to_pyrange_map[sample].intersect(interval_collection.pyrange)
            interval_collection_with_refs = interval_collection.pyrange.copy()\
                .assign("Genotype", lambda df: pd.Series([EventType.REF]).repeat(len(df)))\
                .assign("Quality", lambda df: pd.Series([0]).repeat(len(df)))
            non_overlapping_intervals = interval_collection_with_refs.subtract(overlapping_intervals_genotypes)
            overlapping_intervals_genotypes = overlapping_intervals_genotypes.drop(like="NumBins")
            joined_genotypes = pr.PyRanges(overlapping_intervals_genotypes.df.append(non_overlapping_intervals.df, ignore_index=True)).sort()
            sample_by_interval_genotype_matrix[index] = joined_genotypes.df["Genotype"].map(lambda e: e.value).values
            intervals_by_sample_quality_matrix[index] = joined_genotypes.df["Quality"].values

        return CallsetMatrixView(sample_list, sample_by_interval_genotype_matrix, intervals_by_sample_quality_matrix)

    def get_overlapping_events_for_sample(self, interval: Interval, sample: str) -> List[Event]:
        """
        Find all events in the callset that overlap a given interval for a given sample

        :param interval: interval
        :param sample: sample
        :return: list of all overlapping (non reference) events
        """
        sample_pyrange = self.sample_to_pyrange_map[sample]
        intersecting_calls = sample_pyrange[interval.chrom, interval.start:interval.end]
        return [Event(Interval(row['Chromosome'], row['Start'], row['End']), sample, row['Genotype'],
                      {'NumBins': row['NumBins'], 'Quality': row['Quality']})
                for index, row in intersecting_calls.df.iterrows() if row["Genotype"] != EventType.REF]


class TruthCallset(Callset):
    JOINT_CALLSET_COLUMNS = ['Chromosome', 'Start', 'End', 'Name', 'Genotype', 'Samples', 'Frequency', 'NumBins']
    JOINT_CALLSET_COLUMN_TYPES = {"Chromosome": "category", "Start": "int32", "End": "int32", "Name": "object",
                                  "Genotype": "object", "Samples": "object", "Frequency": "float64", "NumBins": "int64"}

    def __init__(self, truth_callset_pyrange: pr.PyRanges, sample_to_pyrange_map: dict):
        super().__init__(sample_to_pyrange_map)
        self.truth_callset_pyrange = truth_callset_pyrange
        self.sample_to_pyrange_map = sample_to_pyrange_map
        self.filtered_out_events = None

    @classmethod
    def read_in_callset(cls, truth_callset_bed_file: str, interval_collection: IntervalCollection, samples_to_keep: set = None):
        truth_callset_df = pd.read_csv(open(truth_callset_bed_file, 'r'), sep='\t', comment='#', names=TruthCallset.JOINT_CALLSET_COLUMNS)
        truth_callset_df['Samples'] = truth_callset_df['Samples'].map(lambda s: s.split(',')).map(frozenset)
        truth_callset_df['Genotype'] = truth_callset_df['Genotype'].map(lambda g: EventType.get_event_type_from_svtype(g))
        sample_set = frozenset.union(*truth_callset_df['Samples'].values)
        truth_callset_df['Frequency'] = truth_callset_df['Samples'].map(lambda s: len(s)/len(sample_set))
        truth_callset_df['NumBins'] = 0

        if samples_to_keep:
            truth_callset_df['Samples'] = truth_callset_df['Samples'].map(lambda s: s.intersection(samples_to_keep))
            sample_set = frozenset.union(*truth_callset_df['Samples'].values)

        # Annotate callset with count of overlapping bins
        truth_callset_df = truth_callset_df.astype(TruthCallset.JOINT_CALLSET_COLUMN_TYPES)
        truth_callset_pyrange = pr.PyRanges(truth_callset_df)
        truth_callset_with_coverage_pyrange = truth_callset_pyrange.coverage(interval_collection.pyrange,
                                                                             overlap_col="C", fraction_col="F")
        truth_callset_pyrange.NumBins = truth_callset_with_coverage_pyrange.df["C"]

        # Construct a map from samples to corresponding sample level callsets
        sample_to_pyrange_map = TruthCallset._construct_sample_to_pyrange_map(truth_callset_pyrange, sample_set)

        return cls(truth_callset_pyrange, sample_to_pyrange_map)

    @staticmethod
    def _construct_sample_to_pyrange_map(truth_callset_pyrange: pr.PyRanges, sample_set: FrozenSet):
        sample_to_pyrange_map = {}
        sample_to_events_list_map = {}
        for sample in sample_set:
            sample_to_events_list_map[sample] = []
        for k, df in truth_callset_pyrange:
            for index, truth_event in df.iterrows():
                for sample in truth_event['Samples']:
                    if sample in sample_set:
                        event = {'Chromosome': truth_event['Chromosome'], 'Start': truth_event['Start'], "End": truth_event['End'],
                                 "Genotype": truth_event['Genotype'], "NumBins": truth_event['NumBins'], "Quality": 0}
                        sample_to_events_list_map[sample].append(event)

        for sample in sample_set:
            events_df = pd.DataFrame(sample_to_events_list_map[sample])
            events_df = events_df.astype(Callset.CALLSET_COLUMN_TYPES)
            sample_to_pyrange_map[sample] = pr.PyRanges(events_df)
        return sample_to_pyrange_map

    def filter_out_uncovered_events(self, interval_collection: IntervalCollection,
                                    overall_reciprocal_overlap: float = 0.3):
        truth_callset_with_coverage = self.truth_callset_pyrange.coverage(interval_collection.pyrange)
        indices = truth_callset_with_coverage.df['FractionOverlaps'] > overall_reciprocal_overlap
        self.truth_callset_pyrange = truth_callset_with_coverage[indices].drop(['NumberOverlaps', 'FractionOverlaps'])
        self.filtered_out_events = truth_callset_with_coverage[~indices].drop(['NumberOverlaps', 'FractionOverlaps'])
        self.sample_to_pyrange_map = TruthCallset._construct_sample_to_pyrange_map(self.truth_callset_pyrange, self.sample_set)

    def get_multiallelic_regions(self) -> IntervalCollection:
        # TODO implement
        pass


class GCNVCallset(Callset):

    def __init__(self, sample_to_pyrange_map: dict):
        super().__init__(sample_to_pyrange_map)
        self.sample_to_pyrange_map = sample_to_pyrange_map
        self.sample_set = set(sample_to_pyrange_map.keys())

    @classmethod
    def read_in_callset(cls, gcnv_segment_vcfs: List[str]):
        # Handle gzip-ed VCF gcnv files
        new_gcnv_vcfs_list = []
        for vcf_file in gcnv_segment_vcfs:
            if vcf_file.endswith(".gz"):
                new_vcf_file = vcf_file[:-3]
                with open(new_vcf_file, 'wt') as f_out, gzip.open(vcf_file, 'rt') as f_in:
                    shutil.copyfileobj(f_in, f_out)
                    vcf_file = new_vcf_file
            new_gcnv_vcfs_list.append(vcf_file)
        gcnv_segment_vcfs = new_gcnv_vcfs_list

        sample_to_pyrange_map = {}
        for vcf_file in gcnv_segment_vcfs:
            vcf_reader = vcf.Reader(open(vcf_file, 'r'))
            assert len(vcf_reader.samples) == 1
            sample_name = vcf_reader.samples[0]
            events_df = pd.DataFrame(columns=Callset.CALLSET_COLUMNS)
            for record in vcf_reader:
                event_type = EventType.gcnv_call_to_event_type(int(record.genotype(sample_name)['GT']))
                event = (record.CHROM, int(record.POS), int(record.INFO['END']), event_type,
                         int(record.genotype(sample_name)['NP']), int(record.genotype(sample_name)['QS']))
                events_df.loc[len(events_df)] = event

            events_df = events_df.astype(Callset.CALLSET_COLUMN_TYPES)
            events_pr = pr.PyRanges(events_df)
            sample_to_pyrange_map[sample_name] = events_pr
        return cls(sample_to_pyrange_map)

