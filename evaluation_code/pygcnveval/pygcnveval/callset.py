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
        self.samples_by_intervals_call_matrix = intervals_by_sample_call_matrix
        self.samples_by_intervals_quality_matrix = intervals_by_sample_quality_matrix

    def __eq__(self, other):
        if not isinstance(other, CallsetMatrixView):
            return False

        if set(self.sample_list) != set(set(other.sample_list)):
            return False
        for sample in self.sample_list:
            index_self = self.sample_list.index(sample)
            index_other = other.sample_list.index(sample)
            if (self.samples_by_intervals_call_matrix[index_self] != other.samples_by_intervals_call_matrix[index_other]).all() \
                    or (self.samples_by_intervals_quality_matrix[index_self] != other.samples_by_intervals_quality_matrix[index_other]).all():
                return False

        return True


class Callset:
    CALLSET_COLUMNS = ["Chromosome", "Start", "End", "Genotype", "NumBins", "Quality"]
    CALLSET_COLUMN_TYPES = {"Chromosome": "category", "Start": "int32", "End": "int32",
                            "Genotype": "object", "NumBins": "int32", "Quality": "int32"}

    JOINT_CALLSET_COLUMNS = ['Chromosome', 'Start', 'End', 'Name', 'Genotype', 'Samples', 'Frequency', 'NumBins']
    JOINT_CALLSET_COLUMN_TYPES = {"Chromosome": "category", "Start": "int32", "End": "int32", "Name": "object",
                                  "Genotype": "object", "Samples": "object", "Frequency": "float64", "NumBins": "int64"}

    def __init__(self, sample_to_pyrange_map: dict, joint_callset: pr.PyRanges, interval_collection: IntervalCollection):
        self.sample_to_pyrange_map = sample_to_pyrange_map
        self.joint_callset = joint_callset
        self.sample_set = frozenset(sample_to_pyrange_map.keys())
        self.cached_interval_to_targets_map = self._create_cached_intervals_to_target_map(interval_collection)

    def get_event_generator(self, sample_list: List[str], min_quality_threshold=None) -> Generator[Event, None, None]:
        """
        Get a generator that yields all events in the callset, i.e. for each site it will output all non-reference
        events containing all samples at this site

        :return: per event generator
        """
        for sample in sample_list:
            assert sample in self.sample_set, "Queried sample (%s) is not in the callset." % sample
            for k, df in self.sample_to_pyrange_map[sample]:
                for index, event in df.iterrows():
                    if event['Genotype'] == EventType.REF:
                        continue
                    interval = Interval(event['Chromosome'], event['Start'], event['End'])
                    overlapping_target_set = self.cached_interval_to_targets_map[interval]
                    event_type = event['Genotype']
                    call_attributes = {'NumBins': event['NumBins'], 'Quality': event['Quality']}
                    if min_quality_threshold is None:
                        yield Event(interval, sample, event_type, call_attributes, overlapping_target_set)
                    if min_quality_threshold is not None and call_attributes['Quality'] >= min_quality_threshold:
                        yield Event(interval, sample, event_type, call_attributes, overlapping_target_set)

    def get_callset_matrix_view(self, interval_collection: IntervalCollection, sample_list: List[str]) -> CallsetMatrixView:
        """
        Construct a matrix of dimensions (number_of_samples x number_of_intervals), that contains genotypes for
        corresponding positions

        :param interval_collection: interval collection
        :param sample_list: sample subset from the callset
        :return: CallsetMatrixView containing event by interval matrix
        """
        sample_by_interval_genotype_matrix = np.zeros((len(sample_list), len(interval_collection)), dtype=np.uint8)
        sample_by_interval_quality_matrix = np.zeros((len(sample_list), len(interval_collection)), dtype=np.uint16)
        reciprocal_overlaps_matrix = np.zeros((len(sample_list), len(interval_collection)))

        indexed_interval_collection_pyrange = interval_collection.pyrange.insert(pd.Series(
            data=interval_collection.pyrange.df.index.tolist(), name="Index"))
        for i, sample in enumerate(sample_list):
            for k, df in self.sample_to_pyrange_map[sample]:
                for index, event in df.iterrows():
                    event_interval = Interval(event["Chromosome"], event["Start"], event["End"])
                    intersecting_intervals = indexed_interval_collection_pyrange[event_interval.chrom, event_interval.start:event_interval.end]

                    df_ = intersecting_intervals.df
                    if not df_.empty:
                        df_["Overlap"] = Interval.calculate_reciprocal_overlap(df_["Start"], df_["End"], event_interval.start, event_interval.end)
                        indices = df_["Index"][df_["Overlap"] > reciprocal_overlaps_matrix[i][df_["Index"]]].tolist()
                        reciprocal_overlaps_matrix[i][indices] = df_["Overlap"][df_["Index"].isin(indices)]

                        sample_by_interval_genotype_matrix[i, indices] = event["Genotype"].value
                        sample_by_interval_quality_matrix[i, indices] = event["Quality"]

        return CallsetMatrixView(sample_list, sample_by_interval_genotype_matrix, sample_by_interval_quality_matrix)

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
                      {'NumBins': row['NumBins'], 'Quality': row['Quality']},
                      self.cached_interval_to_targets_map[Interval(row['Chromosome'], row['Start'], row['End'])])
                for index, row in intersecting_calls.df.iterrows() if row["Genotype"] != EventType.REF]

    def get_callset_num_events_distribution(self) -> List[int]:
        """
        :return: list containing number of events for each sample
        """
        num_events_distr = []
        for sample in self.sample_set:
            num_events_distr.append(len(self.sample_to_pyrange_map[sample]))
        return num_events_distr

    def get_callset_event_size_distribution(self, max_allele_fraction: float = 1.0) -> List[int]:
        """
        :return: list containing event size (in overlapping bins) for each event
        """
        event_size_distr = []
        for event in self.get_event_generator(list(self.sample_set)):
            if 'Frequency' not in event.call_attributes:
                event_size_distr.append(event.call_attributes["NumBins"])
                continue
            if event.call_attributes['Frequency'] < max_allele_fraction:
                event_size_distr.append(event.call_attributes["NumBins"])
        return event_size_distr

    def _create_cached_intervals_to_target_map(self, interval_collection: IntervalCollection):
        joined_sample_level_callsets_df = pd.concat([self.sample_to_pyrange_map[s].df for s in self.sample_set])
        joined_sample_level_callset_pr = pr.PyRanges(joined_sample_level_callsets_df)

        joined_intervals = joined_sample_level_callset_pr.join(interval_collection.pyrange)
        joined_intervals_df = joined_intervals.df
        joined_intervals_df['Chromosome'] = joined_intervals_df['Chromosome'].astype('object')

        def construct_interval_to_targets_list_map(group: pd.DataFrame):
            return (Interval(group.iloc[0]["Chromosome"], group.iloc[0]['Start'], group.iloc[0]['End']),
                               (group["Chromosome"] + ":" + group["Start_b"].apply(str) + "-" + group["End_b"].apply(str)).tolist())

        s = joined_intervals_df.groupby(["Chromosome", "Start", "End"]).apply(construct_interval_to_targets_list_map)
        interval_to_targets_map = dict(s.tolist())
        for interval in interval_to_targets_map:
           interval_to_targets_map[interval] = set(interval_to_targets_map[interval])

        return interval_to_targets_map


class TruthCallset(Callset):
    JOINT_CALLSET_COLUMNS = ['Chromosome', 'Start', 'End', 'Name', 'Genotype', 'Samples', 'Frequency', 'NumBins']
    JOINT_CALLSET_COLUMN_TYPES = {"Chromosome": "category", "Start": "int32", "End": "int32", "Name": "object",
                                  "Genotype": "object", "Samples": "object", "Frequency": "float64", "NumBins": "int64"}

    def __init__(self, truth_callset_pyrange: pr.PyRanges, sample_to_pyrange_map: dict, interval_collection: IntervalCollection):
        super().__init__(sample_to_pyrange_map, truth_callset_pyrange, interval_collection)
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

        return cls(truth_callset_pyrange, sample_to_pyrange_map, interval_collection)

    # def get_callset_matrix_view(self, interval_collection: IntervalCollection, sample_list: List[str]) -> CallsetMatrixView:
    #     """
    #     We override superclass's method to boost performance, taking advantage of monolithic `truth_callset_pyrange`
    #     """
    #     def get_genotype_array_from_intersecting_events(interval_: Interval, samples_to_index_: dict,
    #                                                     intersecting_events_: pr.PyRanges):
    #         genotype_array = np.zeros(len(samples_to_index))
    #         reciprocal_overlap_array = np.zeros(len(samples_to_index))
    #         if intersecting_events_.empty:
    #             return genotype_array
    #         for k_, df_ in intersecting_events_:
    #             for index_, row_ in df_.iterrows():
    #                 ro = interval_.get_reciprocal_overlap(Interval(row_["Chromosome"], row_["Start"], row_["End"]))
    #                 for sample in row_["Samples"]:
    #                     if ro >= reciprocal_overlap_array[samples_to_index_[sample]]:
    #                         genotype_array[samples_to_index_[sample]] = row_["Genotype"].value
    #                         reciprocal_overlap_array[samples_to_index_[sample]] = ro
    #         return genotype_array
    #
    #     samples_to_index = {sample_list[i]: i for i in range(len(sample_list))}
    #     sample_by_interval_genotype_matrix = np.zeros((len(self.sample_set), len(interval_collection)), dtype=np.uint8)
    #     intervals_by_sample_quality_matrix = np.zeros((len(self.sample_set), len(interval_collection)), dtype=np.uint16)
    #     i = 0
    #     for k, df in interval_collection.pyrange:
    #         for index, row in df.iterrows():
    #             interval = Interval(row["Chromosome"], row["Start"], row["End"])
    #             intersecting_events = self.truth_callset_pyrange[interval.chrom, interval.start:interval.end]
    #             sample_by_interval_genotype_matrix[:, i] = get_genotype_array_from_intersecting_events(interval,
    #                                                                                                    samples_to_index,
    #                                                                                                    intersecting_events)
    #             i += 1
    #
    #     return CallsetMatrixView(sample_list, sample_by_interval_genotype_matrix, intervals_by_sample_quality_matrix)

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
                                    min_overlap_fraction: float = 0.00):
        truth_callset_with_coverage = self.truth_callset_pyrange.coverage(interval_collection.pyrange)
        indices = truth_callset_with_coverage.df['FractionOverlaps'] > min_overlap_fraction
        self.truth_callset_pyrange = truth_callset_with_coverage[indices].drop(['NumberOverlaps', 'FractionOverlaps'])
        self.filtered_out_events = truth_callset_with_coverage[~indices].drop(['NumberOverlaps', 'FractionOverlaps'])
        self.sample_to_pyrange_map = TruthCallset._construct_sample_to_pyrange_map(self.truth_callset_pyrange, self.sample_set)

    def subset_intervals_to_rare_regions(self, intervals: IntervalCollection, max_allelic_fraction: float) -> IntervalCollection:
        """
        Subset a given interval collection to those intervals who only overlap with rare events, as defined by
        max_allele_fraction

        :param intervals: interval collection to subset
        :param max_allelic_fraction: events below that threshold are considered rare
        :return: subset of rare intervals
        """
        # TODO Handle the case of an interval from collection overlapping multiple events. In that case we should only consider
        # allelic fraction of the event with the largest reciprocal overlap
        intervals_pr = intervals.pyrange
        overlaps = self.truth_callset_pyrange.intersect(intervals_pr)
        joined_pr = intervals_pr.join(overlaps, how="left")
        rare_intervals_pr = joined_pr.subset(lambda df: df["Frequency"] <= max_allelic_fraction)
        rare_intervals_pr = rare_intervals_pr.merge()
        rare_intervals_pr = pr.PyRanges(rare_intervals_pr.df)  # This is to resolve a dtypes issue
        return IntervalCollection(rare_intervals_pr)

    def get_callset_allele_size_distribution(self, max_allele_fraction: int = 1.0) -> List[int]:
        """
        :return: distribution of allele sizes (in bp)
        """
        allelic_size_distr = []
        for k, df in self.truth_callset_pyrange:
            for index, truth_event in df.iterrows():
                if truth_event["Frequency"] < max_allele_fraction:
                    allelic_size_distr.append(truth_event['NumBins'])
        return allelic_size_distr


class GCNVCallset(Callset):

    def __init__(self, sample_to_pyrange_map: dict, joint_callset: pr.PyRanges, interval_collection: IntervalCollection):
        super().__init__(sample_to_pyrange_map, joint_callset, interval_collection)

    @classmethod
    def read_in_callset(cls, gcnv_segment_vcfs: List[str], interval_collection: IntervalCollection):
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
            if len(events_df) > 120:
                continue
            events_pr = pr.PyRanges(events_df)
            sample_to_pyrange_map[sample_name] = events_pr
        return cls(sample_to_pyrange_map, None, interval_collection)


class XHMMCallset(Callset):
    def __init__(self, sample_to_pyrange_map: dict, joint_callset: pr.PyRanges, interval_collection: IntervalCollection):
        super().__init__(sample_to_pyrange_map, joint_callset,  interval_collection)

    @classmethod
    def read_in_callset(cls, xhmm_vcf: str, interval_collection: IntervalCollection):
        # Handle gzip-ed VCF gcnv files

        vcf_reader = vcf.Reader(open(xhmm_vcf, 'r'))
        sample_list = vcf_reader.samples

        events_df = pd.DataFrame(columns=Callset.CALLSET_COLUMNS)
        events_list = []
        for record in vcf_reader:
            del_call_samples = []
            dup_call_samples = []

            for sample in sample_list:
                event_type = EventType.gcnv_call_to_event_type(int(record.genotype(sample)['GT']))
                if event_type == EventType.DEL:
                    del_call_samples.append(sample)
                if event_type == EventType.DUP:
                    dup_call_samples.append(sample)

            chrom, start, end, num_bins = record.CHROM, int(record.POS), int(record.INFO['END']), int(record.genotype(sample)['NUMT'])
            if del_call_samples:
                events_list.append([chrom, start, end, EventType.DEL, num_bins, frozenset(del_call_samples)])
            if dup_call_samples:
                events_list.append([chrom, start, end, EventType.DEL, num_bins, frozenset(dup_call_samples)])

        events_df = pd.DataFrame(events_list, columns=Callset.JOINT_CALLSET_COLUMNS)
