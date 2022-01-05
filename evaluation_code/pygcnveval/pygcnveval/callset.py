import pyranges as pr
import numpy as np
from typing import FrozenSet, List, Generator, Optional, Callable
from abc import ABC, abstractmethod
import vcf
import pandas as pd
import shutil
import gzip
from collections import defaultdict

from event import EventType, Event, Allele
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


class Callset(ABC):
    CALLSET_COLUMNS = ["Chromosome", "Start", "End", "Genotype", "NumBins", "Quality", "Frequency"]
    CALLSET_COLUMN_TYPES = {"Chromosome": "category", "Start": "int32", "End": "int32",
                            "Genotype": "object", "NumBins": "int32", "Quality": "int32", "Frequency": "float64"}

    JOINT_CALLSET_COLUMNS = ['Chromosome', 'Start', 'End', 'Name', 'Genotype', 'Samples', 'Frequency', 'NumBins']
    JOINT_CALLSET_COLUMN_TYPES = {"Chromosome": "category", "Start": "int32", "End": "int32", "Name": "object",
                                  "Genotype": "object", "Samples": "object", "Frequency": "float64", "NumBins": "int64"}

    def __init__(self, sample_to_pyrange_map: dict, joint_callset: pr.PyRanges, interval_collection: IntervalCollection):
        self.sample_to_pyrange_map = sample_to_pyrange_map
        self.joint_callset = joint_callset
        self.sample_set = frozenset(sample_to_pyrange_map.keys())
        self.cached_interval_to_targets_map = self._create_cached_intervals_to_target_map(interval_collection)
        # TODO make cached interval map for both callset views in one go
        #self.cached_interval_to_targets_map = self._create_cached_intervals_to_target_map_joint(interval_collection)

    def get_event_generator(self, sample_list: List[str], event_filter: Callable[[Event], bool] = None) -> Generator[Event, None, None]:
        """
        Get a generator that yields all events in the callset, i.e. for each site it will output all non-reference
        events containing all samples at this site

        :return: per event generator
        """
        for sample in sample_list:
            assert sample in self.sample_set, "Queried sample (%s) is not in the callset." % sample
            for k, df in self.sample_to_pyrange_map[sample]:
                for index, row in df.iterrows():
                    if row['Genotype'] == EventType.REF:
                        continue

                    interval = Interval(row['Chromosome'], row['Start'], row['End'])

                    if interval not in self.cached_interval_to_targets_map:
                        continue

                    overlapping_target_set = self.cached_interval_to_targets_map[interval]

                    event_type = row['Genotype']
                    call_attributes = {'NumBins': row['NumBins'], 'Quality': row['Quality'], 'Frequency': row['Frequency']}
                    event = Event(interval, sample, event_type, call_attributes, overlapping_target_set)
                    if event_filter is None:
                        yield event
                    if event_filter is not None and event_filter(event):
                        yield event

    def get_allele_generator(self) -> Generator[Allele, None, None]:
        """
        TODO
        :return:
        """
        for k, df in self.joint_callset:
            for index, event in df.iterrows():
                if event['Genotype'] == EventType.REF:
                    continue

                interval = Interval(event['Chromosome'], event['Start'], event['End'])
                overlapping_target_set = self.cached_interval_to_targets_map[interval]
                allele_attributes = {'NumBins': event['NumBins'], 'Frequency': event['Frequency']}
                event_type = event['Genotype']

                sample_set = event['Samples']

                yield Allele(interval, event_type, allele_attributes, overlapping_target_set, sample_set)

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
                      {'NumBins': row['NumBins'], 'Quality': row['Quality'], 'Frequency': row['Frequency']},
                      self.cached_interval_to_targets_map[Interval(row['Chromosome'], row['Start'], row['End'])])
                for index, row in intersecting_calls.df.iterrows() if row["Genotype"] != EventType.REF and Interval(row['Chromosome'], row['Start'], row['End']) in self.cached_interval_to_targets_map]

    def get_overlapping_alleles(self, interval: Interval) -> List[Allele]:
        """
        TODO
        :param interval:
        :return:
        """
        intersecting_calls = self.joint_callset[interval.chrom, interval.start:interval.end]
        return [Allele(Interval(row['Chromosome'], row['Start'], row['End']),
                       row['Genotype'],
                       {'NumBins': row['NumBins'], 'Frequency': row['Frequency']},
                       self.cached_interval_to_targets_map[Interval(row['Chromosome'], row['Start'], row['End'])],
                       row['Samples'])
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

    def filter_out_uncovered_events_from_joint_callset(self, interval_collection: IntervalCollection,
                                                       min_overlap_fraction: float = 0.00):
        callset_with_coverage = self.joint_callset.coverage(interval_collection.pyrange)
        indices = callset_with_coverage.df['FractionOverlaps'] > min_overlap_fraction
        self.joint_callset = callset_with_coverage[indices].drop(['NumberOverlaps', 'FractionOverlaps'])
        self.sample_to_pyrange_map = Callset._construct_sample_to_pyrange_map(self.joint_callset, self.sample_set)
        #self.filtered_out_events = truth_callset_with_coverage[~indices].drop(['NumberOverlaps', 'FractionOverlaps'])
        #self.sample_to_pyrange_map = TruthCallset._construct_sample_to_pyrange_map(self.joint_callset, self.sample_set) TODO undo this comment

    def filter_out_uncovered_events(self, interval_collection: IntervalCollection,
                                                       min_overlap_fraction: float = 0.00): # TODO check reason for difference in results between this and filter_out_uncovered_events_from_joint_callset
        for sample in self.sample_set:
            callset_with_coverage = self.sample_to_pyrange_map[sample].coverage(interval_collection.pyrange)
            indices = callset_with_coverage.df['FractionOverlaps'] > min_overlap_fraction
            self.sample_to_pyrange_map[sample] = callset_with_coverage[indices].drop(['NumberOverlaps', 'FractionOverlaps'])
        #self.filtered_out_events = truth_callset_with_coverage[~indices].drop(['NumberOverlaps', 'FractionOverlaps'])
        #self.sample_to_pyrange_map = TruthCallset._construct_sample_to_pyrange_map(self.joint_callset, self.sample_set) TODO undo this comment

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

    def _create_cached_intervals_to_target_map_joint(self, interval_collection: IntervalCollection):
        joined_intervals = self.joint_callset.join(interval_collection.pyrange)
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

    @staticmethod
    def _construct_sample_to_pyrange_map(callset_pyrange: pr.PyRanges, sample_set: FrozenSet):
        sample_to_pyrange_map = {}
        sample_to_events_list_map = {}
        for sample in sample_set:
            sample_to_events_list_map[sample] = []
        for k, df in callset_pyrange:
            for index, truth_event in df.iterrows():
                for sample in truth_event['Samples']:
                    if sample in sample_set:
                        event = {'Chromosome': truth_event['Chromosome'], 'Start': truth_event['Start'], "End": truth_event['End'],
                                 "Genotype": truth_event['Genotype'], "NumBins": truth_event['NumBins'], "Quality": 0, "Frequency": truth_event['Frequency']}
                        sample_to_events_list_map[sample].append(event)

        for sample in sample_set:
            events_df = pd.DataFrame(sample_to_events_list_map[sample])
            events_df = events_df.astype(Callset.CALLSET_COLUMN_TYPES)
            sample_to_pyrange_map[sample] = pr.PyRanges(events_df)
        return sample_to_pyrange_map

    @abstractmethod
    def get_name(self):
        pass


class TruthCallset(Callset):
    JOINT_CALLSET_COLUMNS = ['Chromosome', 'Start', 'End', 'Name', 'Genotype', 'Samples', 'Frequency', 'NumBins']
    JOINT_CALLSET_COLUMN_TYPES = {"Chromosome": "category", "Start": "int32", "End": "int32", "Name": "object",
                                  "Genotype": "object", "Samples": "object", "Frequency": "float64", "NumBins": "int64"}

    def __init__(self, sample_to_pyrange_map: dict, truth_callset_pyrange: pr.PyRanges, interval_collection: IntervalCollection):
        super().__init__(sample_to_pyrange_map, truth_callset_pyrange, interval_collection)
        self.truth_callset_pyrange = truth_callset_pyrange
        self.sample_to_pyrange_map = sample_to_pyrange_map

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

        truth_callset_df = truth_callset_df[truth_callset_df['Samples'] != frozenset()]
        # Annotate callset with count of overlapping bins
        truth_callset_df = truth_callset_df.astype(TruthCallset.JOINT_CALLSET_COLUMN_TYPES)
        truth_callset_pyrange = pr.PyRanges(truth_callset_df)
        truth_callset_with_coverage_pyrange = truth_callset_pyrange.coverage(interval_collection.pyrange,
                                                                             overlap_col="C", fraction_col="F")
        truth_callset_pyrange.NumBins = truth_callset_with_coverage_pyrange.df["C"]

        # Construct a map from samples to corresponding sample level callsets
        sample_to_pyrange_map = Callset._construct_sample_to_pyrange_map(truth_callset_pyrange, sample_set)

        return cls(sample_to_pyrange_map, truth_callset_pyrange, interval_collection)

    # def filter_out_uncovered_events(self, interval_collection: IntervalCollection,
    #                                 min_overlap_fraction: float = 0.00):
    #     truth_callset_with_coverage = self.truth_callset_pyrange.coverage(interval_collection.pyrange)
    #     indices = truth_callset_with_coverage.df['FractionOverlaps'] > min_overlap_fraction
    #     self.truth_callset_pyrange = truth_callset_with_coverage[indices].drop(['NumberOverlaps', 'FractionOverlaps'])
    #     self.filtered_out_events = truth_callset_with_coverage[~indices].drop(['NumberOverlaps', 'FractionOverlaps'])
    #     self.sample_to_pyrange_map = TruthCallset._construct_sample_to_pyrange_map(self.truth_callset_pyrange, self.sample_set)

    def subset_intervals_to_rare_regions(self, intervals: IntervalCollection, max_allelic_fraction: float) -> IntervalCollection:
        """
        Subset a given interval collection to those intervals who only overlap with rare events, as defined by
        max_allele_fraction

        :param intervals: interval collection to subset
        :param max_allelic_fraction: events below that threshold are considered rare
        :return: subset of rare intervals
        """
        # TODO Handle the case of an interval from collection overlapping multiple events. In that case we should only consider
        # TODO allelic fraction of the event with the largest reciprocal overlap
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

    def get_name(self):
        return "Truth_callset"


class GCNVCallset(Callset):

    def __init__(self, sample_to_pyrange_map: dict, joint_callset: pr.PyRanges, interval_collection: IntervalCollection):
        super().__init__(sample_to_pyrange_map, joint_callset, interval_collection)

    @classmethod
    def read_in_callset(cls, gcnv_segment_vcfs: List[str], gcnv_callset_tsv: Optional[str], gcnv_joint_vcf: Optional[str],
                        interval_collection: IntervalCollection, max_events_allowed: int = 100):

        sample_to_pyrange_map = {}

        if gcnv_callset_tsv:
            callset_df = pd.read_csv(open(gcnv_callset_tsv), sep='\t')

            def _get_attribute_list(call, num_exon, qs, sf, hg38_str):
                split_str = hg38_str.split(":")
                contig = split_str[0]
                start = split_str[1].split("-")[0]
                end = split_str[1].split("-")[1]
                event_type = EventType.get_event_type_from_svtype(call)
                return [contig, start, end, event_type, num_exon, qs, sf]

            sample_to_event_tuples = [(sample, _get_attribute_list(call, num_exon, qs, sf, hg38_str)) for
                                      sample, call, num_exon, qs, sf, hg38_str in
                                      zip(callset_df['sample'], callset_df['call'], callset_df['num_exon'],
                                          callset_df['QS'], callset_df['sf'], callset_df['hg38']) if hg38_str != "NONE" and num_exon != 0]

            sample_to_event_list_map = {}
            for sample, event in sample_to_event_tuples:
                sample_to_event_list_map.setdefault(sample, []).append(event)

            for sample in sample_to_event_list_map.keys():
                events_df = pd.DataFrame(sample_to_event_list_map[sample], columns=Callset.CALLSET_COLUMNS)
                sample_to_pyrange_map[sample] = pr.PyRanges(events_df)
        else:
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

            # Read in single sample vcfs

            for vcf_file in gcnv_segment_vcfs:
                vcf_reader = vcf.Reader(open(vcf_file, 'r'))
                assert len(vcf_reader.samples) == 1
                sample_name = vcf_reader.samples[0]
                events_df_lists = []
                for record in vcf_reader:
                    #event_type = EventType.gcnv_call_to_event_type(record.genotype(sample_name)['GT'])
                    cn = record.genotype(sample_name)['CN']
                    event_type = EventType.REF if cn == 2 else (EventType.DUP if cn > 2 else EventType.DEL)
                    if event_type == EventType.REF:
                        continue
                    events_df_lists.append([record.CHROM, int(record.POS), int(record.INFO['END']), event_type,
                             int(record.genotype(sample_name)['NP']), int(record.genotype(sample_name)['QS']), 0.0])

                events_df = pd.DataFrame(events_df_lists, columns=Callset.CALLSET_COLUMNS)
                if len(events_df) > max_events_allowed:
                   continue
                events_pr = pr.PyRanges(events_df)
                sample_to_pyrange_map[sample_name] = events_pr

        # Annotate callset with AF: clustering method
        # callset_df_list = [p.df.assign(Sample=s) for s, p in sample_to_pyrange_map.items()]
        # joined_callset_pr = pr.PyRanges(pd.concat(callset_df_list))
        # clustered_joined_callset = joined_callset_pr.cluster()
        # mapping = clustered_joined_callset.df.Cluster.value_counts()
        # sample_set = set(sample_to_pyrange_map.keys())
        # num_samples = len(sample_to_pyrange_map)
        # clustered_joined_callset = pr.PyRanges(clustered_joined_callset.df.assign(
        #     Frequency=clustered_joined_callset.df.Cluster.apply(lambda x: mapping[x] / num_samples)))
        # for s in sample_set:
        #     sample_to_pyrange_map[s] = pr.PyRanges(
        #         clustered_joined_callset.df.loc[clustered_joined_callset.df['Sample'] == s].drop(columns=["Sample", "Cluster"]))

            # Annotate callset with AF: calculate overlaps method
            overlaps_pr = pr.count_overlaps(sample_to_pyrange_map)
            print(overlaps_pr.df)
            num_samples = len(sample_to_pyrange_map)
            sample_set = set(sample_to_pyrange_map.keys())
            for sample in sample_set:
                event_frequencies = []
                for index, event in sample_to_pyrange_map[sample].df.iterrows():
                    overlapping_segments_pr = overlaps_pr[event["Chromosome"], event["Start"]:event["End"]]
                    sample_to_length_map = defaultdict(lambda: 0)
                    for i, e in overlapping_segments_pr.df.iterrows():
                        overlap_length = e["End"] - e["Start"]
                        for s in sample_set:
                            if e[s]:
                                sample_to_length_map[s] += overlap_length
                    length = event["End"] - event["Start"]
                    number_overlaps = sum(sample_to_length_map[s] / length > 0.5 for s in sample_set)
                    af = number_overlaps / num_samples
                    event_frequencies.append(af)
                    #print(sample_to_pyrange_map[sample].Frequency[index])
                sample_to_pyrange_map[sample].Frequency = pd.Series(event_frequencies)
            print("Done calculating AF")

        # Read in joint vcf
        joint_callset_pr = None
        if gcnv_joint_vcf is not None:
            joint_vcf_reader = vcf.Reader(open(gcnv_joint_vcf, 'r'))
            joint_sample_list = list(joint_vcf_reader.samples)
            events_list = []
            for record in joint_vcf_reader:
                if record.FILTER != []:
                    continue
                del_call_samples = []
                dup_call_samples = []
                if record.ALT[0] is None:
                    continue
                for sample in joint_sample_list:
                    # TODO this does not work properly for allosomal contigs but we ignore them for now
                    if record.genotype(sample)['CN']:  # check that there is a value stored in the field
                        if record.genotype(sample)['CN'] < 2:
                            if record.genotype(sample)['QS'] >= 100:
                                del_call_samples.append(sample)
                        elif record.genotype(sample)['CN'] > 2:
                            if record.genotype(sample)['QS'] >= 50:
                                dup_call_samples.append(sample)
                chrom, start, end, num_bins = record.CHROM, int(record.POS), int(record.INFO['END']), int(
                    record.genotype(joint_sample_list[0])['NP'])
                name_prefix = str(chrom) + "_" + str(start) + "_" + str(end) + "_"
                if del_call_samples:
                    events_list.append([chrom, start, end, name_prefix + "DEL", EventType.DEL, frozenset(del_call_samples), 0.0, num_bins])
                if dup_call_samples:
                    events_list.append([chrom, start, end, name_prefix + "DUP", EventType.DUP, frozenset(dup_call_samples), 0.0, num_bins])
            joint_events_df = pd.DataFrame(events_list, columns=Callset.JOINT_CALLSET_COLUMNS)
            joint_events_df.astype(Callset.JOINT_CALLSET_COLUMN_TYPES)
            joint_callset_pr = pr.PyRanges(joint_events_df)

            #sample_to_pyrange_map = {s: None for s in joint_sample_list}

        #sample_to_pyrange_map = Callset._construct_sample_to_pyrange_map(joint_callset_pr, frozenset(joint_sample_list))

        return cls(sample_to_pyrange_map, joint_callset_pr, interval_collection)

    def get_name(self):
        return "gCNV_callset"


class XHMMCallset(Callset):
    def __init__(self, sample_to_pyrange_map: dict, joint_callset: pr.PyRanges, interval_collection: IntervalCollection):
        super().__init__(sample_to_pyrange_map, joint_callset,  interval_collection)

    @classmethod
    def read_in_callset(cls, xhmm_vcfs: List[str], interval_collection: IntervalCollection, samples_to_keep: set = None):
        assert len(xhmm_vcfs) > 0

        events_list = []
        joint_sample_list = []
        sample_to_df_lists = {}
        for xhmm_vcf in xhmm_vcfs:
            vcf_reader = vcf.Reader(open(xhmm_vcf, 'r'))
            sample_list = vcf_reader.samples
            joint_sample_list.extend(sample_list)
            if samples_to_keep:
                list(set(sample_list).intersection(samples_to_keep))
            sample_to_df_lists.update({sample: [] for sample in sample_list})  # TODO refactor

            for record in vcf_reader:
                del_call_samples = []
                dup_call_samples = []

                chrom, start, end, num_bins = record.CHROM, int(record.POS), int(record.INFO['END']), int(record.INFO['NUMT'])
                af = sum([float(af) for af in record.INFO['AF']])

                for sample in sample_list:
                    event_type = EventType.gcnv_call_to_event_type(record.genotype(sample)['GT'])
                    if event_type == EventType.DEL:
                        del_call_samples.append(sample)
                        sample_to_df_lists[sample].append([chrom, start, end, event_type, num_bins, int(record.genotype(sample)['NDQ']), af]) # TODO refactor
                    if event_type == EventType.DUP:
                        dup_call_samples.append(sample)
                        sample_to_df_lists[sample].append([chrom, start, end, event_type, num_bins, int(record.genotype(sample)['NDQ']), af])  # TODO refactor

                name_prefix = str(chrom) + "_" + str(start) + "_" + str(end) + "_"
                af = (len(del_call_samples) + len(dup_call_samples)) / len(sample_list)  # TODO refactor

                if del_call_samples:
                    events_list.append([chrom, start, end, name_prefix + "DEL", EventType.DEL, frozenset(del_call_samples), af, num_bins])
                if dup_call_samples:
                    events_list.append([chrom, start, end, name_prefix + "DUP", EventType.DUP, frozenset(dup_call_samples), af, num_bins])

        events_df = pd.DataFrame(events_list, columns=Callset.JOINT_CALLSET_COLUMNS)
        events_df.astype(Callset.JOINT_CALLSET_COLUMN_TYPES)
        joint_callset_pr = pr.PyRanges(events_df)

        sample_to_pyrange_map = {}  # TODO refactor
        for sample in joint_sample_list:
            sample_df = pd.DataFrame(sample_to_df_lists[sample], columns=Callset.CALLSET_COLUMNS)
            sample_df.astype(Callset.CALLSET_COLUMN_TYPES)
            sample_to_pyrange_map[sample] = pr.PyRanges(sample_df) # TODO refactor
        #sample_to_pyrange_map = Callset._construct_sample_to_pyrange_map(joint_callset_pr, frozenset(sample_list))
        return cls(sample_to_pyrange_map, joint_callset_pr, interval_collection)

    def get_name(self):
        return "XHMM_callset"
