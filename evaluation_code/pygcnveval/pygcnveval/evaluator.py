import numpy as np
from typing import List

from event import EventType
from callset import TruthCallset, GCNVCallset, Callset
from evaluation_result import PerEventEvaluationResult, PerBinEvaluationResult
from interval_collection import IntervalCollection


class PerEventEvaluator:

    def __init__(self, truth_callset: TruthCallset, callset: Callset, sample_list_to_evaluate: List):
        self.truth_callset = truth_callset
        self.callset = callset
        assert set(sample_list_to_evaluate).issubset(self.truth_callset.sample_set)
        assert set(sample_list_to_evaluate).issubset(self.callset.sample_set)
        self.sample_list_to_eval = sample_list_to_evaluate

    def evaluate_callset_against_the_truth(self, minimum_overlap: float = 0.2,
                                           gcnv_sq_min_del: int = 100, gcnv_sq_min_dup: int = 50) -> PerEventEvaluationResult:
        evaluation_result = PerEventEvaluationResult(self.callset.get_name())

        # construct callset filter
        def gcnv_calls_filter(e):
            is_allosome = e.interval.chrom == "chrX" or e.interval.chrom == "chrY" \
                          or e.interval.chrom == "X" or e.interval.chrom == "Y"
            return not is_allosome and (e.call_attributes['Frequency'] <= 0.02) \
                   and ((e.call_attributes['Quality'] >= gcnv_sq_min_del and e.event_type == EventType.DEL)
                        or (e.call_attributes['Quality'] >= gcnv_sq_min_dup and e.event_type == EventType.DUP))

        def xhmm_calls_filter(e):
            is_allosome = e.interval.chrom == "chrX" or e.interval.chrom == "chrY" \
                          or e.interval.chrom == "X" or e.interval.chrom == "Y"
            return not is_allosome and (e.call_attributes['Frequency'] <= 0.02) and (e.call_attributes['Quality'] >= 60)

        calls_filter = None
        if self.callset.get_name() == "gCNV_callset":
            calls_filter = gcnv_calls_filter
        elif self.callset.get_name() == "XHMM_callset":
            calls_filter = xhmm_calls_filter

        # Calculate precision
        for validated_event in self.callset.get_event_generator(self.sample_list_to_eval, calls_filter):

            # if self.callset.get_name() == "gCNV_callset" and validated_event.call_attributes['Quality'] < 50 and validated_event.event_type == EventType.DUP:
            #     continue
            # if self.callset.get_name() == "gCNV_callset" and validated_event.call_attributes['Quality'] < 100 and validated_event.event_type == EventType.DEL:
            #     continue

            # if self.callset.get_name() == "XHMM_callset" and validated_event.call_attributes['Quality'] < 60:
            #     continue
            #
            # if validated_event.interval.chrom == "chrX" or validated_event.interval.chrom == "chrY" or validated_event.interval.chrom == "X" or validated_event.interval.chrom == "Y":
            #     continue
            #
            # if self.callset.get_name() == "XHMM_callset" and validated_event.call_attributes['Frequency'] > 0.02:
            #     continue
            # # check if event is common
            # if self.callset.get_name() == "gCNV_callset" and validated_event.call_attributes['Frequency'] > 0.05:
            #     continue

            overlapping_truth_events = self.truth_callset.get_overlapping_events_for_sample(validated_event.interval,
                                                                                            validated_event.sample)
            overlapping_truth_event_best_match = validated_event.find_event_with_largest_overlap(overlapping_truth_events)

            if overlapping_truth_event_best_match:
                if overlapping_truth_event_best_match.call_attributes['Frequency'] > 0.01:
                    evaluation_result.update_precision(validated_event.call_attributes['NumBins'], None,
                                                       validated_event.event_type, validated_event.interval,
                                                       "FILTERED_HIGH_TRUTH_AF", list([validated_event.sample]))
                    continue
                event_validates, rejection_reason = validated_event.compare_to(overlapping_truth_event_best_match, minimum_overlap)
                evaluation_result.update_precision(validated_event.call_attributes['NumBins'], event_validates,
                                                   validated_event.event_type, validated_event.interval,
                                                   rejection_reason, list([validated_event.sample]))
            else:
                evaluation_result.update_precision(validated_event.call_attributes['NumBins'], False,
                                                   validated_event.event_type, validated_event.interval,
                                                   "NO_OVERLAPPING_EVENT_FOUND", list([validated_event.sample]))

        # Calculate recall
        for truth_event in self.truth_callset.get_event_generator(self.sample_list_to_eval):
            if truth_event.interval.chrom == "chrX" or truth_event.interval.chrom == "chrY" or truth_event.interval.chrom == "X" or truth_event.interval.chrom == "Y":
                continue
            if truth_event.call_attributes['Frequency'] > 0.01:
                continue
            overlapping_gcnv_events = self.callset.get_overlapping_events_for_sample(truth_event.interval, truth_event.sample)

            if not overlapping_gcnv_events:
                evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False,
                                                truth_event.event_type, truth_event.interval,
                                                "NO_OVERLAPPING_EVENT_FOUND", list([truth_event.sample]))
            else:
                overlapping_gcnv_event_best_match = truth_event.find_event_with_largest_overlap(overlapping_gcnv_events)
                if self.callset.get_name() == "gCNV_callset" and overlapping_gcnv_event_best_match.call_attributes[
                        'Quality'] < gcnv_sq_min_dup and overlapping_gcnv_event_best_match.event_type == EventType.DUP:
                    evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False,
                                                    truth_event.event_type, truth_event.interval,
                                                    "FILTERED_LOW_QUAL", list([truth_event.sample]))
                    continue
                if self.callset.get_name() == "gCNV_callset" and overlapping_gcnv_event_best_match.call_attributes[
                        'Quality'] < gcnv_sq_min_del and overlapping_gcnv_event_best_match.event_type == EventType.DEL:
                    evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False,
                                                    truth_event.event_type, truth_event.interval,
                                                    "FILTERED_LOW_QUAL", list([truth_event.sample]))
                    continue

                if self.callset.get_name() == "XHMM_callset" and overlapping_gcnv_event_best_match.call_attributes['Frequency'] > 0.02:
                    evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False,
                                                    truth_event.event_type, truth_event.interval,
                                                    "FILTERED_HIGH_AF", list([truth_event.sample])) # TODO refactor
                    continue
                if self.callset.get_name() == "XHMM_callset" and overlapping_gcnv_event_best_match.call_attributes['Quality'] < 60:
                    evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False,
                                                    truth_event.event_type, truth_event.interval,
                                                    "FILTERED_LOW_QUAL", list([truth_event.sample]))  # TODO refactor
                    continue

                if self.callset.get_name() == "gCNV_callset" and overlapping_gcnv_event_best_match.call_attributes['Frequency'] > 0.02:
                    evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False, truth_event.event_type,
                                                    truth_event.interval, "FILTERED_HIGH_AF", list([truth_event.sample]))
                    continue
                event_validates, rejection_reason = truth_event.compare_to(overlapping_gcnv_event_best_match, minimum_overlap)
                evaluation_result.update_recall(truth_event.call_attributes['NumBins'], event_validates,
                                                truth_event.event_type, truth_event.interval,
                                                rejection_reason, list([truth_event.sample]))

        return evaluation_result


class PerSiteEvaluator:

    def __init__(self, truth_callset: TruthCallset, callset: Callset):
        self.truth_callset = truth_callset
        self.callset = callset
        self.sample_list_to_eval = list(self.truth_callset.sample_set & self.callset.sample_set)

    def evaluate_callset_against_the_truth(self) -> PerEventEvaluationResult:
        evaluation_result = PerEventEvaluationResult(self.callset.get_name())
        overall_events = 0
        # Calculate precision
        for validated_allele in self.callset.get_allele_generator():

            if len(validated_allele.sample_set.intersection(set(self.sample_list_to_eval))) == 0:
                continue
            if validated_allele.interval.chrom == "chrX" or validated_allele.interval.chrom == "chrY":
                continue
            overall_events += len(validated_allele.sample_set)
            overlapping_truth_alleles = self.truth_callset.get_overlapping_alleles(validated_allele.interval)
            overlapping_truth_allele_best_match = validated_allele.find_allele_with_largest_overlap(
                overlapping_truth_alleles)

            if overlapping_truth_allele_best_match:
                if overlapping_truth_allele_best_match.allele_attributes['Frequency'] > 0.01:
                    evaluation_result.update_precision(validated_allele.allele_attributes['NumBins'], None,
                                                       validated_allele.event_type, validated_allele.interval,
                                                       "FILTERED_HIGH_TRUTH_AF", list(validated_allele.sample_set))
                    continue
                allele_validates, rejection_reason = validated_allele.compare_to(overlapping_truth_allele_best_match, 0.2, 0.2, set(self.sample_list_to_eval))
                evaluation_result.update_precision(validated_allele.allele_attributes['NumBins'], allele_validates,
                                                   validated_allele.event_type, validated_allele.interval,
                                                   rejection_reason, list(validated_allele.sample_set))
            else:
                evaluation_result.update_precision(validated_allele.allele_attributes['NumBins'], False,
                                                   validated_allele.event_type, validated_allele.interval,
                                                   "NO_OVERLAPPING_EVENT_FOUND", list(validated_allele.sample_set))

        # Calculate recall
        for truth_allele in self.truth_callset.get_allele_generator():
            if len(truth_allele.sample_set.intersection(set(self.sample_list_to_eval))) == 0:
                continue
            if truth_allele.interval.chrom == "chrX" or truth_allele.interval.chrom == "chrY":
                continue
            if truth_allele.allele_attributes['Frequency'] > 0.01:
                continue

            overlapping_gcnv_alleles = self.callset.get_overlapping_alleles(truth_allele.interval)

            if not overlapping_gcnv_alleles:
                evaluation_result.update_recall(truth_allele.allele_attributes['NumBins'], False,
                                                truth_allele.event_type, truth_allele.interval,
                                                "NO_OVERLAPPING_EVENT_FOUND", list(truth_allele.sample_set))
            else:
                overlapping_gcnv_allele_best_match = truth_allele.find_allele_with_largest_overlap(overlapping_gcnv_alleles)
                event_validates, rejection_reason = truth_allele.compare_to(overlapping_gcnv_allele_best_match, 0.2, 0.2, set(self.sample_list_to_eval))
                evaluation_result.update_recall(truth_allele.allele_attributes['NumBins'], event_validates,
                                                truth_allele.event_type, truth_allele.interval,
                                                rejection_reason, list(truth_allele.sample_set))

        print(overall_events)
        return evaluation_result


class PerBinEvaluator:

    def __init__(self, truth_callset: TruthCallset, gcnv_callset: GCNVCallset, interval_collection: IntervalCollection):
        self.truth_callset = truth_callset
        self.gcnv_callset = gcnv_callset
        self.sample_list_to_eval = list(self.truth_callset.sample_set & self.gcnv_callset.sample_set)
        self.interval_collection = interval_collection

    def evaluate_callset_against_the_truth(self, gcnv_callset: GCNVCallset, number_quality_thresholds: int = 30):
        print("Creating interval by sample genotype matrix for truth callset...", flush=True)
        truth_matrix_view = self.truth_callset.get_callset_matrix_view(self.interval_collection, self.sample_list_to_eval)
        print("Creating interval by sample genotype matrix for gCNV callset...", flush=True)
        gcnv_matrix_view = gcnv_callset.get_callset_matrix_view(self.interval_collection, self.sample_list_to_eval)
        truth_callset_matrix = truth_matrix_view.samples_by_intervals_call_matrix
        gcnv_callset_matrix = gcnv_matrix_view.samples_by_intervals_call_matrix
        gcnv_quality_matrix = gcnv_matrix_view.samples_by_intervals_quality_matrix
        evaluation_result = PerBinEvaluationResult(gcnv_matrix_view.samples_by_intervals_quality_matrix,
                                                   number_quality_thresholds)
        print(evaluation_result.quality_thresholds)
        for threshold in evaluation_result.quality_thresholds:
            true_positives = np.sum(np.logical_and(gcnv_quality_matrix >= threshold,
                                                   np.logical_and(truth_callset_matrix == gcnv_callset_matrix,
                                                                  gcnv_callset_matrix != EventType.REF.value)))
            false_positives = np.sum(np.logical_and(gcnv_quality_matrix >= threshold,
                                                    np.logical_and(truth_callset_matrix != gcnv_callset_matrix,
                                                                   gcnv_callset_matrix != EventType.REF.value)))

            positives = np.sum(truth_callset_matrix != EventType.REF.value)

            # We do not want to double count false positives where both calls are non-ref, but disagree as false negatives
            false_negatives_adjust = np.sum(np.logical_and(gcnv_quality_matrix >= threshold,
                                                           np.logical_and(truth_callset_matrix != gcnv_callset_matrix,
                                                                          np.logical_and(truth_callset_matrix != EventType.REF.value,
                                                                                         gcnv_callset_matrix != EventType.REF.value))))

            false_negatives = (positives - true_positives) - false_negatives_adjust
            evaluation_result.set_values_for_quality_cutoff(threshold, true_positives, false_positives, false_negatives)

        return evaluation_result
