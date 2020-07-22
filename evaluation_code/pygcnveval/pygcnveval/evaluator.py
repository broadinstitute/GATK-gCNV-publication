import numpy as np

from event import EventType
from callset import TruthCallset, GCNVCallset
from evaluation_result import PerEventEvaluationResult, PerBinEvaluationResult
from interval_collection import IntervalCollection


class PerEventEvaluator:

    def __init__(self, truth_callset: TruthCallset, gcnv_callset: GCNVCallset):
        self.truth_callset = truth_callset
        self.gcnv_callset = gcnv_callset
        self.sample_list_to_eval = list(self.truth_callset.sample_set & self.gcnv_callset.sample_set)

    def evaluate_callset_against_the_truth(self, gcnv_callset: GCNVCallset, regions_to_ignore: IntervalCollection = None,
                                           minimum_reciprocal_overlap: float = 0.3) -> PerEventEvaluationResult:
        evaluation_result = PerEventEvaluationResult()

        # Calculate precision
        for gcnv_event in gcnv_callset.get_event_generator(self.sample_list_to_eval, 60):
            if regions_to_ignore and \
                    len(regions_to_ignore.pyrange[gcnv_event.interval.chrom, gcnv_event.interval.start:gcnv_event.interval.end]) > 0:
                continue  # skip common regions
            overlapping_truth_events = self.truth_callset.get_overlapping_events_for_sample(gcnv_event.interval,
                                                                                            gcnv_event.sample)
            overlapping_truth_event_best_match = gcnv_event.find_event_with_largest_reciprocal_overlap(overlapping_truth_events)
            if overlapping_truth_event_best_match:
                event_validates = gcnv_event.compare_to(overlapping_truth_event_best_match, minimum_reciprocal_overlap)
                evaluation_result.update_precision(gcnv_event.call_attributes['NumBins'], event_validates)
            else:
                evaluation_result.update_precision(gcnv_event.call_attributes['NumBins'], False)

        # Calculate recall
        for truth_event in self.truth_callset.get_event_generator(self.sample_list_to_eval):
            if regions_to_ignore and \
                    len(regions_to_ignore.pyrange[truth_event.interval.chrom, truth_event.interval.start:truth_event.interval.end]) > 0:
                continue  # skip common regions
            overlapping_gcnv_events = gcnv_callset.get_overlapping_events_for_sample(truth_event.interval, truth_event.sample)
            overlapping_gcnv_events = [e for e in overlapping_gcnv_events if e.call_attributes['Quality'] > 60]
            overlapping_gcnv_event_best_match = truth_event.find_event_with_largest_reciprocal_overlap(overlapping_gcnv_events)
            if overlapping_gcnv_event_best_match:
                event_validates = truth_event.compare_to(overlapping_gcnv_event_best_match, minimum_reciprocal_overlap)
                evaluation_result.update_recall(truth_event.call_attributes['NumBins'], event_validates)
            else:
                evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False)

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

        for threshold in evaluation_result.quality_thresholds:
            true_positives = np.sum(np.logical_and(gcnv_quality_matrix >= threshold,
                                                   np.logical_and(truth_callset_matrix == gcnv_callset_matrix,
                                                                  gcnv_callset_matrix != EventType.REF.value)))
            false_positives = np.sum(np.logical_and(gcnv_quality_matrix >= threshold,
                                                    np.logical_and(truth_callset_matrix != gcnv_callset_matrix,
                                                                   gcnv_callset_matrix != EventType.REF.value)))

            positives = np.sum(truth_callset_matrix != EventType.REF.value)

            # We want to not double count false positives where both calls are non-ref, but disagree as false negatives
            false_negatives_adjust = np.sum(np.logical_and(gcnv_quality_matrix >= threshold,
                                                           np.logical_and(truth_callset_matrix != gcnv_callset_matrix,
                                                                          np.logical_and(truth_callset_matrix != EventType.REF.value,
                                                                                         gcnv_callset_matrix != EventType.REF.value))))

            false_negatives = (positives - true_positives) - false_negatives_adjust
            evaluation_result.set_values_for_quality_cutoff(threshold, true_positives, false_positives, false_negatives)

        return evaluation_result
