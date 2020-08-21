import numpy as np

from event import EventType
from callset import TruthCallset, GCNVCallset, Callset
from evaluation_result import PerEventEvaluationResult, PerBinEvaluationResult
from interval_collection import IntervalCollection


class PerEventEvaluator:

    def __init__(self, truth_callset: TruthCallset, callset: Callset):
        self.truth_callset = truth_callset
        self.callset = callset
        self.sample_list_to_eval = list(self.truth_callset.sample_set & self.callset.sample_set)

    def evaluate_callset_against_the_truth(self, minimum_overlap: float = 0.2,
                                           min_quality_threshold: int = 30) -> PerEventEvaluationResult:
        evaluation_result = PerEventEvaluationResult(self.callset.get_name())

        # Calculate precision
        for validated_event in self.callset.get_event_generator(self.sample_list_to_eval, min_quality_threshold):
            overlapping_truth_events = self.truth_callset.get_overlapping_events_for_sample(validated_event.interval,
                                                                                            validated_event.sample)
            overlapping_truth_event_best_match = validated_event.find_event_with_largest_overlap(overlapping_truth_events)
            if overlapping_truth_event_best_match:
                if overlapping_truth_event_best_match.call_attributes['Frequency'] > 0.01:
                    continue
                event_validates = validated_event.compare_to(overlapping_truth_event_best_match, minimum_overlap)
                evaluation_result.update_precision(validated_event.call_attributes['NumBins'], event_validates)
            else:
                evaluation_result.update_precision(validated_event.call_attributes['NumBins'], False)

        # Calculate recall
        for truth_event in self.truth_callset.get_event_generator(self.sample_list_to_eval):
            if truth_event.call_attributes['Frequency'] > 0.01:
                continue
            overlapping_gcnv_events = self.callset.get_overlapping_events_for_sample(truth_event.interval, truth_event.sample)
            overlapping_gcnv_events = [e for e in overlapping_gcnv_events if e.call_attributes['Quality'] >= min_quality_threshold]
            if not overlapping_gcnv_events:
                evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False)
            else:
                overlapping_gcnv_event_best_match = truth_event.find_event_with_largest_overlap(overlapping_gcnv_events)
                if overlapping_gcnv_event_best_match.event_type == truth_event.event_type:
                    evaluation_result.update_recall(truth_event.call_attributes['NumBins'], True)
                else:
                    evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False)
            overlapping_gcnv_event_best_match = truth_event.find_event_with_largest_overlap(overlapping_gcnv_events)
            # if overlapping_gcnv_event_best_match:
            #     event_validates = truth_event.compare_to(overlapping_gcnv_event_best_match, minimum_overlap)
            #     evaluation_result.update_recall(truth_event.call_attributes['NumBins'], event_validates)
            # else:
            #     evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False)

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
