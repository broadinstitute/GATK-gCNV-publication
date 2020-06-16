import random
import numpy as np

from event import EventType
from callset import TruthCallset, GCNVCallset
from evaluation_result import PerEventEvaluationResult, PerBinEvaluationResult
from interval_collection import IntervalCollection


class PerEventEvaluator:

    def __init__(self, truth_callset: TruthCallset):
        self.truth_callset = truth_callset

    def evaluate_callset_against_the_truth(self, gcnv_callset: GCNVCallset, regions_to_ignore: IntervalCollection = None,
                                           minimum_reciprocal_overlap: float = 0.3) -> PerEventEvaluationResult:
        evaluation_result = PerEventEvaluationResult()
        # TODO skip any region intersecting with intervals in `regions_to_ignore`

        # Calculate precision
        for gcnv_event in gcnv_callset.get_event_generator():
            overlapping_truth_events = self.truth_callset.get_overlapping_events_for_sample(gcnv_event.interval, gcnv_event.sample)
            if overlapping_truth_events:
                event_validates = gcnv_event.compare_to(random.choice(overlapping_truth_events), minimum_reciprocal_overlap) # TODO replace random with something more sensible
                evaluation_result.update_precision(gcnv_event.call_attributes['NumBins'], event_validates)
            else:
                evaluation_result.update_precision(gcnv_event.call_attributes['NumBins'], False)

        # Calculate recall
        for truth_event in self.truth_callset.get_event_generator():
            overlapping_gcnv_events = gcnv_callset.get_overlapping_events_for_sample(truth_event.interval, truth_event.sample)
            if overlapping_gcnv_events:
                event_validates = truth_event.compare_to(random.choice(overlapping_gcnv_events), minimum_reciprocal_overlap) # TODO replace random with something more sensible
                evaluation_result.update_recall(truth_event.call_attributes['NumBins'], event_validates)
            else:
                evaluation_result.update_recall(truth_event.call_attributes['NumBins'], False)

        return evaluation_result


class PerBinEvaluator:

    def __init__(self, truth_callset: TruthCallset, interval_collection: IntervalCollection):
        self.truth_callset = truth_callset
        self.interval_collection = interval_collection

    def evaluate_callset_against_the_truth(self, gcnv_callset: GCNVCallset, number_quality_thresholds: int = 30):
        truth_matrix_view = self.truth_callset.get_callset_matrix_view(self.interval_collection, list(gcnv_callset.sample_set))
        gcnv_matrix_view = gcnv_callset.get_callset_matrix_view(self.interval_collection, list(gcnv_callset.sample_set))
        truth_callset_matrix = truth_matrix_view.intervals_by_sample_call_matrix
        gcnv_callset_matrix = gcnv_matrix_view.intervals_by_sample_call_matrix
        gcnv_quality_matrix = gcnv_matrix_view.intervals_by_sample_quality_matrix
        evaluation_result = PerBinEvaluationResult(gcnv_matrix_view.intervals_by_sample_quality_matrix,
                                                   number_quality_thresholds)

        for threshold in evaluation_result.quality_thresholds:
            true_positives = np.sum(np.logical_and(gcnv_quality_matrix <= threshold,
                                                   np.logical_and(truth_callset_matrix == gcnv_callset_matrix,
                                                                  gcnv_callset_matrix != EventType.REF.value)))
            false_positives = np.sum(np.logical_and(gcnv_quality_matrix <= threshold,
                                                    np.logical_and(truth_callset_matrix != gcnv_callset_matrix,
                                                                   gcnv_callset_matrix != EventType.REF.value)))
            false_negatives = np.sum(np.logical_and(gcnv_quality_matrix <= threshold,
                                                    np.logical_and(truth_callset_matrix != gcnv_callset_matrix,
                                                                   gcnv_callset_matrix == EventType.REF.value)))
            evaluation_result.set_values_for_quality_cutoff(threshold, true_positives, false_positives, false_negatives)

        return evaluation_result
