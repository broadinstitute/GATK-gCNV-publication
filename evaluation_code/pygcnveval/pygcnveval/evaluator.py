from callset import TruthCallset, GCNVCallset
from evaluation_result import EvaluationResult
import random


class PerEventEvaluator:

    def __init__(self, truth_callset: TruthCallset):
        self.truth_callset = truth_callset

    def evaluate_callset_against_the_truth(self, gcnv_callset: GCNVCallset, minimum_reciprocal_overlap: float = 0.3)\
            -> EvaluationResult:
        evaluation_result = EvaluationResult()

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
