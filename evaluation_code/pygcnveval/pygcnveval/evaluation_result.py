import numpy as np


class PerEventEvaluationResult:
    DEFAULT_EVENT_SIZE_BINS = list(range(25))

    def __init__(self, tool_name: str, event_size_bins=DEFAULT_EVENT_SIZE_BINS):
        self.event_size_bins = event_size_bins
        self.tool_name = tool_name
        self.precision_size_to_tp = {i: 0 for i in event_size_bins}
        self.precision_size_to_fp = {i: 0 for i in event_size_bins}

        self.recall_size_to_tp = {i: 0 for i in event_size_bins}
        self.recall_size_to_fn = {i: 0 for i in event_size_bins}

    def update_precision(self, num_intervals: int, event_valid: bool):
        index = min(num_intervals - 1, len(self.event_size_bins) - 1)
        if event_valid:
            self.precision_size_to_tp[index] += 1
        else:
            self.precision_size_to_fp[index] += 1

    def update_recall(self, num_intervals: int, event_valid: bool):
        index = min(num_intervals - 1, len(self.event_size_bins) - 1)
        if event_valid:
            self.recall_size_to_tp[index] += 1
        else:
            self.recall_size_to_fn[index] += 1

    def __eq__(self, other):
        if isinstance(other, PerEventEvaluationResult):
            return self.tool_name == other.tool_name \
                   and self.precision_size_to_tp == other.precision_size_to_tp \
                   and self.precision_size_to_fp == other.precision_size_to_fp \
                   and self.recall_size_to_tp == other.recall_size_to_tp \
                   and self.recall_size_to_fn == other.recall_size_to_fn
        return False


class PerBinEvaluationResult:

    def __init__(self, quality_matrix: np.ndarray, number_points):
        """
        TODO
        :param quality_matrix:
        :param number_points:
        """
        self.quality_thresholds = [0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000]
        #self.quality_thresholds = sorted(list(set(np.percentile(
        #    quality_matrix, [(i / number_points) * 100 for i in range(number_points + 1)]))))
        self.true_positives = {q: 0 for q in self.quality_thresholds}
        self.false_positives = {q: 0 for q in self.quality_thresholds}
        self.false_negatives = {q: 0 for q in self.quality_thresholds}

    def set_values_for_quality_cutoff(self, threshold, true_positives, false_positives, false_negatives):
        self.true_positives[threshold] = true_positives
        self.false_positives[threshold] = false_positives
        self.false_negatives[threshold] = false_negatives

    def get_precision_array(self):
        return np.array([self.true_positives[q] / (self.true_positives[q] + self.false_positives[q])
                         if self.true_positives[q] != 0. else 0
                         for q in self.quality_thresholds])

    def get_recall_array(self):
        return np.array([self.true_positives[q] / (self.true_positives[q] + self.false_negatives[q])
                         if self.true_positives[q] != 0. else 0
                         for q in self.quality_thresholds])

    def __eq__(self, other):
        if isinstance(other, PerBinEvaluationResult):
            return self.true_positives == other.true_positives \
                   and self.false_positives == other.false_positives \
                   and self.false_negatives == other.false_negatives
        return False
