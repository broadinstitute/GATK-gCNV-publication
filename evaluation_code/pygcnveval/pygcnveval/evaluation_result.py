import numpy as np
import pandas as pd
import os
from typing import List, Optional

from interval import Interval
from event import EventType


class PerEventEvaluationResult:
    DEFAULT_EVENT_SIZE_BINS = list(range(25))
    RESULT_COLUMNS = ["Chromosome", "Start", "End", "Genotype", "NumBins", "ClassificationStatus", "ClassificationReason", "Samples"]

    def __init__(self, tool_name: str, event_size_bins=DEFAULT_EVENT_SIZE_BINS):
        self.event_size_bins = event_size_bins
        self.tool_name = tool_name
        self.precision_results_df = pd.DataFrame(columns=PerEventEvaluationResult.RESULT_COLUMNS)
        self.recall_results_df = pd.DataFrame(columns=PerEventEvaluationResult.RESULT_COLUMNS)

        self.precision_size_to_tp = {i: 0 for i in event_size_bins}
        self.precision_size_to_fp = {i: 0 for i in event_size_bins}

        self.recall_size_to_tp = {i: 0 for i in event_size_bins}
        self.recall_size_to_fn = {i: 0 for i in event_size_bins}

    def update_precision(self, num_bins: int, event_valid: Optional[bool], event_type: EventType, interval: Interval, reason: str, sample_list: List[str]):
        classification_status = "NA" if event_valid is None else ("TP" if event_valid else "FP")
        classification_reason = reason if reason is not None else "NA"
        self.precision_results_df = self.precision_results_df.append({"Chromosome": interval.chrom,
                                                                      "Start": interval.start,
                                                                      "End": interval.end,
                                                                      "Genotype": event_type.name,
                                                                      "NumBins": num_bins,
                                                                      "ClassificationStatus": classification_status,
                                                                      "ClassificationReason": classification_reason,
                                                                      "Samples": ",".join(sample_list)},
                                                                     ignore_index=True)
        if event_valid is not None:
            index = min(num_bins - 1, len(self.event_size_bins) - 1)
            if event_valid:
                self.precision_size_to_tp[index] += 1
            else:
                self.precision_size_to_fp[index] += 1

    def update_recall(self, num_bins: int, event_valid: Optional[bool], event_type: EventType, interval: Interval, reason: str, sample_list: List[str]):
        classification_status = "NA" if event_valid is not None else ("TP" if event_valid else "FP")
        classification_reason = reason if reason is not None else "NA"
        self.recall_results_df = self.recall_results_df.append({"Chromosome": interval.chrom,
                                                                "Start": interval.start,
                                                                "End": interval.end,
                                                                "Genotype": event_type.name,
                                                                "NumBins": num_bins,
                                                                "ClassificationStatus": classification_status,
                                                                "ClassificationReason": classification_reason,
                                                                "Samples": ",".join(sample_list)},
                                                               ignore_index=True)
        if event_valid is not None:
            index = min(num_bins - 1, len(self.event_size_bins) - 1)
            if event_valid:
                self.recall_size_to_tp[index] += 1
            else:
                self.recall_size_to_fn[index] += 1

    def write_to_file(self, path: str):
        precision_file_output = os.path.join(path, self.tool_name + ".precision.confusion_table.tsv")
        self.precision_results_df.to_csv(path_or_buf=open(precision_file_output, 'w'), sep='\t', index=False)
        recall_file_output = os.path.join(path, self.tool_name + ".recall.confusion_table.tsv")
        self.recall_results_df.to_csv(path_or_buf=open(recall_file_output, 'w'), sep='\t', index=False)

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
