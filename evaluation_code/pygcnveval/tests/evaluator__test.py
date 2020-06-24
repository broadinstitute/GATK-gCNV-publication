from callset import GCNVCallset, TruthCallset
from interval_collection import IntervalCollection
from evaluator import PerEventEvaluator, PerBinEvaluator
from evaluation_result import PerEventEvaluationResult, PerBinEvaluationResult
import pandas as pd
import numpy as np

pd.set_option("display.max_rows", None, "display.max_columns", None)

# analyzed intervals
ANALYZED_INTERVALS = "./test_files/analyzed_intervals.interval_list"

# gCNV test callset resources
GCNV_CALLSET_TEST_VCF_LIST = ["./test_files/GCNV_SAMPLE_0.vcf", "./test_files/GCNV_SAMPLE_1.vcf",
                              "./test_files/GCNV_SAMPLE_2.vcf"]

# Truth test callset resources
TRUTH_CALLSET_TEST_BED = "./test_files/truth.bed"


# Per event evaluator expected results
PER_EVENT_EVALUATION_RESULT_EXPECTED = PerEventEvaluationResult()
PER_EVENT_EVALUATION_RESULT_EXPECTED.precision_size_to_tp[1] = 3
PER_EVENT_EVALUATION_RESULT_EXPECTED.precision_size_to_fp[0] = 2
PER_EVENT_EVALUATION_RESULT_EXPECTED.precision_size_to_fp[1] = 1

PER_EVENT_EVALUATION_RESULT_EXPECTED.recall_size_to_tp[2] = 3
PER_EVENT_EVALUATION_RESULT_EXPECTED.recall_size_to_fn[0] = 3
PER_EVENT_EVALUATION_RESULT_EXPECTED.recall_size_to_fn[1] = 3

# Per bin evaluator expected results
PER_BIN_EVALUATION_RESULT_EXPECTED = PerBinEvaluationResult(np.zeros(1), 1)  # initialize with fake values
PER_BIN_EVALUATION_RESULT_EXPECTED.quality_thresholds = [0.0, 20.0, 60.0, 100.0]
PER_BIN_EVALUATION_RESULT_EXPECTED.true_positives = {0.0: 7, 20.0: 7, 60.0: 7, 100.0: 1}
PER_BIN_EVALUATION_RESULT_EXPECTED.false_positives = {0.0: 3, 20.0: 2, 60.0: 0, 100.0: 0}
PER_BIN_EVALUATION_RESULT_EXPECTED.false_negatives = {0.0: 8, 20.0: 8, 60.0: 10, 100.0: 16}


def test_event_evaluator():
    interval_collection = IntervalCollection.read_interval_list(ANALYZED_INTERVALS)
    gcnv_callset = GCNVCallset.read_in_callset(gcnv_segment_vcfs=GCNV_CALLSET_TEST_VCF_LIST)
    truth_callset = TruthCallset.read_in_callset(truth_callset_bed_file=TRUTH_CALLSET_TEST_BED,
                                                 interval_collection=interval_collection,
                                                 samples_to_keep=gcnv_callset.sample_set)
    truth_callset.filter_out_uncovered_events(interval_collection, min_overlap_fraction=0.3)
    evaluator = PerEventEvaluator(truth_callset=truth_callset)
    evaluation_result_actual = evaluator.evaluate_callset_against_the_truth(gcnv_callset=gcnv_callset, minimum_reciprocal_overlap=0.4)
    assert evaluation_result_actual == PER_EVENT_EVALUATION_RESULT_EXPECTED


def test_bin_evaluator():
    interval_collection = IntervalCollection.read_interval_list(ANALYZED_INTERVALS)
    gcnv_callset = GCNVCallset.read_in_callset(gcnv_segment_vcfs=GCNV_CALLSET_TEST_VCF_LIST)
    truth_callset = TruthCallset.read_in_callset(truth_callset_bed_file=TRUTH_CALLSET_TEST_BED,
                                                 interval_collection=interval_collection,
                                                 samples_to_keep=gcnv_callset.sample_set)

    evaluator = PerBinEvaluator(truth_callset=truth_callset, interval_collection=interval_collection)
    evaluation_result_actual = evaluator.evaluate_callset_against_the_truth(gcnv_callset, 4)
    assert evaluation_result_actual == PER_BIN_EVALUATION_RESULT_EXPECTED


def main():
    test_event_evaluator()
    test_bin_evaluator()


if __name__ == '__main__':
    main()
