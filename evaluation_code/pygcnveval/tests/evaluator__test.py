from callset import GCNVCallset, TruthCallset
from interval_collection import IntervalCollection
from evaluator import PerEventEvaluator
from evaluation_result import EvaluationResult
import pandas as pd

pd.set_option("display.max_rows", None, "display.max_columns", None)

# analyzed intervals
ANALYZED_INTERVALS = "./test_files/analyzed_intervals.interval_list"

# gCNV test callset resources
GCNV_CALLSET_TEST_VCF_LIST = ["./test_files/GCNV_SAMPLE_0.vcf", "./test_files/GCNV_SAMPLE_1.vcf",
                              "./test_files/GCNV_SAMPLE_2.vcf"]

# Truth test callset resources
TRUTH_CALLSET_TEST_BED = "./test_files/truth.bed"


# What the results should be
EVALUATION_RESULT_EXPECTED = EvaluationResult()
EVALUATION_RESULT_EXPECTED.precision_size_to_tp[1] = 3
EVALUATION_RESULT_EXPECTED.precision_size_to_fp[0] = 2
EVALUATION_RESULT_EXPECTED.precision_size_to_fp[1] = 1

EVALUATION_RESULT_EXPECTED.recall_size_to_tp[2] = 3
EVALUATION_RESULT_EXPECTED.recall_size_to_fn[1] = 3


def test_evaluator():
    interval_collection = IntervalCollection.read_interval_list(ANALYZED_INTERVALS)
    truth_callset = TruthCallset.read_in_callset(truth_callset_bed=TRUTH_CALLSET_TEST_BED,
                                                 interval_collection=interval_collection)
    truth_callset.filter_out_uncovered_events(interval_collection, overall_reciprocal_overlap=0.3, padding=0)
    gcnv_callset = GCNVCallset.read_in_callset(gcnv_segment_vcfs=GCNV_CALLSET_TEST_VCF_LIST)
    evaluator = PerEventEvaluator(truth_callset=truth_callset)
    evaluation_result_actual = evaluator.evaluate_callset_against_the_truth(gcnv_callset=gcnv_callset, minimum_reciprocal_overlap=0.4)
    assert evaluation_result_actual == EVALUATION_RESULT_EXPECTED


def main():
    test_evaluator()


if __name__ == '__main__':
    main()
