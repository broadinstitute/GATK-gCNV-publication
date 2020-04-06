from callset import GCNVCallset, TruthCallset
from interval_collection import IntervalCollection
from event import EventType
import pandas as pd
import pyranges as pr

pd.set_option("display.max_rows", None, "display.max_columns", None)

# analyzed intervals
ANALYZED_INTERVALS = "./test_files/analyzed_intervals.interval_list"

# gCNV test callset resources
GCNV_CALLSET_TEST_VCF = "./test_files/GCNV_SAMPLE_1.vcf"
GCNV_CALLSET_TEST_VALUES = [('1', 1001, 3000, EventType.DEL, 2, 60),
                            ('1', 3001, 10000, EventType.REF, 4, 100),
                            ('2', 4001, 5000, EventType.DUP, 1, 100),
                            ('2', 6001, 7000, EventType.REF, 1, 20)]
GCNV_CALLSET_SAMPLE_NAME = "SAMPLE_1"
GCNV_CALLSET_TEST_DF = pd.DataFrame(GCNV_CALLSET_TEST_VALUES, columns=GCNVCallset.vcf_columns)
GCNV_CALLSET_TEST_DF = GCNV_CALLSET_TEST_DF.astype(GCNVCallset.vcf_column_types)
GCNV_CALLSET_TEST_PYRANGE_EXPECTED = pr.PyRanges(GCNV_CALLSET_TEST_DF)

# Truth test callset resources
TRUTH_CALLSET_TEST_BED = "./test_files/truth.bed"
TRUTH_CALLSET_VALUES = [('1', 501, 4500, 'DEL_chr1_1', EventType.DEL, frozenset(['SAMPLE_0', 'SAMPLE_1', 'SAMPLE_2']), 1.0, 3),
                        ('1', 7001, 10000, 'DEL_chr1_2', EventType.DEL, frozenset(['SAMPLE_0']), 1./3, 2),
                        ('2', 1001, 3000, 'DUP_chr2_1', EventType.DUP, frozenset(['SAMPLE_0']), 1./3, 0),
                        ('2', 4001, 7000, 'DUP_chr2_2', EventType.DUP, frozenset(['SAMPLE_0', 'SAMPLE_1']), 2./3, 2)]
TRUTH_CALLSET_TEST_DF = pd.DataFrame(TRUTH_CALLSET_VALUES, columns=TruthCallset.bed_columns)
TRUTH_CALLSET_TEST_DF = TRUTH_CALLSET_TEST_DF.astype(TruthCallset.bed_column_types)
TRUTH_CALLSET_TEST_PYRANGE_EXPECTED = pr.PyRanges(TRUTH_CALLSET_TEST_DF)
TRUTH_CALLSET_TEST_FILTERED_PYRANGE_EXPECTED = pr.PyRanges(TRUTH_CALLSET_TEST_DF.iloc[[0, 1, 3]])


def test_gcnv_callset():
    gcnv_callset_actual = GCNVCallset.read_in_callset(gcnv_segment_vcfs=[GCNV_CALLSET_TEST_VCF])
    assert gcnv_callset_actual.sample_to_pyrange_map[GCNV_CALLSET_SAMPLE_NAME].df.equals(GCNV_CALLSET_TEST_PYRANGE_EXPECTED.df)


def test_truth_callset():
    interval_collection = IntervalCollection.read_interval_list(ANALYZED_INTERVALS)
    # Test callset parsing
    truth_callset_actual = TruthCallset.read_in_callset(truth_callset_bed=TRUTH_CALLSET_TEST_BED,
                                                        interval_collection=interval_collection)
    assert truth_callset_actual.truth_callset_pyrange.df.equals(TRUTH_CALLSET_TEST_PYRANGE_EXPECTED.df)
    # Test filtering

    truth_callset_actual.filter_out_uncovered_events(interval_collection, overall_reciprocal_overlap=0.3, padding=0)
    assert truth_callset_actual.truth_callset_pyrange.df.equals(TRUTH_CALLSET_TEST_FILTERED_PYRANGE_EXPECTED.df)


def main():
    test_gcnv_callset()
    test_truth_callset()


if __name__ == '__main__':
    main()
