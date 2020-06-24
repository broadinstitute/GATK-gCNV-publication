import argparse
from typing import List

from callset import TruthCallset, GCNVCallset
from interval_collection import IntervalCollection
from evaluator import PerEventEvaluator, PerBinEvaluator
import plotting


def evaluate_cnv_callsets_and_plot_results(analyzed_intervals: str, truth_callset_bed: str, gcnv_vcfs: List[str],
                                           output_directory: str):
    print("Reading in interval list...", flush=True)
    interval_collection = IntervalCollection.read_interval_list(analyzed_intervals)
    print("Reading in gCNV callset...", flush=True)
    gcnv_callset = GCNVCallset.read_in_callset(gcnv_segment_vcfs=gcnv_vcfs)
    print("Reading in truth callset...", flush=True)
    truth_callset = TruthCallset.read_in_callset(truth_callset_bed_file=truth_callset_bed,
                                                 interval_collection=interval_collection,
                                                 samples_to_keep=gcnv_callset.sample_set)
    print("Filtering truth callset...", flush=True)
    truth_callset.filter_out_uncovered_events(interval_collection)

    print("Performing per event evaluation...", flush=True)
    per_event_evaluator = PerEventEvaluator(truth_callset=truth_callset)
    per_event_evaluation_result = per_event_evaluator.evaluate_callset_against_the_truth(gcnv_callset=gcnv_callset)
    plotting.plot_and_save_per_event_evaluation_results(per_event_evaluation_result, output_directory)

    print("Performing per event evaluation...", flush=True)
    rare_intervals_subset = truth_callset.subset_intervals_to_rare_regions(interval_collection, max_allelic_fraction=0.01)
    per_bin_evaluator = PerBinEvaluator(truth_callset=truth_callset, interval_collection=rare_intervals_subset)
    # TODO pass an optional number of PR curve points parameter
    per_bin_evaluation_result = per_bin_evaluator.evaluate_callset_against_the_truth(gcnv_callset)
    plotting.plot_and_save_per_bin_evaluation_results(per_bin_evaluation_result, output_directory)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--output_dir', metavar='OutputDirectory', type=str,
                        help='output directory', required=True)

    parser.add_argument('--gcnv_segment_vcfs', metavar='gCNVSegmentVCF', type=str, nargs='+',
                        help='Segment VCFs output by gCNV', required=True)

    parser.add_argument('--sorted_truth_calls_bed', metavar='SortedTruthCallsBed', type=str,
                        help='Sorted bed file that contains truth calls on superset of samples that are being evaluated',
                        required=True)

    parser.add_argument('--analyzed_intervals', metavar='AnalyzedIntervalsFile', type=str,
                        help='Interval list file that was given as input to the CNV calling tools '
                             '(in case of gCNV that is the padded, but not filtered interval list', required=True)

    args = parser.parse_args()
    output_dir = args.output_dir
    gcnv_segment_vcfs = args.gcnv_segment_vcfs
    truth_callset = args.sorted_truth_calls_bed
    analyzed_intervals = args.analyzed_intervals

    evaluate_cnv_callsets_and_plot_results(analyzed_intervals, truth_callset, gcnv_segment_vcfs, output_dir)


if __name__ == '__main__':
    main()