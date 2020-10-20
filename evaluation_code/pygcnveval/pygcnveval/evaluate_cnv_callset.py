import argparse
from typing import List

from callset import TruthCallset, GCNVCallset, XHMMCallset
from interval_collection import IntervalCollection
from evaluator import PerEventEvaluator, PerBinEvaluator
import plotting


def evaluate_cnv_callsets_and_plot_results(analyzed_intervals: str, truth_callset_bed: str, gcnv_vcfs: List[str],
                                           xhmm_vcf: str, output_directory: str, minimum_overlap: float,
                                           min_sq_threshold: int):
    print("Reading in interval list...", flush=True)
    interval_collection = IntervalCollection.read_interval_list(analyzed_intervals)
    print("Reading in gCNV callset...", flush=True)
    gcnv_callset = GCNVCallset.read_in_callset(gcnv_segment_vcfs=gcnv_vcfs, interval_collection=interval_collection)
    print("Reading in XHMM callset", flush=True)
    xhmm_callset = XHMMCallset.read_in_callset(xhmm_vcf=xhmm_vcf,
                                               interval_collection=interval_collection,
                                               samples_to_keep=gcnv_callset.sample_set)
    print("Reading in truth callset...", flush=True)
    truth_callset = TruthCallset.read_in_callset(truth_callset_bed_file=truth_callset_bed,
                                                 interval_collection=interval_collection,
                                                 samples_to_keep=gcnv_callset.sample_set)
    print("Filtering truth callset...", flush=True)
    truth_callset.filter_out_uncovered_events(interval_collection)
    print("Filtering XHMM callset...", flush=True)
    xhmm_callset.filter_out_uncovered_events(interval_collection)
    plotting.plot_and_save_callset_event_distribution_plots(gcnv_callset, output_directory)
    plotting.plot_and_save_callset_event_distribution_plots(truth_callset, output_directory)
    plotting.plot_and_save_callset_event_distribution_plots(xhmm_callset, output_directory)


    rare_intervals_subset = truth_callset.subset_intervals_to_rare_regions(interval_collection,
                                                                           max_allelic_fraction=0.01)
    common_intervals_subset = IntervalCollection(interval_collection.pyrange.subtract(rare_intervals_subset.pyrange))
    print("Performing gCNV per event evaluation...", flush=True)
    per_event_evaluator_gcnv = PerEventEvaluator(truth_callset=truth_callset, callset=gcnv_callset)
    print("%s/%s matching samples are found in the truth callset" %
          (len(per_event_evaluator_gcnv.sample_list_to_eval), len(xhmm_callset.sample_set)))
    per_event_evaluation_result_gcnv = per_event_evaluator_gcnv.evaluate_callset_against_the_truth(minimum_overlap=minimum_overlap,
                                                                                                   min_quality_threshold=min_sq_threshold)
    print("Performing XHMM per event evaluation...", flush=True)
    per_event_evaluator_xhmm = PerEventEvaluator(truth_callset=truth_callset, callset=xhmm_callset)
    per_event_evaluation_result_xhmm = per_event_evaluator_xhmm.evaluate_callset_against_the_truth(minimum_overlap=minimum_overlap,
                                                                                                   min_quality_threshold=min_sq_threshold)
    plotting.plot_and_save_per_event_evaluation_results([per_event_evaluation_result_gcnv,
                                                         per_event_evaluation_result_xhmm], output_directory)


    #print("Performing per bin evaluation...", flush=True)
    #per_bin_evaluator = PerBinEvaluator(truth_callset=truth_callset, gcnv_callset=gcnv_callset,
    #                                   interval_collection=rare_intervals_subset)
    # TODO pass an optional number of PR curve points parameter
    #per_bin_evaluation_result = per_bin_evaluator.evaluate_callset_against_the_truth(gcnv_callset)
    #plotting.plot_and_save_per_bin_evaluation_results(per_bin_evaluation_result, output_directory)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--output_dir', metavar='OutputDirectory', type=str,
                        help='Output directory.', required=True)

    parser.add_argument('--gcnv_segment_vcfs', metavar='gCNVSegmentVCF', type=str, nargs='+',
                        help='Segment VCFs output by gCNV.', required=True)

    parser.add_argument('--xhmm_vcf', metavar='xhmmVCF', type=str,
                        help='Path to XHMM VCF.')

    parser.add_argument('--sorted_truth_calls_bed', metavar='SortedTruthCallsBed', type=str,
                        help='Sorted bed file that contains truth calls on superset of samples that are being evaluated',
                        required=True)

    parser.add_argument('--analyzed_intervals', metavar='AnalyzedIntervalsFile', type=str,
                        help='Interval list file that was given as input to the CNV calling tools '
                             '(in case of gCNV that is the padded, but not filtered interval list.', required=True)

    parser.add_argument('--min_required_overlap', metavar='MinimumRequiredOverlap', type=float,
                        help='Minimum required overlap (non-reciprocal) to validate an event.', required=True)

    parser.add_argument('--min_sq_threshold', metavar='MinimumSQThreshold', type=int,
                        help='SQ threshold to filter gCNV events on.', required=True)

    args = parser.parse_args()
    output_dir = args.output_dir
    gcnv_segment_vcfs = args.gcnv_segment_vcfs
    xhmm_vcf = args.xhmm_vcf
    truth_callset = args.sorted_truth_calls_bed
    analyzed_intervals = args.analyzed_intervals
    min_required_overlap = args.min_required_overlap
    min_sq_threshold = args.min_sq_threshold

    evaluate_cnv_callsets_and_plot_results(analyzed_intervals, truth_callset, gcnv_segment_vcfs, xhmm_vcf, output_dir,
                                           min_required_overlap, min_sq_threshold)


if __name__ == '__main__':
    main()