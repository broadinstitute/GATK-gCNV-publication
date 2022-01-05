import argparse
from typing import List

from callset import TruthCallset, GCNVCallset, XHMMCallset
from interval_collection import IntervalCollection
from evaluator import PerEventEvaluator, PerBinEvaluator, PerSiteEvaluator
import plotting


def evaluate_cnv_callsets_and_plot_results(analyzed_intervals: str,
                                           truth_callset_bed: str,
                                           gcnv_vcfs: List[str],
                                           gcnv_callset_tsv: str,
                                           gcnv_max_event_number: int,
                                           gcnv_joint_vcf: str,
                                           xhmm_vcfs: List[str],
                                           output_directory: str,
                                           minimum_overlap: float,
                                           gcnv_sq_min_del: int,
                                           gcnv_sq_min_dup: int,
                                           samples_to_evaluate_path: str):
    perform_per_site_eval = False

    print("Reading in interval list...", flush=True)
    interval_collection = IntervalCollection.read_interval_list(analyzed_intervals)
    callsets_to_evaluate = []
    if gcnv_vcfs or gcnv_callset_tsv:
        print("Reading in gCNV callset...", flush=True)
        gcnv_callset = GCNVCallset.read_in_callset(gcnv_segment_vcfs=gcnv_vcfs,
                                                   gcnv_callset_tsv=gcnv_callset_tsv,
                                                   gcnv_joint_vcf=gcnv_joint_vcf,
                                                   interval_collection=interval_collection,
                                                   max_events_allowed=gcnv_max_event_number)
        callsets_to_evaluate.append(gcnv_callset)

    if xhmm_vcfs:
        print("Reading in XHMM callset", flush=True)
        xhmm_callset = XHMMCallset.read_in_callset(xhmm_vcfs=xhmm_vcfs,
                                                   interval_collection=interval_collection,
                                                   samples_to_keep=gcnv_callset.sample_set)
        callsets_to_evaluate.append(xhmm_callset)

    if samples_to_evaluate_path:
        with open(samples_to_evaluate_path, 'r') as f:
            sample_list_to_evaluate = [s.strip() for s in f.readlines()]
    else:
        sample_list_to_evaluate = set.intersection(*[set(c.sample_set) for c in callsets_to_evaluate])

    print("Reading in truth callset...", flush=True)

    truth_callset = TruthCallset.read_in_callset(truth_callset_bed_file=truth_callset_bed,
                                                 interval_collection=interval_collection,
                                                 samples_to_keep=sample_list_to_evaluate)
    sample_list_to_evaluate = list(set(sample_list_to_evaluate).intersection(truth_callset.sample_set))
    print("Filtering truth callset...", flush=True)
    truth_callset.filter_out_uncovered_events_from_joint_callset(interval_collection)

    if perform_per_site_eval:
        print("Performing gCNV per site evaluation...", flush=True)
        per_site_evaluator_gcnv = PerSiteEvaluator(truth_callset=truth_callset, callset=gcnv_callset)
        print("%s/%s matching samples are found in the truth callset" %
              (len(per_site_evaluator_gcnv.sample_list_to_eval), len(gcnv_callset.sample_set)))
        per_site_evaluation_result_gcnv = per_site_evaluator_gcnv.evaluate_callset_against_the_truth()

        plotting.plot_and_save_per_event_evaluation_results([per_site_evaluation_result_gcnv], output_directory)
        per_site_evaluation_result_gcnv.write_to_file(output_directory)

    #sample_list_to_evaluate = set(sample_list_to_evaluate)
    per_event_evaluation_results = []
    for callset in callsets_to_evaluate:
        print("Performing per event evaluation on callset: {0}...".format(callset.get_name()), flush=True)
        per_event_evaluator = PerEventEvaluator(truth_callset=truth_callset,
                                                callset=callset,
                                                sample_list_to_evaluate=sample_list_to_evaluate)
        print("%s/%s matching samples from {0} are found in the truth callset".format(callset.get_name()) %
              (len(per_event_evaluator.sample_list_to_eval), len(callset.sample_set)))
        per_event_evaluation_result = per_event_evaluator.evaluate_callset_against_the_truth(minimum_overlap=minimum_overlap,
                                                                                             gcnv_sq_min_del=gcnv_sq_min_del,
                                                                                             gcnv_sq_min_dup=gcnv_sq_min_dup)
        per_event_evaluation_results.append(per_event_evaluation_result)

        per_event_evaluation_result.write_to_file(output_directory)

    plotting.plot_and_save_per_event_evaluation_results(per_event_evaluation_results, output_directory)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--output_dir', metavar='OutputDirectory', type=str,
                        help='Output directory.', required=True)

    parser.add_argument('--gcnv_segment_vcfs', metavar='gCNVSegmentVCF', type=str, nargs='+',
                        help='Segment VCFs output by gCNV.', required=False)

    parser.add_argument('--gcnv_callset_tsv', metavar='gCNVCallsetTSV', type=str,
                        help='gCNV callset in custom TSV format', required=False)

    parser.add_argument('--gcnv_max_event_number', metavar='gCNVMaxEventNumber', type=int,
                        help='Cutoff for maximum number of events that are kept for samples in the gCNV callset',
                        required=True)

    parser.add_argument('--gcnv_joint_vcf', metavar='gCNVJointVCF', type=str,
                        help='Jointly genotyped gCNV callset.', required=False)

    parser.add_argument('--xhmm_vcfs', metavar='xhmmVCF', type=str, nargs='+',
                        help='Paths to the XHMM VCFs.')

    parser.add_argument('--sorted_truth_calls_bed', metavar='SortedTruthCallsBed', type=str,
                        help='Sorted bed file that contains truth calls on superset of samples that are being evaluated',
                        required=True)

    parser.add_argument('--analyzed_intervals', metavar='AnalyzedIntervalsFile', type=str,
                        help='Interval list file that was given as input to the CNV calling tools '
                             '(in case of gCNV that is the padded, but not filtered interval list.', required=True)

    parser.add_argument('--min_required_overlap', metavar='MinimumRequiredOverlap', type=float,
                        help='Minimum required overlap (non-reciprocal) to validate an event.', required=True)

    parser.add_argument('--gcnv_min_sq_del_threshold', metavar='MinimumSQDelThreshold', type=int,
                        help='SQ threshold to filter gCNV deletion events on.', required=True)

    parser.add_argument('--gcnv_min_sq_dup_threshold', metavar='MinimumSQDelThreshold', type=int,
                        help='SQ threshold to filter gCNV duplication events on', required=True)

    parser.add_argument('--samples_to_evaluate_path', metavar='SamplesToEvaluate', type=str,
                        help='A file containing the set of samples to evaluate, one sample per line.',
                        required=False)

    args = parser.parse_args()
    output_dir = args.output_dir
    gcnv_segment_vcfs = args.gcnv_segment_vcfs
    gcnv_callset_tsv = args.gcnv_callset_tsv
    gcnv_max_event_number = args.gcnv_max_event_number
    gcnv_joint_vcf = args.gcnv_joint_vcf
    xhmm_vcfs = args.xhmm_vcfs
    truth_callset = args.sorted_truth_calls_bed
    analyzed_intervals = args.analyzed_intervals
    min_required_overlap = args.min_required_overlap
    gcnv_sq_min_del = args.gcnv_min_sq_del_threshold
    gcnv_sq_min_dup = args.gcnv_min_sq_dup_threshold
    samples_to_evaluate_path = args.samples_to_evaluate_path

    assert (gcnv_segment_vcfs is None) ^ (gcnv_callset_tsv is None),\
        "Exactly one of the gCNV segment VCF list or gCNV TSV callset must be defined"

    evaluate_cnv_callsets_and_plot_results(analyzed_intervals, truth_callset, gcnv_segment_vcfs, gcnv_callset_tsv, gcnv_max_event_number,
                                           gcnv_joint_vcf, xhmm_vcfs, output_dir, min_required_overlap,
                                           gcnv_sq_min_del, gcnv_sq_min_dup, samples_to_evaluate_path)


if __name__ == '__main__':
    main()