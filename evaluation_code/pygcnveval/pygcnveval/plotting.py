import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from cycler import cycler
import pandas as pd
import os
import numpy as np
import math
from typing import List

from evaluation_result import PerEventEvaluationResult, PerBinEvaluationResult
from callset import Callset

plt.style.use('seaborn')
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = '5'
mpl.rcParams['figure.dpi'] = 80
mpl.rcParams['axes.labelsize'] = 13.0
mpl.rcParams['axes.titlesize'] = 17.0
mpl.rcParams['axes.prop_cycle'] = cycler('color', ['#5d769c', '#55A868', '#C44E52', '#8172B2', '#CCB974', '#64B5CD'])
mpl.rcParams['figure.figsize'] = [9., 7.]


def plot_and_save_per_event_evaluation_results(evaluation_result_list: List[PerEventEvaluationResult], output_directory: str):
    bins = evaluation_result_list[0].event_size_bins
    for evaluation_result in evaluation_result_list:
        assert evaluation_result.event_size_bins == bins, "Evaluation results are stratified by bin length differently"

    def sum_over_dict(d: dict, last_index: int, bins: list):
        return sum([d[i] for i in bins[last_index:]])

    # calculate precision, recall and f-1 scores
    precisions_for_bins, recall_for_bins, f_1_for_bins = [], [], []
    for index, evaluation_result in enumerate(evaluation_result_list):
        tp = evaluation_result.precision_size_to_tp
        fp = evaluation_result.precision_size_to_fp
        print(tp)
        print(fp)
        precisions_for_bins.append([sum_over_dict(tp, i, bins) / max(1., (sum_over_dict(tp, i, bins)+sum_over_dict(fp, i, bins))) for i in bins])
        tp = evaluation_result.recall_size_to_tp
        fn = evaluation_result.recall_size_to_fn
        print(tp)
        print(fn)
        recall_for_bins.append([sum_over_dict(tp, i, bins) / max(1., (sum_over_dict(tp, i, bins) + sum_over_dict(fn, i, bins))) for i in bins])
        f_1_for_bins.append([2 / ((1. / recall_for_bins[index][i]) + (1. / precisions_for_bins[index][i])) for i in bins])

    # plot precision
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[1, 1.8])
    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, :])
    num_results = len(evaluation_result_list)
    width = 0.8 / num_results
    offsets = np.linspace(-width / 2, width / 2, num_results) if num_results > 1 else [0]

    for index, evaluation_result in enumerate(evaluation_result_list):
        tp = evaluation_result.precision_size_to_tp
        fp = evaluation_result.precision_size_to_fp

        called_events_in_bins = [tp[i] + fp[i] for i in bins]
        ax1.bar(np.array(bins) + offsets[index], called_events_in_bins, width=width)
        ax2.scatter(bins, precisions_for_bins[index], label=evaluation_result.tool_name)

    ax1.set_title("Precision stratified by number of overlapping bins")
    ax1.set_ylabel("Number called events")
    ax2.set_xlabel(">= number of bins")
    ax2.set_ylabel("Precision")
    ax2.set_ylim(0, 1.)
    ax2.legend()

    plt.savefig(os.path.join(output_directory, "precision.png"))
    plt.close()

    # Plot recall
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[1, 1.8])
    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, :])

    for index, evaluation_result in enumerate(evaluation_result_list):
        ax2.scatter(bins, recall_for_bins[index], label=evaluation_result.tool_name)

    true_events_in_bins = [evaluation_result_list[0].recall_size_to_tp[i]
                           + evaluation_result_list[0].recall_size_to_fn[i] for i in bins]
    # for evaluation_result in evaluation_result_list:
    #     assert true_events_in_bins == [evaluation_result.recall_size_to_tp[i]
    #                                    + evaluation_result.recall_size_to_fn[i] for i in bins],\
    #         "Subset of considered truth events was not same for all evaluated tools."

    ax1.bar(bins, true_events_in_bins)
    ax1.set_title("Recall stratified by number of overlapping bins")
    ax1.set_ylabel("Number of true events")
    ax2.set_xlabel(">= number of bins")
    ax2.set_ylabel("Recall")
    ax2.set_ylim(0, 1.)
    ax2.legend()

    plt.savefig(os.path.join(output_directory, "recall.png"))
    plt.close()

    # Plot f-1 score
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(gs[0, :])
    for index, evaluation_result in enumerate(evaluation_result_list):
        ax1.scatter(bins, f_1_for_bins[index], label=evaluation_result.tool_name)

    ax1.set_title("f-1 score stratified by number of overlapping bins")
    ax1.set_xlabel(">= number of bins")
    ax1.set_ylabel("f-1 score")
    ax1.set_ylim(0, 1.)
    ax1.legend()

    plt.savefig(os.path.join(output_directory, "f_1.png"))
    plt.close()


def plot_and_save_callset_event_distribution_plots(callset_: Callset, output_directory: str):
    callset_name = callset_.get_name()
    callset_num_events_distr = callset_.get_callset_num_events_distribution()
    plt.hist(callset_num_events_distr, bins=(int(max(callset_num_events_distr) / 2)))
    plt.xlabel('Number of events')
    plt.title(callset_name + ", number of events distribution")
    plt.savefig(os.path.join(output_directory, callset_name + "_num_events_distribution.png"))
    plt.close()

    callset_event_size_distr = callset_.get_callset_event_size_distribution()
    plt.hist(callset_event_size_distr, bins=(int(max(callset_event_size_distr) / 2)))
    plt.xlabel('Event size in bins')
    plt.title(callset_name + ", event size distribution")
    plt.savefig(os.path.join(output_directory, callset_name + "_event_size_distribution.png"))
    plt.close()


def plot_and_save_per_bin_evaluation_results(evaluation_result: PerBinEvaluationResult, output_directory: str):
    plt.plot(evaluation_result.get_recall_array(), evaluation_result.get_precision_array(), 'o')
    plt.title("Precision-Recall curve")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.1, 1.1)
    plt.savefig(os.path.join(output_directory, "precision_recall_curve.png"))
    plt.close()


def plot_number_of_events_distribution(output_file: str, gcnv_segment_vcfs: list):
    variant_num_list = []
    for gcnv_vcf in gcnv_segment_vcfs:
        variant_num_list.append(count_num_variants(gcnv_vcf))
    fig = plt.figure(figsize=(10,6))
    plt.title("Event number distribution")
    plt.hist(variant_num_list, bins=100, edgecolor='black', linewidth=1)
    fig.savefig(output_file, dpi=120)
    print("Mean number of segments in each sample: ")
    print(np.mean(np.array(variant_num_list)))
    print("\n")


def count_num_variants(vcf_file):
    num_vars = 0
    with open(vcf_file, 'r') as vcf:
        for line in vcf.readlines():
            if not line.startswith("#"):
                num_vars += 1
    return num_vars


def plot_ard_components(models_dir, num_shards):
    num_plots_per_row = 4
    num_rows = math.ceil(num_shards/num_plots_per_row)
    plt.figure(figsize=(20, 20))
    for i in range(num_shards):
        ard_file = models_dir + "/shard-" + str(i) + "/mu_ard_u_log__.tsv"
        ard_df = pd.read_csv(open(ard_file, 'r'), comment="@", sep='\t', header=None)
        ard_vector = np.transpose(ard_df.as_matrix())[0]
        plt.subplot(num_rows, num_plots_per_row, i + 1)
        plt.plot(range(ard_vector.size), np.sort(ard_vector))
        plt.plot(range(ard_vector.size), np.zeros(ard_vector.size), 'r--')
        plt.ylim(-2, 2)
        plt.title("Shard " + str(i) + " ARD mean values")
    plt.tight_layout()
