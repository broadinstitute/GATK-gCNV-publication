import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import math

from evaluation_result import PerEventEvaluationResult, PerBinEvaluationResult


def plot_and_save_per_event_evaluation_results(evaluation_result: PerEventEvaluationResult, output_directory: str):
    bins = evaluation_result.event_size_bins

    # Plot precision
    tp = evaluation_result.precision_size_to_tp
    fp = evaluation_result.precision_size_to_fp
    precisions_for_bins = [tp[i]/max((tp[i]+fp[i]), 1) for i in bins]
    plt.plot(bins, precisions_for_bins)
    plt.title("Precision stratified by exon number")
    plt.xlabel("Number of exons")
    plt.ylabel("Precision")
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0, 1.)
    plt.savefig(os.path.join(output_directory, "precision.png"))
    plt.close()
    # Plot recall
    tp = evaluation_result.recall_size_to_tp
    fn = evaluation_result.recall_size_to_fn
    recall_for_bins = [tp[i]/max((tp[i]+fn[i]), 1) for i in bins]
    plt.plot(bins, recall_for_bins)
    plt.title("Recall stratified by interval number")
    plt.xlabel("Number of exons")
    plt.ylabel("Recall")
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0, 1.)
    plt.savefig(os.path.join(output_directory, "recall.png"))
    plt.close()


def plot_and_save_per_bin_evaluation_results(evaluation_result: PerBinEvaluationResult, output_directory: str):
    plt.plot(evaluation_result.get_recall_array(), evaluation_result.get_precision_array(), 'o')
    plt.title("Precision stratified by exon number")
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
    plt.figure(figsize=(20,20))
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
