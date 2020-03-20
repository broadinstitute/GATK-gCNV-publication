# Workflow for collecting coverage profiles on a list of BAMs 
#
#############

import "https://api.firecloud.org/ga4gh/v1/tools/broad-firecloud-dsde:GATK_Germline_CNV_Validation_Common_Tasks/versions/16/plain-WDL/descriptor" as CNVTasks

workflow CNVGermlineCaseWorkflow {

    ##################################
    #### required basic arguments ####
    ##################################
    File intervals
    Array[String]+ normal_bams
    Array[String]+ normal_bais
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    String gatk_docker

    ##################################
    #### optional basic arguments ####
    ##################################
    File? gatk4_jar_override
    Int? preemptible_attempts

    ####################################################
    #### optional arguments for PreprocessIntervals ####
    ####################################################
    Int? padding
    Int? bin_length

    ##############################################
    #### optional arguments for CollectCounts ####
    ##############################################
    String? collect_counts_format
    Int? mem_gb_for_collect_counts
    
    Array[Pair[String, String]] normal_bams_and_bais = zip(normal_bams, normal_bais)
    
    call CNVTasks.PreprocessIntervals {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            padding = padding,
            bin_length = bin_length,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    scatter (normal_bam_and_bai in normal_bams_and_bais) {
        call CNVTasks.CollectCounts {
            input:
                intervals = PreprocessIntervals.preprocessed_intervals,
                bam = normal_bam_and_bai.left,
                bam_idx = normal_bam_and_bai.right,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                format = collect_counts_format,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_collect_counts,
                preemptible_attempts = preemptible_attempts
        }
    }
    
    output {
        File preprocessed_intervals = PreprocessIntervals.preprocessed_intervals
        Array[File] read_counts_entity_ids = CollectCounts.entity_id
        Array[File] read_counts = CollectCounts.counts
   }
}
