# Evaluate performance of a CNV callset against a truth callset
####

workflow EvaluateCNVCallset {

    ##################################
    #### required basic arguments ####
    ##################################
    Array[File]+ genotyped_segments_vcfs
    File truth_bed
    File analyzed_intervals
    String gcnv_evaluation_docker

    ##################################
    #### optional basic arguments ####
    ##################################
    Int? preemptible_attempts
    String gcnv_eval_script = "/root/evaluate_cnv_callset.py"

    call EvaluateCalls {
        input:
            genotyped_segments_vcfs = genotyped_segments_vcfs,
            truth_bed = truth_bed,
            analyzed_intervals = analyzed_intervals,
            gcnv_evaluation_docker = gcnv_evaluation_docker,
            gcnv_eval_script = gcnv_eval_script,
            preemptible_attempts = preemptible_attempts
    }

    output {
        Array[File] performance_plots = EvaluateCalls.performance_plots
    }
}

task EvaluateCalls {
    Array[File]+ genotyped_segments_vcfs
    File truth_bed
    File analyzed_intervals
    String gcnv_eval_script

    #Runtime parameters
    String gcnv_evaluation_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    command <<<
        mkdir plots

        python ${gcnv_eval_script} \
          --output_dir "./plots/" \
          --gcnv_segment_vcfs ${sep=' ' genotyped_segments_vcfs} \
          --sorted_truth_calls_bed ${truth_bed} \
          --analyzed_intervals ${analyzed_intervals}
    >>>

    runtime {
        docker: "${gcnv_evaluation_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        Array[File] performance_plots = glob("./plots/*")
    }
}