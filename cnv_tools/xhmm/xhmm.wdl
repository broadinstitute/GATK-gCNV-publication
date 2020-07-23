#Workflow for running XHMM pipeline

workflow XHMM { 
    Array[String]+ normal_bams
    Array[String]+ normal_bais
    File interval_list
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String xhmm_docker
    String samtools_docker_image
    String cohort_id 
    File discover_parameters_file

    # Filtering parameters and resources
    File? seq_db
    Array[Float]? gc_low_high_filter_params
    Float? interval_complexity_filter_threshold
    Int min_target_size_filter
    Int max_target_size_filter
    Int min_mean_target_rd
    Int max_mean_target_rd
    Int min_mean_sample_rd
    Int max_mean_sample_rd
    Int max_sd_sample_rd

    # Location of gatk, xhmm and plinkseq in the docker
    String gatk_local_jar = "/root/GenomeAnalysisTK-3.8-0.jar"
    String xhmm_executable = "/root/statgen-xhmm-cc14e528d909/xhmm"
    String plinkseq = "/root/plinkseq-0.10/pseq"

    scatter (index in range(length(normal_bams))) {
        if (sub(normal_bams[index], "cram$", "") != normal_bams[index]) {
            call CramToBam {
                input:
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    ref_fasta_dict = ref_fasta_dict,
                    cram = normal_bams[index],
                    cram_idx = normal_bais[index],
                    samtools_docker_image = samtools_docker_image
            }
        }

        call CollectCoverage {
            input:
                bam = select_first([CramToBam.output_bam, normal_bams[index]]),
                bam_idx = select_first([CramToBam.output_bam_index, normal_bais[index]]),
                interval_list = interval_list,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                xhmm_docker = xhmm_docker,
                gatk_local_jar = gatk_local_jar
        }
    }

    call MergeDepths {
        input:
            depth_files = CollectCoverage.read_depth_file,
            xhmm_docker = xhmm_docker,
            xhmm_executable = xhmm_executable,
            cohort_id = cohort_id
    }

    if (defined(gc_low_high_filter_params)) {
        if (length(select_first([gc_low_high_filter_params])) == 2) {
            call GetExtremeGCContentTargets {
                input:
                    interval_list = interval_list,
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    ref_fasta_dict = ref_fasta_dict,
                    gc_low_high_filter_params_final = select_first([gc_low_high_filter_params]),
                    cohort_id = cohort_id,
                    xhmm_docker = xhmm_docker,
                    gatk_local_jar = gatk_local_jar
            }
        }
    }

    if (defined(seq_db)) {
        call GetLowComplexityTargets {
            input:
                seq_db = select_first([seq_db]),
                interval_list = interval_list,
                interval_complexity_filter_threshold = select_first([interval_complexity_filter_threshold, 0.25]),
                cohort_id = cohort_id,
                xhmm_docker = xhmm_docker
        }
    }

    call FilterTargetsAndMeanCenter {
        input:
            merged_depth_file = MergeDepths.merged_depth_file,
            extreme_gc_targets = GetExtremeGCContentTargets.extreme_gc_targets,
            low_complexity_targets = GetLowComplexityTargets.low_complexity_targets,
            xhmm_docker = xhmm_docker,
            xhmm_executable = xhmm_executable,
            cohort_id = cohort_id,
            min_target_size_filter = min_target_size_filter,
            max_target_size_filter = max_target_size_filter,
            min_mean_target_rd = min_mean_target_rd,
            max_mean_target_rd = max_mean_target_rd,
            min_mean_sample_rd = min_mean_sample_rd,
            max_mean_sample_rd = max_mean_sample_rd,
            max_sd_sample_rd = max_sd_sample_rd
    }

    call RunPCA {
        input:
            filtered_centered_counts = FilterTargetsAndMeanCenter.filtered_centered_counts,
            xhmm_docker = xhmm_docker,
            xhmm_executable = xhmm_executable,
            cohort_id = cohort_id
    }

    call NormalizeCounts {
        input:
            panel_of_normals = RunPCA.panel_of_normals,
            panel_of_normals_sd = RunPCA.panel_of_normals_sd,
            filtered_centered_counts = FilterTargetsAndMeanCenter.filtered_centered_counts,
            xhmm_docker = xhmm_docker,
            xhmm_executable = xhmm_executable,
            cohort_id = cohort_id
    }

    call FilterAndZScoreCenter {
        input:
            normalized_counts = NormalizeCounts.normalized_counts,
            xhmm_docker = xhmm_docker,
            xhmm_executable = xhmm_executable,
            cohort_id = cohort_id
    }

    call FilterOriginalReadDepth {
        input:
            merged_depth_file = MergeDepths.merged_depth_file,
            centered_filtered_targets = FilterTargetsAndMeanCenter.mean_center_step_excluded_targets,
            pca_normalized_filtered_targets = FilterAndZScoreCenter.z_score_step_excluded_targets,
            centered_filtered_samples = FilterTargetsAndMeanCenter.mean_center_step_excluded_samples,
            pca_normalized_filtered_samples = FilterAndZScoreCenter.z_score_step_excluded_samples,
            xhmm_docker = xhmm_docker,
            xhmm_executable = xhmm_executable,
            cohort_id = cohort_id
    }

    call DiscoverCNV {
        input:
            normalized_filtered_zcores = FilterAndZScoreCenter.normalized_filtered_zcores,
            original_depth_filtered = FilterOriginalReadDepth.original_depth_filtered,
            discover_parameters_file = discover_parameters_file,
            xhmm_docker = xhmm_docker,
            xhmm_executable = xhmm_executable,
            cohort_id = cohort_id
    }

    call GenotypeCNV {
        input:
            normalized_filtered_zcores = FilterAndZScoreCenter.normalized_filtered_zcores,
            original_depth_filtered = FilterOriginalReadDepth.original_depth_filtered,
            discover_file = DiscoverCNV.discovered_file,
            discover_parameters_file = discover_parameters_file,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            xhmm_docker = xhmm_docker,
            xhmm_executable = xhmm_executable,
            cohort_id = cohort_id
    }

    output {
        Array[File] coverage_files = CollectCoverage.read_depth_file
        File merged_depth_file = MergeDepths.merged_depth_file
        File? extreme_gc_targets = GetExtremeGCContentTargets.extreme_gc_targets
        File? low_complexity_targets = GetLowComplexityTargets.low_complexity_targets
        File filtered_centered_counts = FilterTargetsAndMeanCenter.filtered_centered_counts
        File mean_center_step_excluded_targets = FilterTargetsAndMeanCenter.mean_center_step_excluded_targets
        File mean_center_step_excluded_samples = FilterTargetsAndMeanCenter.mean_center_step_excluded_samples
        File panel_of_normals = RunPCA.panel_of_normals
        File panel_of_normals_sd = RunPCA.panel_of_normals_sd
        File normalized_counts  = NormalizeCounts.normalized_counts
        File normalized_filtered_zcores = FilterAndZScoreCenter.normalized_filtered_zcores
        File z_score_step_excluded_targets = FilterAndZScoreCenter.z_score_step_excluded_targets
        File z_score_step_excluded_samples = FilterAndZScoreCenter.z_score_step_excluded_samples
        File original_depth_filtered = FilterOriginalReadDepth.original_depth_filtered
        File discovered_file = DiscoverCNV.discovered_file
        File genotyped_vcf_file = GenotypeCNV.genotyped_vcf_file
    }
}

task CramToBam {
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File cram
    File cram_idx
    
    String samtools_docker_image
    Int? mem_size
    Int? disk_size
    Int? preemptible_attempts

    String sample_id = basename(cram, ".cram")
    command <<<

        set -e
        set -o pipefail

        ln -vs ${ref_fasta} reference.fasta
        ln -vs ${ref_fasta_fai} reference.fasta.fai
        ln -vs ${ref_fasta_dict} reference.dict

        samtools view -h -T reference.fasta ${cram} |
        samtools view -b -o ${sample_id}.bam -
        samtools index -b ${sample_id}.bam
        mv ${sample_id}.bam.bai ${sample_id}.bai
    >>>
    runtime {
        docker: "${samtools_docker_image}"
        memory: select_first([mem_size, 4]) + " GB"
        disks: "local-disk " + select_first([disk_size, 200]) + " HDD"
        preemptible: select_first([preemptible_attempts, 5])
    }
    output {
        File output_bam = "${sample_id}.bam"
        File output_bam_index = "${sample_id}.bai"
    }
}

task CollectCoverage {
    File bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File interval_list

    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    String gatk_local_jar
    String base_filename = basename(bam, ".bam")
 
    command <<<
        set -e
        java -jar -Xmx${command_mem_mb}m ${gatk_local_jar} \
            -T DepthOfCoverage \
            -R ${ref_fasta} \
            -L ${interval_list} \
            -I ${bam} \
            -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
            --minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
            --includeRefNSites \
            --countType COUNT_FRAGMENTS \
            -o ${base_filename}
    >>>

    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File read_depth_file = "${base_filename}.sample_interval_summary"
    }
}

# Combines GATK Depth-of-Coverage outputs for multiple samples (at same loci)
task MergeDepths {
    Array[File] depth_files
    String xhmm_executable
    String cohort_id

    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
 
    command <<<
        ${xhmm_executable} --mergeGATKdepths -o ${cohort_id}.RD.txt \
            --GATKdepths ${sep=" --GATKdepths " depth_files}
    >>>

    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File merged_depth_file = "${cohort_id}.RD.txt"
    }
}

# Run GATK to calculate the per-target GC content and create a list of the targets with extreme GC content:
task GetExtremeGCContentTargets {
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File interval_list
    Array[Float] gc_low_high_filter_params_final
    String cohort_id

    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 3]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    String gatk_local_jar = "/root/GenomeAnalysisTK-3.8-0.jar"

    command <<<
        java -Xmx${command_mem_mb}m -jar ${gatk_local_jar} \
            -T GCContentByInterval -L ${interval_list} \
            -R ${ref_fasta} \
            -o ./${cohort_id}.locus_GC.txt

        cat ./${cohort_id}.locus_GC.txt | awk '{if ($2 < ${gc_low_high_filter_params_final[0]} || $2 > ${gc_low_high_filter_params_final[1]}) print $1}' \
        > ./extreme_gc_targets.txt
    >>>

    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File extreme_gc_targets = "extreme_gc_targets.txt"
    }
}

# calculate the fraction of repeat-masked bases in each target and create a list of those to filter out
task GetLowComplexityTargets {
    File seq_db
    File interval_list
    Float interval_complexity_filter_threshold
    String cohort_id

    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    command <<<
        /root/statgen-xhmm-cc14e528d909/sources/scripts/interval_list_to_pseq_reg ${interval_list} > ./${interval_list}.targets.reg

        pseq . loc-load --locdb ${seq_db} --file ./${interval_list}.targets.reg --group targets --out ./EXOME.targets.LOCDB.loc-load

        pseq . loc-stats --locdb ${seq_db} --group targets --seqdb ./seqdb | \
        awk '{if (NR > 1) print $_}' | sort -k1 -g | awk '{print $10}' | paste ./${interval_list}.interval_list - | \
        awk '{print $1"\t"$2}' \
        > ./${cohort_id}.locus_complexity.txt

        cat ./${cohort_id}.locus_complexity.txt | awk '{if ($2 > ${interval_complexity_filter_threshold}) print $1}' \
        > ./low_complexity_targets.txt
    >>>

    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File low_complexity_targets = "low_complexity_targets.txt"
    }
}

# Filters samples and targets and then mean-centers the targets
task FilterTargetsAndMeanCenter {
    File merged_depth_file
    File? extreme_gc_targets
    File? low_complexity_targets
    String xhmm_executable
    String cohort_id

    #Filtering parameters
    Int min_target_size_filter
    Int max_target_size_filter
    Int min_mean_target_rd
    Int max_mean_target_rd
    Int min_mean_sample_rd
    Int max_mean_sample_rd
    Int max_sd_sample_rd

    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    command <<<
        ${xhmm_executable} --matrix -r ${merged_depth_file} --centerData --centerType target \
            -o ./${cohort_id}.filtered_centered.RD.txt \
            --outputExcludedTargets ./${cohort_id}.filtered_centered.RD.txt.filtered_targets.txt \
            --outputExcludedSamples ./${cohort_id}.filtered_centered.RD.txt.filtered_samples.txt \
            ${"--excludeTargets " + extreme_gc_targets} ${"--excludeTargets " + low_complexity_targets} \
            --minTargetSize ${min_target_size_filter} --maxTargetSize ${max_target_size_filter} \
            --minMeanTargetRD ${min_mean_target_rd} --maxMeanTargetRD ${max_mean_target_rd} \
            --minMeanSampleRD ${min_mean_sample_rd} --maxMeanSampleRD ${max_mean_sample_rd} \
            --maxSdSampleRD ${max_sd_sample_rd}
    >>>

    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File filtered_centered_counts = "${cohort_id}.filtered_centered.RD.txt"
        File mean_center_step_excluded_targets = "${cohort_id}.filtered_centered.RD.txt.filtered_targets.txt"
        File mean_center_step_excluded_samples = "${cohort_id}.filtered_centered.RD.txt.filtered_samples.txt"
    }
}

#Runs PCA on mean-centered data
task RunPCA {
    File filtered_centered_counts
    String xhmm_executable
    String cohort_id

    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    command <<<
        ${xhmm_executable} --PCA -r ${filtered_centered_counts} --PCAfiles ./${cohort_id}.RD_PCA
    >>>
    
    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File panel_of_normals = "${cohort_id}.RD_PCA.PC.txt"
        File panel_of_normals_sd = "${cohort_id}.RD_PCA.PC_SD.txt"
    }
}

# Normalizes mean-centered data using PCA information
task NormalizeCounts {
    File panel_of_normals
    File panel_of_normals_sd
    File filtered_centered_counts
    String xhmm_executable
    String cohort_id
 
    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    command <<<
        ${xhmm_executable} --normalize -r ${filtered_centered_counts} --PCAfiles $(dirname "${panel_of_normals}")/${cohort_id}.RD_PCA \
            --normalizeOutput ./${cohort_id}.PCA_normalized.txt \
            --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
    >>>

    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }
    
    output {
        File normalized_counts  = "${cohort_id}.PCA_normalized.txt"
    }
}

# Filters and z-score centers (by sample) the PCA-normalized data
task FilterAndZScoreCenter {
    File normalized_counts
    String xhmm_executable
    String cohort_id

    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    command <<<
        ${xhmm_executable} --matrix -r ${normalized_counts} --centerData --centerType sample --zScoreData \
            -o ./${cohort_id}.PCA_normalized.filtered.sample_zscores.RD.txt \
            --outputExcludedTargets ./${cohort_id}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
            --outputExcludedSamples ./${cohort_id}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
            --maxSdTargetRD 50
    >>>
    
    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File normalized_filtered_zcores = "${cohort_id}.PCA_normalized.filtered.sample_zscores.RD.txt"
        File z_score_step_excluded_targets = "${cohort_id}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt"
        File z_score_step_excluded_samples = "${cohort_id}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt"
    }
}

# Filters original read-depth data to be the same as filtered, normalized data
task FilterOriginalReadDepth {
    File merged_depth_file
    File centered_filtered_targets
    File pca_normalized_filtered_targets
    File centered_filtered_samples
    File pca_normalized_filtered_samples
    String xhmm_executable
    String cohort_id

    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    command <<<
    ${xhmm_executable} --matrix -r ${merged_depth_file} \
        --excludeTargets ${centered_filtered_targets} \
        --excludeTargets ${pca_normalized_filtered_targets} \
        --excludeSamples ${centered_filtered_samples} \
        --excludeSamples ${pca_normalized_filtered_samples} \
        -o ./${cohort_id}.same_filtered.RD.txt
    >>>

    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File original_depth_filtered = "${cohort_id}.same_filtered.RD.txt"
    }
}

# Discovers CNVs in normalized data
task DiscoverCNV {
    File normalized_filtered_zcores
    File original_depth_filtered
    File discover_parameters_file
    String xhmm_executable
    String cohort_id

    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    command <<<
        ${xhmm_executable} --discover -p ${discover_parameters_file} \
            -r ${normalized_filtered_zcores} -R ${original_depth_filtered} \
            -c ${cohort_id}.xcnv -a ./${cohort_id}.aux_xcnv -s ./${cohort_id}
    >>>

    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File discovered_file = "${cohort_id}.xcnv"
    }
}

# Genotypes discovered CNVs in all samples
task GenotypeCNV {
    File normalized_filtered_zcores
    File original_depth_filtered
    File discover_file
    File discover_parameters_file
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String xhmm_executable
    String cohort_id

    #Runtime parameters
    String xhmm_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    command <<<
        ${xhmm_executable} --genotype -p ${discover_parameters_file} \
            -r ${normalized_filtered_zcores} -R ${original_depth_filtered} \
            -g ${discover_file} -F ${ref_fasta} \
            -v ./${cohort_id}.vcf
    >>>

    runtime {
        docker: "${xhmm_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File genotyped_vcf_file = "${cohort_id}.vcf"
    }   
}