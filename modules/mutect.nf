process MUTECT2_FILTER_MUTECT_CALLS {
    
    cpus 4
    memory "16GB"
    
    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple val(somatic_meta), path(unfiltered_vcf), path(orientation_model), path(stats), path(contamination_table)
        tuple path(reference_fa), path(reference_index), path(reference_dict)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_mutect_filtered.vcf.gz"), path("${somatic_meta.somatic_name}_mutect_filtered.vcf.gz.tbi"), emit: filtered_vcf


    script:
    """
    gatk FilterMutectCalls \
        -R $reference_fa \
        -V $unfiltered_vcf \
        --orientation-bias-artifact-priors $orientation_model \
        --contamination-table $contamination_table \
        -stats $stats \
        --create-output-variant-index \
        -O "${somatic_meta.somatic_name}_mutect_filtered.vcf.gz"
    """
}

process MUTECT2_MERGE_STATS {
    
    cpus 4
    memory "8GB"
    
    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple val(somatic_meta), path(stats)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_merged.stats")

    script:
    
    def stat_as_input = stats.collect {stat ->
            "-stats ${stat}"
        }.join(' ')

    """
     gatk MergeMutectStats \
        $stat_as_input \
        -O "${somatic_meta.somatic_name}_merged.stats"
    """
}

process MUTECT2_LEARN_READ_ORIENTATION {
    
    cpus 4
    memory "8GB"
    
    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple val(somatic_meta), path(f1r2s)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_orientmodel.tar.gz")

    script:
    
    def f1r2_as_input = f1r2s.collect { f1r2 ->
            "-I ${f1r2}"
        }.join(' ')

    """
    gatk LearnReadOrientationModel \
        $f1r2_as_input \
        -O "${somatic_meta.somatic_name}_orientmodel.tar.gz"
    """
}

process MUTECT2_CALCULATE_CONTAMINATION {

    cpus 4
    memory "8GB"
    
    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple val(somatic_meta), path(tumor_pileups), path(normal_pileups)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_contamination.table")

    script:

    """
    gatk CalculateContamination \
        -I $tumor_pileups \
        -matched $normal_pileups \
        -O "${somatic_meta.somatic_name}_contamination.table"
    """

}



process MUTECT2_GATHER_VCFS {

    cpus 2
    memory "8GB"
    
    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple val(somatic_meta), path(vcfs)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_merged.vcf")

    script:
    
    def sorted_vcfs = vcfs.sort { a, b -> a.name <=> b.name }
    def vcf_as_input = sorted_vcfs.collect { vcf ->
            "--INPUT ${vcf}"
        }.join(' ')

    """
    gatk GatherVcfs \
        $vcf_as_input \
        -O "${somatic_meta.somatic_name}_merged.vcf"
    
    """

}

process MUTECT2_GATHER_SELECT_VARIANTS {

    cpus 2
    memory "8GB"
    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple val(somatic_meta), path(vcf), path(vcf_index), val(interval_index), path(interval_shard)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${interval_index}.vcf.gz")

    script:
    
    """
    gatk SelectVariants \
        -V $vcf \
        -L $interval_shard \
        -O "${somatic_meta.somatic_name}_${interval_index}.vcf.gz"
    """

}


process MUTECT2_SCATTER {


    /*

        Call mutect 2 on a single interval shard. This process creates a new metadata block specific to
        the somatic caller of the form [somatic_name, somatic_caller, tumor_metadata, normal_metadata].

    */

    cpus 4
    memory "12GB"
    cache "lenient"

    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple val(somatic_meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(interval_index), path(interval_shard)
        tuple path(reference_fa), path(reference_index), path(reference_dict)
        tuple path(germline_resource), path(germline_resource_index)
        tuple path(pon), path(pon_index)
        val interval_padding

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${interval_index}_mutect.vcf.gz"), 
        path("${somatic_meta.somatic_name}_${interval_index}_mutect.vcf.gz.tbi"), val(interval_index), path(interval_shard), emit: vcf
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${interval_index}_mutect_f1r2.tar.gz"), emit: f1r2
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${interval_index}_mutect.vcf.gz.stats"), emit: stats

    script:

        """
        gatk Mutect2 \
            -R "${reference_fa}" \
            -I ${tumor_bam} \
            -I ${normal_bam} \
            -normal ${somatic_meta.normal_meta.sample_name} \
            --germline-resource $germline_resource \
            --panel-of-normals $pon \
            --f1r2-tar-gz "${somatic_meta.somatic_name}_${interval_index}_mutect_f1r2.tar.gz" \
            -L $interval_shard \
            -ip $interval_padding \
            -O "${somatic_meta.somatic_name}_${interval_index}_mutect.vcf.gz" \
            --native-pair-hmm-threads $task.cpus
        """
}

