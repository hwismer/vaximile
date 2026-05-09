process POSTPROCESS_GERMLINE {

    /*

        Germline variant post-processing similar to VT_SOMATIC_POSTPROCESS using
        variant normalization, decomposition, sort, and duplicate removal and indexing.

    */
   
    cpus 2
    memory "8GB"
    container "staphb/bcftools:1.23.1" 

    tag "Normalizing and filtering germline calls on ${meta.sample_name}"

    input:
        tuple val(meta), path(germline_vcf), path(germline_vcf_index)
        tuple path(reference_fa), path(reference_index), path(reference_dict)

    output:
        tuple val(meta), path("${meta.sample_name}_germline.vcf.gz"), path("${meta.sample_name}_germline.vcf.gz.tbi"), emit: germline_vcf
    script:
        """
        bcftools view -f PASS $germline_vcf -Oz -o filtered_germline.vcf.gz
        bcftools norm -m -any -f $reference_fa filtered_germline.vcf.gz -Oz -o ${meta.sample_name}_germline.vcf.gz
        bcftools index -t --threads $task.cpus ${meta.sample_name}_germline.vcf.gz
     """
}

process HAPLOTYPE_CALLER_FILTER_VARIANTS {
    
    cpus 4
    memory "32GB"

    container "broadinstitute/gatk:4.3.0.0"

    tag "Filtering germline variants in ${meta.sample_name}"

    input:
        tuple val(meta), path(vcf)
        tuple path(hapmap), path(hapmap_index)
        tuple path(mills), path(mills_index)

    output:
        tuple val(meta), path("${meta.sample_name}_germline_filtered.vcf.gz"), path("${meta.sample_name}_germline_filtered.vcf.gz.tbi"), emit: germline_vcf
    
    script:
    """
    gatk FilterVariantTranches \
        -V $vcf \
        --resource $hapmap \
        --resource $mills \
        --info-key CNN_1D \
        --snp-tranche 99.95 \
        --indel-tranche 99.4 \
        -O ${meta.sample_name}_germline_filtered.vcf.gz \
        --create-output-variant-index
    """
}


process HAPLOTYPE_CALLER_GATHER_VCFS {

    cpus 2
    memory "8GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Gathering VCFs for ${sample_meta.sample_name}"

    input:
        tuple val(sample_meta), path(vcfs)

    output:
        tuple val(sample_meta), path("${sample_meta.sample_name}_merged.vcf")

    script:
    
    def sorted_vcfs = vcfs.sort { a, b -> a.name <=> b.name }
    def vcf_as_input = sorted_vcfs.collect { vcf ->
            "--INPUT ${vcf}"
        }.join(' ')

    """
    gatk GatherVcfs \
        $vcf_as_input \
        -O "${sample_meta.sample_name}_merged.vcf"
    
    """

}

process HAPLOTYPE_CALLER_GATHER_SELECT_VARIANTS {

    cpus 2
    memory "8GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Select variants from ${sample_meta.sample_name} within ${interval_index}"

    input:
        tuple val(sample_meta), path(vcf), path(vcf_index), val(interval_index), path(interval_shard)

    output:
        tuple val(sample_meta), path("${sample_meta.sample_name}_${interval_index}.vcf.gz")

    script:
    
    """
    gatk SelectVariants \
        -V $vcf \
        -L $interval_shard \
        -O "${sample_meta.sample_name}_${interval_index}.vcf.gz"
    """
}


process HAPLOTYPE_CALLER_SCATTER {

    /*

    Use HaplotypeCaller on a single scattered interval. Post processes with CNNScoreVariants.

    */

    cpus 4
    memory "16GB"

    container "broadinstitute/gatk:4.3.0.0"

    tag "Calling germline variants on ${meta.sample_name} on ${interval_index}"

    input:
        tuple val(meta), path(bam), path(bai), val(interval_index), path(interval_shard)
        tuple path(reference_fa), path(reference_index), path(reference_dict)
        val(interval_padding)

    output:
        tuple val(meta), path("${meta.sample_name}_${interval_index}.vcf.gz"), path("${meta.sample_name}_${interval_index}.vcf.gz.tbi"), val(interval_index), path(interval_shard)

    script:
        """
        gatk HaplotypeCaller \
            -R $reference_fa \
            -I $bam \
            -L $interval_shard \
            -O "${meta.sample_name}_${interval_index}.vcf.gz" \
            -ERC NONE \
            --native-pair-hmm-threads $task.cpus \
            -ip $interval_padding \
            --create-output-variant-index
        """
}

process HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS {
    
    /*

    Use HaplotypeCaller on a single scattered interval. Post processes with CNNScoreVariants.

    */

    cpus 4
    memory "16GB"

    container "broadinstitute/gatk:4.3.0.0"

    tag "Scoring variants from ${meta.sample_name} on ${interval_index}"

    input:
        tuple val(meta), path(vcf), path(vcf_index),val(interval_index), path(interval_shard)
        tuple path(reference_fa), path(reference_index), path(reference_dict)
        val(interval_padding)

    output:
        tuple val(meta), path("${meta.sample_name}_${interval_index}_CNN.vcf.gz"), path("${meta.sample_name}_${interval_index}_CNN.vcf.gz.tbi"), val(interval_index), path(interval_shard)

    script:
        """
        gatk CNNScoreVariants \
            -V $vcf \
            -L $interval_shard \
            -ip $interval_padding \
            -R $reference_fa \
            --create-output-variant-index \
            -O "${meta.sample_name}_${interval_index}_CNN.vcf.gz"
        """
}

