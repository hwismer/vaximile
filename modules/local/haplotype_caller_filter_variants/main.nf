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
        tuple val(meta), val("haplotypecaller"), path("${meta.sample_name}_germline_filtered.vcf.gz"), path("${meta.sample_name}_germline_filtered.vcf.gz.tbi"), emit: germline_vcf
    
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
