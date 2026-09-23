process HAPLOTYPE_CALLER_FILTER_VARIANTS {
    
    label 'process_low_memory'
    // GATK pinned to 4.3.0.0: CNNScoreVariants is not in newer GATK4. Do not bump.
    conda "bioconda::gatk4=4.3.0.0"
    container "broadinstitute/gatk:4.3.0.0"

    tag "Filtering germline variants in ${meta.sample_name}"

    input:
        tuple val(meta), path(vcf)
        tuple path(hapmap), path(hapmap_index)
        tuple path(mills), path(mills_index)

    output:
        tuple val(meta), val("haplotypecaller"), path("*_germline_filtered.vcf.gz"), path("*_germline_filtered.vcf.gz.tbi"), emit: germline_vcf

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    gatk FilterVariantTranches \
        -V $vcf \
        --resource $hapmap \
        --resource $mills \
        $args \
        -O ${prefix}_germline_filtered.vcf.gz \
        --create-output-variant-index
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    touch ${prefix}_germline_filtered.vcf.gz
    touch ${prefix}_germline_filtered.vcf.gz.tbi
    """
}
