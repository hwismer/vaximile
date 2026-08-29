process HAPLOTYPE_CALLER_FILTER_VARIANTS {
    
    label 'process_medium'
    // GATK PINNED TO 4.3.0.0 - DO NOT BUMP TO MATCH THE 4.6.1.0 MODULES.
    // This module is part of the CNNScoreVariants germline chain
    // (HaplotypeCaller -> CNNScoreVariants -> FilterVariantTranches). CNNScoreVariants
    // was deprecated in favour of NVScoreVariants and is not available in current GATK4,
    // so the three modules in this chain must stay on a release that still ships it.
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
