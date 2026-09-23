process ANNOTATE_VCF_TRANSCRIPT_EXPRESSION {

    // Annotate a VCF with transcript expression from salmon's quant.sf, via the custom parser.
    label 'process_low'

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_name), val(somatic_meta), val(sample_name), path(vcf), val(sample_meta), path(tx_abundance)

    output:
        tuple val(somatic_meta), path("*_tx_expression.vcf"), emit: vcf

    script:
        def prefix = task.ext.prefix ?: "${somatic_name}"
        """
        vcf-expression-annotator \
            $vcf \
            -s ${sample_name} \
            $tx_abundance \
            custom transcript \
            -i Name \
            -e TPM \
            --ignore-ensembl-id-version \
            -o ${prefix}_tx_expression.vcf
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_name}"
        """
        touch ${prefix}_tx_expression.vcf
        """

}
