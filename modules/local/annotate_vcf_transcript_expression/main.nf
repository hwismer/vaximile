process ANNOTATE_VCF_TRANSCRIPT_EXPRESSION {

    /*

    Use the abundance estimates from kallist to annotate transcript expression in a vcf file.

    */
    label 'process_low'

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_name), val(somatic_meta), path(vcf), val(sample_meta), path(tx_abundance)

    output:
        tuple val(somatic_meta), path("*_tx_expression.vcf")

    script:
        def prefix = task.ext.prefix ?: "${somatic_name}"
        """
        vcf-expression-annotator \
            $vcf \
            -s ${sample_meta.sample_name} \
            $tx_abundance \
            kallisto transcript \
            -o ${prefix}_tx_expression.vcf
        """

}
