process ANNOTATE_VCF_TRANSCRIPT_EXPRESSION {

    /*

    Use the abundance estimates from kallist to annotate transcript expression in a vcf file.

    */
    cpus 2
    memory "16GB"

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_name), val(somatic_meta), path(vcf), val(sample_meta), path(tx_abundance)

    output:
        tuple val(somatic_meta), path("${somatic_name}_tx_expression.vcf")

    script:

        """
        vcf-expression-annotator \
            $vcf \
            -s ${sample_meta.sample_name} \
            $tx_abundance \
            kallisto transcript \
            -o ${somatic_name}_tx_expression.vcf
        """

}
