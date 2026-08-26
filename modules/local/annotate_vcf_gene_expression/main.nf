process ANNOTATE_VCF_GENE_EXPRESSION {

    /*

    Use the abundance estimates from kallist to annotate transcript expression in a vcf file.

    */
    label 'process_low'

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_name), val(somatic_meta), path(vcf), val(sample_meta), path(gene_abundance)


    output:
        tuple val(somatic_meta), path("*_gene_expression.vcf")

    script:
        def prefix = task.ext.prefix ?: "${somatic_name}"

        """
        vcf-expression-annotator \
            $vcf \
            $gene_abundance \
            custom gene \
            -i ENSEMBLID \
            -e TPM \
            -s ${sample_meta.sample_name} \
            --ignore-ensembl-id-version \
            -o "${prefix}_gene_expression.vcf"
        """

}
