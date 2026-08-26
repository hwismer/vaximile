process ANNOTATE_VCF_GENE_EXPRESSION {

    /*

    Use the abundance estimates from kallist to annotate transcript expression in a vcf file.

    */
    cpus 2
    memory "16GB"

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_name), val(somatic_meta), path(vcf), val(sample_meta), path(gene_abundance)


    output:
        tuple val(somatic_meta), path("${somatic_name}_gene_expression.vcf")

    script:

        """
        vcf-expression-annotator \
            $vcf \
            $gene_abundance \
            custom gene \
            -i ENSEMBLID \
            -e TPM \
            -s ${sample_meta.sample_name} \
            --ignore-ensembl-id-version \
            -o "${somatic_name}_gene_expression.vcf"
        """

}
