process ANNOTATE_VCF_GENE_EXPRESSION {

    // Annotate a VCF with gene expression from salmon's quant.genes.sf.
    label 'process_low'

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_name), val(somatic_meta), val(sample_name), path(vcf), val(sample_meta), path(gene_abundance)


    output:
        tuple val(somatic_meta), path("*_gene_expression.vcf"), emit: vcf

    script:
        def prefix = task.ext.prefix ?: "${somatic_name}"

        """
        vcf-expression-annotator \
            $vcf \
            $gene_abundance \
            custom gene \
            -i Name \
            -e TPM \
            -s ${sample_name} \
            --ignore-ensembl-id-version \
            -o "${prefix}_gene_expression.vcf"
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_name}"

        """
        touch ${prefix}_gene_expression.vcf
        """

}
