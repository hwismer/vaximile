process ANNOTATE_VCF_GENE_EXPRESSION {

    /*

    Use the abundance estimates from kallist to annotate gene expression in a vcf file.

    */
    label 'process_low'

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_name), val(somatic_meta), val(sample_name), path(vcf), val(sample_meta), path(gene_abundance)


    output:
        tuple val(somatic_meta), path("*_gene_expression.vcf"), emit: vcf
        path "versions.yml", topic: versions

    script:
        def prefix = task.ext.prefix ?: "${somatic_name}"

        """
        vcf-expression-annotator \
            $vcf \
            $gene_abundance \
            custom gene \
            -i ENSEMBLID \
            -e TPM \
            -s ${sample_name} \
            --ignore-ensembl-id-version \
            -o "${prefix}_gene_expression.vcf"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vatools: 5.2.0
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_name}"

        """
        touch ${prefix}_gene_expression.vcf
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vatools: 5.2.0
        END_VERSIONS
        """

}
