process ANNOTATE_VCF_TRANSCRIPT_EXPRESSION {

    /*

    Annotate transcript expression in a VCF from salmon's quant.sf.

    salmon has no native parser in vcf-expression-annotator, so this goes through `custom`:
    quant.sf is Name/Length/EffectiveLength/TPM/NumReads where the kallisto parser expects
    target_id/.../tpm. --ignore-ensembl-id-version is needed because the Ensembl cDNA FASTA
    salmon was run against carries versioned transcript IDs (ENST...\.N).

    */
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
