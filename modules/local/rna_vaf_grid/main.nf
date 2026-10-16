process RNA_VAF_GRID {

    // RNA vs DNA VAF scatter plots of every pair of a patient, side by side in one image.

    label 'process_low'
    conda "python=3.10 pandas=2.1 matplotlib=3.8"

    tag "RNA vs DNA VAF grid for ${patient}"

    input:
        tuple val(patient), path(tables)
        path(script)

    output:
        tuple val(patient), path("*_somatic_rna_vaf_mqc.png"), emit: png, optional: true

    script:
    """
    python3 ${script} rna-grid --patient ${patient} --tables ${tables}
    """

    stub:
    """
    touch ${patient}_somatic_rna_vaf_mqc.png
    """
}
