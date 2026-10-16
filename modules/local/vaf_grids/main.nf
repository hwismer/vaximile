process VAF_GRIDS {

    // A patient's VAF histograms and scatter plots, one panel per pair or normal, side by side.

    label 'process_low'
    conda "python=3.10 pandas=2.1 matplotlib=3.8"

    tag "VAF plots for ${patient}"

    input:
        tuple val(patient), path(tables)
        path(script)

    output:
        tuple val(patient), path("*_mqc.png"), emit: png, optional: true

    script:
    """
    python3 ${script} grid --patient ${patient} --tables ${tables}
    """

    stub:
    """
    touch ${patient}_somatic_vaf_mqc.png ${patient}_somatic_normal_vaf_mqc.png ${patient}_somatic_rna_vaf_mqc.png ${patient}_germline_vaf_mqc.png
    """
}
