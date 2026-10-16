process TUMOR_COMPARISON_QC {

    // Somatic variants shared between a patient's tumours.

    label 'process_low'
    conda "python=3.10 pandas=2.1"

    tag "Comparing tumours for ${patient}"

    input:
        tuple val(patient), val(pair_names), val(tumor_names), path(vcfs, stageAs: 'pair?/*')
        path(script)

    output:
        tuple val(patient), path("*_mqc.json"), emit: mqc

    script:
    """
    python3 ${script} compare \\
        --patient ${patient} \\
        --names ${pair_names.join(' ')} \\
        --tumors ${tumor_names.join(' ')} \\
        --vcfs ${vcfs}
    """

    stub:
    """
    touch ${patient}_somatic_shared_mqc.json
    """
}
