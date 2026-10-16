process GERMLINE_QC {

    // VAF, het/hom, Ti/Tv and caller overlap plots for one normal.

    label 'process_low'
    conda "python=3.10 pandas=2.1 matplotlib=3.8"

    tag "Germline QC for ${meta.sample_name}"

    input:
        tuple val(meta), path(caller_vcfs), path(consensus_vcf, stageAs: 'consensus/*')
        path(script)

    output:
        tuple val(meta), path("*_mqc.{json,png}"), emit: mqc

    script:
    """
    python3 ${script} germline \\
        --sample ${meta.sample_name}_${meta.molecule} \\
        --sample-column ${meta.sample_name} \\
        --callers ${caller_vcfs} \\
        --consensus ${consensus_vcf}
    """

    stub:
    """
    touch ${meta.sample_name}_${meta.molecule}_germline_metrics_mqc.json ${meta.sample_name}_${meta.molecule}_germline_vaf_mqc.png
    """
}
