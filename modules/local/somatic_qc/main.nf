process SOMATIC_QC {

    // VAF, caller overlap, filtering funnel, spectrum, TMB and pVACseq tier plots for one pair.

    label 'process_low'
    conda "python=3.10 pandas=2.1"

    tag "Somatic QC for ${meta.somatic_name}"

    input:
        tuple val(meta), path(caller_vcfs), path(consensus_vcf, stageAs: 'consensus/*'),
            path(final_vcf, stageAs: 'final/*'), path(contamination), path(capture_bed), path(pvacseq_report)
        path(script)

    output:
        tuple val(meta), path("*_mqc.json"), emit: mqc

    script:
    def pvacseq = pvacseq_report ? "--pvacseq ${pvacseq_report}" : ''
    """
    python3 ${script} somatic \\
        --pair ${meta.somatic_name} \\
        --tumor ${meta.tumor_meta.sample_name} \\
        --normal ${meta.normal_meta.sample_name} \\
        --callers ${caller_vcfs} \\
        --consensus ${consensus_vcf} \\
        --final ${final_vcf} \\
        --contamination ${contamination} \\
        --bed ${capture_bed} \\
        ${pvacseq}
    """

    stub:
    """
    touch ${meta.somatic_name}_somatic_metrics_mqc.json
    """
}
