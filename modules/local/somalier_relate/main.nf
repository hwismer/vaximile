process SOMALIER_RELATE {


    label 'process_medium'
    conda "bioconda::somalier=0.3.2-0 bioconda::htslib=1.23.1"

    tag "somalier relate on patient ${patient}"

    input:
        tuple val(patient), val(metas),  path(files)

    output:
        tuple val(patient), path("*.pairs.tsv"),  emit: pairs
        tuple val(patient), path("*.samples.tsv"), emit: samples
        tuple val(patient), path("*.groups.tsv"), emit: groups
        tuple val(patient), path("*.html"), emit:html

    script:
    def prefix = task.ext.prefix ?: "${patient}"
    """
    export SOMALIER_REPORT_ALL_PAIRS=1
    somalier relate \
        -o $prefix \
        *.somalier
    """

    stub:
    def prefix = task.ext.prefix ?: "${patient}"
    """
    touch ${prefix}.pairs.tsv
    touch ${prefix}.samples.tsv
    touch ${prefix}.groups.tsv
    touch ${prefix}.html
    """
}
