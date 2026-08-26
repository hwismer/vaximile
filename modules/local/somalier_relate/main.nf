process SOMALIER_RELATE {


    cpus 4
    memory "32GB"
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
    """
    export SOMALIER_REPORT_ALL_PAIRS=1
    somalier relate \
        -o $patient \
        *.somalier
    """
}
