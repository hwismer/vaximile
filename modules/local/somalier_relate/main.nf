process SOMALIER_RELATE {


    // somalier has no internal threading - its README parallelises across samples,
    // which Nextflow already does. nf-core's somalier modules are process_low.
    label 'process_low'
    conda "bioconda::somalier=0.3.2-0 bioconda::htslib=1.23.1"

    tag "somalier relate on patient ${patient}"

    input:
        tuple val(patient), val(metas),  path(files)

    output:
        tuple val(patient), path("*.pairs.tsv"),  emit: pairs
        tuple val(patient), path("*.samples.tsv"), emit: samples
        tuple val(patient), path("*.groups.tsv"), emit: groups
        tuple val(patient), path("*.html"), emit:html
        path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: "${patient}"
    """
    export SOMALIER_REPORT_ALL_PAIRS=1
    somalier relate \
        -o $prefix \
        *.somalier
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(somalier --version 2>&1 | grep -Eo '[0-9]+\.[0-9.]+' | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${patient}"
    """
    touch ${prefix}.pairs.tsv
    touch ${prefix}.samples.tsv
    touch ${prefix}.groups.tsv
    touch ${prefix}.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: 0.3.2
    END_VERSIONS
    """
}
