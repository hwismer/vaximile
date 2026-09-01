process NOVOALIGN_HLA_FASTA {

    // process_single: the script body is empty - this only re-emits its inputs and
    // writes versions.yml - so process_max was reserving 32 CPUs to do nothing.
    label 'process_single'

    conda "bioconda::novoalign=4.03.04"

    input:
        tuple path(hla_fasta), path(hla_fai)

    output:
        tuple path(hla_fasta), path(hla_fai), emit: out
        path "versions.yml", topic: versions

    script:
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        novoalign: 4.03.04
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        novoalign: 4.03.04
    END_VERSIONS
    """
}
