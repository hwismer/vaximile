process NOVOALIGN_HLA_FASTA {

    cpus 32
    memory "32GB"

    conda "bioconda::novoalign=4.03.04"

    input:
        tuple path(hla_fasta), path(hla_fai)

    output:
        tuple path(hla_fasta), path(hla_fai)

    script:
    """

    """
}
