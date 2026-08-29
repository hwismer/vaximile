process NOVOALIGN_HLA_FASTA {

    label 'process_max'

    conda "bioconda::novoalign=4.03.04"

    input:
        tuple path(hla_fasta), path(hla_fai)

    output:
        tuple path(hla_fasta), path(hla_fai)

    script:
    """

    """

    stub:
    """
    """
}
