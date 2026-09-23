process NOVOALIGN_HLA_FASTA {

    label 'process_single'

    conda "bioconda::novoalign=3.09.04"

    input:
        tuple path(hla_fasta), path(hla_fai)

    output:
        tuple path(hla_fasta), path(hla_fai), path("${hla_fasta.baseName}.nix"), emit: out

    script:
    """
    novoindex "${hla_fasta.baseName}.nix" $hla_fasta

    """

    stub:
    """
    touch "${hla_fasta.baseName}.nix"

    """
}
