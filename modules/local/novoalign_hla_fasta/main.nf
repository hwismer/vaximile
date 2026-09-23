process NOVOALIGN_HLA_FASTA {

    // process_single: novoindex is single-threaded and the HLA reference is a few hundred
    // kilobases of allele sequence, so this finishes in seconds.
    label 'process_single'

    conda "bioconda::novoalign=3.09.04"

    input:
        tuple path(hla_fasta), path(hla_fai)

    output:
        // Emit the .nix with its FASTA: mhcflow expects the index next to the reference.
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
