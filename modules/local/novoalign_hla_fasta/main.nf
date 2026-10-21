process NOVOALIGN_HLA_FASTA {

    // Index the HLA reference for mhcflow with novoindex from --novoalign_dir.

    label 'process_single'

    input:
        tuple path(hla_fasta), path(hla_fai)
        path(novoalign_dir)

    output:
        tuple path(hla_fasta), path(hla_fai), path("${hla_fasta.baseName}.nix"), emit: out

    script:
    """
    ${novoalign_dir}/novoindex "${hla_fasta.baseName}.nix" $hla_fasta
    """

    stub:
    """
    touch "${hla_fasta.baseName}.nix"
    """
}
