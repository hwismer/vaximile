process NOVOALIGN_HLA_FASTA {

    // process_single: novoindex is single-threaded and the HLA reference is a few hundred
    // kilobases of allele sequence, so this finishes in seconds.
    label 'process_single'

    conda "bioconda::novoalign=3.09.04"

    input:
        tuple path(hla_fasta), path(hla_fai)

    output:
        // The .nix is emitted alongside the FASTA it was built from, because MHCFLOW has
        // to stage all three into one directory: mhcflow does not take the index as an
        // argument, it derives the path from --ref with Path.with_suffix(".nix") and then
        // calls `novoalign -d <that path>`. The index must therefore sit next to the FASTA
        // under exactly the same stem.
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
