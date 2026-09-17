process MHCFLOW {

    label 'process_max'

    conda "${moduleDir}/environment.yml"

    input:
        tuple val(meta), val(sample_name), path(bam), path(bai)
        // The novoalign index comes in with the FASTA rather than being passed to mhcflow:
        // mhcflow derives it from --ref as `ref.with_suffix(".nix")` and never takes it as
        // an argument, so it only has to be staged next to the FASTA. Built by
        // NOVOALIGN_HLA_FASTA; without it mhcflow exits on
        // "Failed to find HLA reference novoalign index file".
        tuple path(hla_fasta), path(hla_fai), path(hla_nix)
        path hla_bed
        path hla_kmers
        path hla_freqs

    output:
        // val(meta), not path(meta): meta is the metadata map, not a file. Declaring it
        // as a path makes Nextflow look for a file literally named by the map's toString,
        // which fails with "Missing output file(s) [somatic_name:..., patient:...]".
        tuple val(meta), path("${meta.sample_name}"), emit: out
        // finalizer/ holds what an LOH analysis needs, and is not the same as realigner/.
        // mhcflow realigns twice: once against the full HLA reference under realigner/, and
        // again under finalizer/ against only the alleles it typed for this sample. LOH
        // compares tumour and normal over one subject-specific reference, so it is the
        // finalizer copy that matters here. The FASTA and its .nix travel together because
        // the tumour's realignment is driven from this sample's reference, and mhcflow
        // derives the index from the FASTA path rather than taking it as an argument.
        tuple val(meta), path("${meta.sample_name}/finalizer/*.hla.fasta"), path("${meta.sample_name}/finalizer/*.hla.nix"), emit: sample_hla_ref
        tuple val(meta), path("${meta.sample_name}/finalizer/*.hla.realn.bam"), path("${meta.sample_name}/finalizer/*.hla.realn.bam.bai"), emit: realn_bam
        path "versions.yml", topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    mhcflow --bam $bam \
        --ref $hla_fasta \
        --bed $hla_bed \
        --tag $hla_kmers \
        --freq $hla_freqs \
        --nproc $task.cpus \
        $args \
        --outdir ${meta.sample_name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unknown: unknown
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.sample_name}/finalizer
    touch ${meta.sample_name}/finalizer/${meta.sample_name}.hla.fasta
    touch ${meta.sample_name}/finalizer/${meta.sample_name}.hla.nix
    touch ${meta.sample_name}/finalizer/${meta.sample_name}.hla.realn.bam
    touch ${meta.sample_name}/finalizer/${meta.sample_name}.hla.realn.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unknown: unknown
    END_VERSIONS
    """
}
