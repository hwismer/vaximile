process MHCFLOW {

    // Runs MHCflow in preparation for LOHHLA

    label 'process_max'

    conda "${moduleDir}/environment.yml"

    input:
        tuple val(meta), val(sample_name), path(bam), path(bai)
        tuple path(hla_fasta), path(hla_fai), path(hla_nix)
        path hla_bed
        path hla_kmers
        path hla_freqs

    output:
        tuple val(meta), path("${meta.sample_name}"), emit: out
        tuple val(meta), path("${meta.sample_name}/finalizer/*.hla.fasta"), path("${meta.sample_name}/finalizer/*.hla.nix"), emit: sample_hla_ref
        tuple val(meta), path("${meta.sample_name}/finalizer/*.hla.realn.bam"), path("${meta.sample_name}/finalizer/*.hla.realn.bam.bai"), emit: realn_bam

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
    """

    stub:
    """
    mkdir -p ${meta.sample_name}/finalizer
    touch ${meta.sample_name}/finalizer/${meta.sample_name}.hla.fasta
    touch ${meta.sample_name}/finalizer/${meta.sample_name}.hla.nix
    touch ${meta.sample_name}/finalizer/${meta.sample_name}.hla.realn.bam
    touch ${meta.sample_name}/finalizer/${meta.sample_name}.hla.realn.bam.bai
    """
}
