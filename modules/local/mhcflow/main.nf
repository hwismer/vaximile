process MHCFLOW {

    cpus 32
    memory "64GB"

    conda "${moduleDir}/environment.yml"

    input:
        tuple val(meta), path(bam), path(bai)
        tuple path(hla_fasta), path(hla_fai)
        path hla_bed
        path hla_kmers
        path hla_freqs

    output:
        tuple path(meta), path("${meta.sample_name}")

    script:
    """
    mhcflow --bam $bam \
        --ref $hla_fasta \
        --bed $hla_bed \
        --tag $hla_kmers \
        --freq $hla_freqs \
        --nproc $task.cpus \
        --min-ecnt 1 \
        --outdir ${meta.sample_name}
	ls
    """
}
