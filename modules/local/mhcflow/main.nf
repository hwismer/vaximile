process MHCFLOW {

    label 'process_max'

    conda "${moduleDir}/environment.yml"

    input:
        tuple val(meta), val(sample_name), path(bam), path(bai)
        tuple path(hla_fasta), path(hla_fai)
        path hla_bed
        path hla_kmers
        path hla_freqs

    output:
        // val(meta), not path(meta): meta is the metadata map, not a file. Declaring it
        // as a path makes Nextflow look for a file literally named by the map's toString,
        // which fails with "Missing output file(s) [somatic_name:..., patient:...]".
        tuple val(meta), path("${meta.sample_name}"), emit: out
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
    mkdir -p ${meta.sample_name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unknown: unknown
    END_VERSIONS
    """
}
