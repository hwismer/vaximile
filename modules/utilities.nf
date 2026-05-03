process COMBINE_FASTQS {

    cpus 2
    memory "8GB"

    input:
        tuple val(meta), path(fastqs_r1), path(fastqs_r2)
    output:
        tuple val(meta), path("${meta.somatic_name}_R1.merged.fastq.gz"), path("${meta.somatic_name}_R2.merged.fastq.gz")


    script:
    """
    cat ${fastqs_r1.join(' ')} > ${meta.somatic_name}_R1.merged.fastq.gz
    cat ${fastqs_r2.join(' ')} > ${meta.somatic_name}_R2.merged.fastq.gz
    """


}
process SPLIT_INTERVALS {

    /*

        Given a file of genomic intervals, split into scatter_count number of shards.

    */

    cpus 2
    memory "8GB"
    cache "lenient"

    tag "Split intervals for ${params.scatter_count} shards"

    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple path(reference_fa), path(reference_fa_index), path(reference_dict)
        path intervals_file
        val scatter_count

    output:
        path "*-scattered.interval_list", emit: interval_shards

    script:
    """
    gatk SplitIntervals \
        -R $reference_fa \
        -L $intervals_file \
        --scatter-count $scatter_count \
        -O .

    """
}


process SPLIT_INTERVALS_PADDED {

    /*

        Given a file of genomic intervals, split into scatter_count number of shards.

    */

    cpus 2
    memory "8GB"
    cache "lenient"

    tag "Split intervals for ${params.scatter_count} shards"

    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple path(reference_fa), path(reference_fa_index), path(reference_dict)
        path intervals_file
        val scatter_count

    output:
        path "*-scattered.interval_list", emit: interval_shards

    script:
    """
    gatk SplitIntervals \
        -R $reference_fa \
        -L $intervals_file \
        --scatter-count $scatter_count \
        --interval-padding 100 \
        -O .

    """
}


process BWA_POSTPROCESS {

    cpus 32
    memory "32GB"
    cache "lenient"

    conda "bioconda::samtools"

    //publishDir "${params.outdir}/${meta.somatic_sample}/alignment/bwa/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), val(bwa_sam)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_sorted.bam.bai"), emit: mapped_bam

    script:
    """

    samtools sort --threads $task.cpus $bwa_sam -o "${meta.sample_name}_${meta.molecule}_sorted.bam"
    samtools index "${meta.sample_name}_${meta.molecule}_sorted.bam"

    """
}
