process STAR_SORT_INDEX_BAM {

    // Sort and index the STAR BAM.

    label 'process_high'
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*_STAR_sorted.bam"), path("*_STAR_sorted.bam.bai"), emit: bam

    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        samtools sort --threads $task.cpus  $bam -o "${prefix}_STAR_sorted.bam"
        samtools index -@ $task.cpus  "${prefix}_STAR_sorted.bam"
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        touch ${prefix}_STAR_sorted.bam
        touch ${prefix}_STAR_sorted.bam.bai
        """

}
