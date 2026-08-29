process STAR_SORT_INDEX_BAM {

    /*

        Index the BAM file from a star process.

    */

    label 'process_high'
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*_STAR_sorted.bam"), path("*_STAR_sorted.bam.bai"), emit: bam
        path "versions.yml", topic: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        samtools sort --threads $task.cpus  $bam -o "${prefix}_STAR_sorted.bam"
        samtools index -@ $task.cpus  "${prefix}_STAR_sorted.bam"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        touch ${prefix}_STAR_sorted.bam
        touch ${prefix}_STAR_sorted.bam.bai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: 1.23.1
        END_VERSIONS
        """

}
