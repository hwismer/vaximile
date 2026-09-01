process SORT_BAM {

    /*
    Sorts a BAM file and indexes it.
    */

    label 'process_high'

    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    tag "Sorting ${meta.sample_name}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*_sorted.bam"), path("*_sorted.bam.bai"), emit: bam
        path "versions.yml", topic: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        samtools sort --threads $task.cpus  $bam -o "${prefix}_sorted.bam"
        samtools index -@ $task.cpus  "${prefix}_sorted.bam"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        touch ${prefix}_sorted.bam
        touch ${prefix}_sorted.bam.bai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: 1.23.1
        END_VERSIONS
        """

}
