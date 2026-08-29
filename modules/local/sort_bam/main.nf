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
        tuple val(meta), path("*_sorted.bam"), path("*_sorted.bam.bai")

    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        samtools sort --threads $task.cpus  $bam -o "${prefix}_sorted.bam"
        samtools index -@ $task.cpus  "${prefix}_sorted.bam"
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        touch ${prefix}_sorted.bam
        touch ${prefix}_sorted.bam.bai
        """

}
