process SORT_BAM {

    /*
    Sorts a BAM file and indexes it.
    */

    cpus 8
    memory "24GB"

    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    tag "Sorting ${meta.sample_name}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_sorted.bam.bai")

    script:
        """
        samtools sort --threads $task.cpus  $bam -o "${meta.sample_name}_${meta.molecule}_sorted.bam"
        samtools index -@ $task.cpus  "${meta.sample_name}_${meta.molecule}_sorted.bam"
        """

}
