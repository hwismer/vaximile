process BAM_TO_FASTQ {

    /*

    */
    
    label 'process_medium'
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*_R1.fastq"), path("*_R2.fastq")


    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        samtools fastq --threads $task.cpus -1 ${prefix}_R1.fastq -2 ${prefix}_R2.fastq -n $bam
        """

}
