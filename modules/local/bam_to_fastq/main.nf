process BAM_TO_FASTQ {

    /*

    */
    
    cpus 4
    memory "16GB"
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_R1.fastq"), path("${meta.sample_name}_${meta.molecule}_R2.fastq")


    script:
        """
        samtools fastq --threads $task.cpus -1 ${meta.sample_name}_${meta.molecule}_R1.fastq -2 ${meta.sample_name}_${meta.molecule}_R2.fastq -n $bam
        """

}
