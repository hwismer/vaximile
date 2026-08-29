process BAM_TO_FASTQ {

    /*

    */
    
    label 'process_medium'
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*_R1.fastq"), path("*_R2.fastq"), emit: reads
        path "versions.yml", topic: versions


    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        samtools fastq --threads $task.cpus -1 ${prefix}_R1.fastq -2 ${prefix}_R2.fastq -n $bam
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        touch ${prefix}_R1.fastq
        touch ${prefix}_R2.fastq
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: 1.23.1
        END_VERSIONS
        """

}
