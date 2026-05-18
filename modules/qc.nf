process FASTP {

    /*

    Run fastp to do initial qc on fastq files.

    Performs automatic adapter trimming and generates qc report.

    Outputs trimmed FASTQ files.


    */
    
    cpus 8
    memory "32GB"
    conda "bioconda::fastp=1.0.1"

    tag "FastP on ${meta.sample_name} w/ ${meta.molecule}"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_R1_fastp.fastq.gz"), path("${meta.sample_name}_${meta.molecule}_R2_fastp.fastq.gz"), emit: fastqs
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}.json"), emit: reports
    

    //publishDir "${params.outdir}/${meta.somatic_sample}/fastp_qc/${meta.sample_name}_${meta.molecule}_fastp/", mode: 'copy'

    script:
        """
        fastp --thread $task.cpus \
              -i ${reads[0]} \
              -I ${reads[1]} \
              -o "${meta.sample_name}_${meta.molecule}_R1_fastp.fastq.gz" \
              -O "${meta.sample_name}_${meta.molecule}_R2_fastp.fastq.gz" \
              -R "${meta.sample_name}_${meta.molecule}_fastp_report" \
              -h "${meta.sample_name}_${meta.molecule}.html" \
              -j "${meta.sample_name}_${meta.molecule}.json" \

        """
}

process SAMTOOLS_FLAGSTAT {

    cpus 8
    memory "24GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"
    
    publishDir "./flagstat"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("${meta.sample_name}_${meta.molecule}.flagstat")

    script:
    """
    samtools flagstat --threads $task.cpus $bam > ${meta.sample_name}_${meta.molecule}.flagstat
    """

}

process SAMTOOLS_COVERAGE {
    
    cpus 8
    memory "24GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("${meta.sample_name}_${meta.molecule}_coverage.tsv")

    script:
    """
    samtools coverage $bam > ${meta.sample_name}_${meta.molecule}_coverage.tsv
    """
}


process SAMTOOLS_IDXSTATS {
    
    cpus 4
    memory "16GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("${meta.sample_name}_${meta.molecule}_idxstats.tsv")

    script:
    """
    samtools idxstats $bam > ${meta.sample_name}_${meta.molecule}_idxstats.tsv
    """
}



