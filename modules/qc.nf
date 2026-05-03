process FASTP {

    /*

    Run fastp to do initial qc on fastq files.

    Performs automatic adapter trimming and generates qc report.

    Outputs trimmed FASTQ files.


    */
    
    cpus 16
    memory "8GB"
    conda "bioconda::fastp=1.0.1"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_R1_fastp.fastq.gz"), path("${meta.sample_name}_${meta.molecule}_R2_fastp.fastq.gz"), emit: fastqs
        tuple val(meta), path("${meta.sample_name}_*{html,json}*"), emit: reports
    
    tag "FastP on ${meta.sample_name} w/ ${meta.molecule}"

    //publishDir "${params.outdir}/${meta.somatic_sample}/fastp_qc/${meta.sample_name}_${meta.molecule}_fastp/", mode: 'copy'

    script:
        """
        fastp --thread $task.cpus \
              -i ${reads[0]} \
              -I ${reads[1]} \
              -o "${meta.sample_name}_${meta.molecule}_R1_fastp.fastq.gz" \
              -O "${meta.sample_name}_${meta.molecule}_R2_fastp.fastq.gz" \
              -R "${meta.sample_name}_${meta.molecule}_fastp_report" \
              -h "${meta.sample_name}_${meta.molecule}_fastp_report.html" \
              -j "${meta.sample_name}_${meta.molecule}_fastp_report.json" \

        """
}
