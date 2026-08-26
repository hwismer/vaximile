process FASTP {

    cpus 8
    memory "32GB"
    conda "bioconda::fastp=1.0.1"

    tag "FastP on ${meta.sample_name} w/ ${meta.molecule}"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_R1_fastp.fastq.gz"), path("${meta.sample_name}_${meta.molecule}_R2_fastp.fastq.gz"), emit: fastqs
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}.json"), emit: reports
    

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
              --detect_adapter_for_pe

        """
}
