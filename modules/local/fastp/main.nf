process FASTP {

    label 'process_high'
    conda "bioconda::fastp=1.0.1"

    tag "FastP on ${meta.sample_name} w/ ${meta.molecule}"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*_R1_fastp.fastq.gz"), path("*_R2_fastp.fastq.gz"), emit: fastqs
        tuple val(meta), path("*.json"), emit: reports


    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        fastp --thread $task.cpus \
              -i ${reads[0]} \
              -I ${reads[1]} \
              -o "${prefix}_R1_fastp.fastq.gz" \
              -O "${prefix}_R2_fastp.fastq.gz" \
              -R "${prefix}_fastp_report" \
              -h "${prefix}.html" \
              -j "${prefix}.json" \
              $args

        """
}
