process MARK_DUPLICATES_SPARK {

    /*
    Part of GATK pre-processing best practices. Takes an aligned bam or sam file and outputs a coordinate-sorted
    BAM file with duplicates marked.
    */

    label 'process_very_high'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"
    clusterOptions '--gres=scratch:600G'

    tag "MarkDuplicatesSpark on ${meta.sample_name}"

    input:
        tuple val(meta), path(aligned_sam)

    output:
        tuple val(meta), path("*_markdup.bam"), path("*_markdup.bam.bai"), emit: bam
        path "versions.yml", topic: versions

    script:

    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    mkdir -p tmp
    gatk MarkDuplicatesSpark \
        -I $aligned_sam \
        -O "${prefix}_markdup.bam" \
        --create-output-bam-index \
        --tmp-dir ./tmp \
        --spark-master local[${task.cpus}] \
        --conf spark.local.dir=./tmp \
        --conf spark.sql.shuffle.partitions=${task.cpus * 3} \
        --conf spark.executor.memory=${(task.memory.toGiga() * 0.8) as int}g \
        --conf spark.driver.memory=8g
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    touch ${prefix}_markdup.bam
    touch ${prefix}_markdup.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """

}
