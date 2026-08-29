process MARK_DUPLICATES_SPARK {

    /*
    Par of GATK pre-processing best practices. Takes an aligned bam or sam file and outputrdinate-sorted
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
        tuple val(meta), path("*_markdup.bam"), path("*_markdup.bam.bai")

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    touch ${prefix}_markdup.bam
    touch ${prefix}_markdup.bam.bai
    """

}
