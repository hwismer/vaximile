process MARK_DUPLICATES_SPARK {

    /*
    Par of GATK pre-processing best practices. Takes an aligned bam or sam file and outputrdinate-sorted
    BAM file with duplicates marked.
    */

    cpus 16
    memory "32GB"
    container "broadinstitute/gatk:4.6.1.0"
    clusterOptions '--gres=scratch:600G'

    tag "MarkDuplicatesSpark on ${meta.sample_name}"

    input:
        tuple val(meta), path(aligned_sam)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_markdup.bam"), path("${meta.sample_name}_${meta.molecule}_markdup.bam.bai")
    
    script:

    """
    mkdir -p tmp
    gatk MarkDuplicatesSpark \
        -I $aligned_sam \
        -O "${meta.sample_name}_${meta.molecule}_markdup.bam" \
        --create-output-bam-index \
        --tmp-dir ./tmp \
        --spark-master local[${task.cpus}] \
        --conf spark.local.dir=./tmp \
        --conf spark.sql.shuffle.partitions=${task.cpus * 3} \
        --conf spark.executor.memory=${(task.memory.toGiga() * 0.8) as int}g \
        --conf spark.driver.memory=8g
    """

}
