process COMBINE_FASTQS {
        
    // Combines FASTQS

    label 'process_low'

    tag "Combing FASTQs from ${meta.sample_name}"

    input:
        tuple val(meta), path(fastqs_r1), path(fastqs_r2)
    output:
        tuple val(meta), path("*_R1.merged.fastq.gz"), path("*_R2.merged.fastq.gz")


    script:
    def prefix = task.ext.prefix ?: "${meta.somatic_name}"
    """
    cat ${fastqs_r1.join(' ')} > ${prefix}_R1.merged.fastq.gz
    cat ${fastqs_r2.join(' ')} > ${prefix}_R2.merged.fastq.gz
    """


}

//***************************************************************************************************************************
// INTERVAL PROCESSING
