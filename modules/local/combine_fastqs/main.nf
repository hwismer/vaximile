process COMBINE_FASTQS {
        
    // Combines FASTQS

    cpus 2
    memory "8GB"
    
    tag "Combing FASTQs from ${meta.sample_name}"

    input:
        tuple val(meta), path(fastqs_r1), path(fastqs_r2)
    output:
        tuple val(meta), path("${meta.somatic_name}_R1.merged.fastq.gz"), path("${meta.somatic_name}_R2.merged.fastq.gz")


    script:
    """
    cat ${fastqs_r1.join(' ')} > ${meta.somatic_name}_R1.merged.fastq.gz
    cat ${fastqs_r2.join(' ')} > ${meta.somatic_name}_R2.merged.fastq.gz
    """


}

//***************************************************************************************************************************
// INTERVAL PROCESSING
