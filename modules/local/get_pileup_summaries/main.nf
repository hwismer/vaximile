process GET_PILEUP_SUMMARIES {

    /*
    Takes a merged post-bqsr bam and a vcf of common germline sites and gets pileup summaries at provided sites
    */
    
    cpus 4
    memory "16GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "GetPileupSummaries on ${meta.sample_name}"

    input:
        tuple val(meta), path(bqsr_bam), path(bqsr_index)
        tuple path(common_germline), path(common_germline_index)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_pileups.table")
    
    script:
    """
    gatk GetPileupSummaries \
        -I $bqsr_bam \
        -V $common_germline \
        -L $common_germline \
        -O "${meta.sample_name}_${meta.molecule}_pileups.table"

    """
}
