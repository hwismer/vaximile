process GET_PILEUP_SUMMARIES {

    /*
    Takes a merged post-bqsr bam and a vcf of common germline sites and gets pileup summaries at provided sites
    */
    
    label 'process_medium'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "GetPileupSummaries on ${meta.sample_name}"

    input:
        tuple val(meta), path(bqsr_bam), path(bqsr_index)
        tuple path(common_germline), path(common_germline_index)

    output:
        tuple val(meta), path("*_pileups.table")

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    gatk GetPileupSummaries \
        -I $bqsr_bam \
        -V $common_germline \
        -L $common_germline \
        -O "${prefix}_pileups.table"

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    touch "${prefix}_pileups.table"
    """
}
