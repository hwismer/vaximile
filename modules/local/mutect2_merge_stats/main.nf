process MUTECT2_MERGE_STATS {
    
    label 'process_low'
    container "broadinstitute/gatk:4.6.1.0"
    
    tag "Merging Mutect stats for ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(stats)

    output:
        tuple val(somatic_meta), path("*_merged.stats")

    script:

    def stat_as_input = stats.collect {stat ->
            "-stats ${stat}"
        }.join(' ')
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"

    """
     gatk MergeMutectStats \
        $stat_as_input \
        -O "${prefix}_merged.stats"
    """
}
