process MUTECT2_MERGE_STATS {
    
    cpus 2
    memory "8GB"
    container "broadinstitute/gatk:4.6.1.0"
    
    tag "Merging Mutect stats for ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(stats)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_merged.stats")

    script:
    
    def stat_as_input = stats.collect {stat ->
            "-stats ${stat}"
        }.join(' ')

    """
     gatk MergeMutectStats \
        $stat_as_input \
        -O "${somatic_meta.somatic_name}_merged.stats"
    """
}
