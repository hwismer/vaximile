process MUTECT2_MERGE_STATS {
    
    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"
    
    tag "Merging Mutect stats for ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(stats)

    output:
        tuple val(somatic_meta), path("*_merged.stats"), emit: stats

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

    stub:
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
    """
    touch ${prefix}_merged.stats
    """
}
