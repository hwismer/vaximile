process MUTECT2_MERGE_STATS {
    
    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"
    
    tag "Merging Mutect stats for ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(stats)

    output:
        tuple val(somatic_meta), path("*_merged.stats"), emit: stats
        path "versions.yml", topic: versions

    script:

    def stat_as_input = stats.collect {stat ->
            "-stats ${stat}"
        }.join(' ')
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"

    """
     gatk MergeMutectStats \
        $stat_as_input \
        -O "${prefix}_merged.stats"
     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
     END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
    """
    touch ${prefix}_merged.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
