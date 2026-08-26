process APPLY_BQSR_GATHER {

    /*
    Gathers all bqsr games to create a final merged bam with adjusted base quality scores.
    */

    label 'process_medium'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "GatherBams on ${meta.sample_name}"

    input:
        tuple val(meta), path(bams)
    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_bqsr.bam")

    script:
    def sorted_bams = bams.sort { it.name }
    """
    gatk GatherBamFiles \
        ${bams.collect { "-I ${it}" }.join(' ')} \
        -O "${meta.sample_name}_${meta.molecule}_bqsr.bam"
    """

}
