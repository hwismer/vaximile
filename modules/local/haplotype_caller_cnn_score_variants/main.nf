process HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS {
    
    /*

    Use HaplotypeCaller on a single scattered interval. Post processes with CNNScoreVariants.

    */

    cpus 4
    memory "16GB"

    container "broadinstitute/gatk:4.3.0.0"

    tag "Scoring haplotypecaller variants from ${meta.sample_name} on ${interval_shard}"

    input:
        tuple val(meta), path(vcf), path(vcf_index), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict
        val(interval_padding)

    output:
        tuple val(meta), path("${meta.sample_name}_${interval_shard}_CNN.vcf.gz"), path("${meta.sample_name}_${interval_shard}_CNN.vcf.gz.tbi"), path(interval_shard)

    script:
        """
        gatk CNNScoreVariants \
            -V $vcf \
            -L $interval_shard \
            -ip $interval_padding \
            -R $reference_fa \
            --create-output-variant-index \
            -O "${meta.sample_name}_${interval_shard}_CNN.vcf.gz"
        """
}
