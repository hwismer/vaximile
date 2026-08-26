process HAPLOTYPE_CALLER_SCATTER {

    /*

    Use HaplotypeCaller on a single scattered interval. Post processes with CNNScoreVariants.

    */

    cpus 4
    memory "16GB"
    container "broadinstitute/gatk:4.3.0.0"

    tag "Haplotype Caller on ${meta.sample_name} on ${interval_shard}"

    input:
        tuple val(meta), path(bam), path(bai), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict
        val(interval_padding)

    output:
        tuple val(meta), path("${meta.sample_name}_${interval_shard}.vcf.gz"), path("${meta.sample_name}_${interval_shard}.vcf.gz.tbi"),  path(interval_shard)

    script:
        """
        gatk HaplotypeCaller \
            -R $reference_fa \
            -I $bam \
            -L $interval_shard \
            -O "${meta.sample_name}_${interval_shard}.vcf.gz" \
            -ERC NONE \
            --native-pair-hmm-threads $task.cpus \
            -ip $interval_padding \
            --create-output-variant-index
        """
}
