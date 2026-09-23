process HAPLOTYPE_CALLER_SCATTER {

    // HaplotypeCaller on one interval shard.

    label 'process_low'
    // GATK pinned to 4.3.0.0: CNNScoreVariants is not in newer GATK4. Do not bump.
    conda "bioconda::gatk4=4.3.0.0"
    container "broadinstitute/gatk:4.3.0.0"

    tag "Haplotype Caller on ${meta.sample_name} on ${interval_shard}"

    input:
        tuple val(meta), path(bam), path(bai), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict
        val(interval_padding)

    output:
        tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"),  path(interval_shard), emit: vcf

    script:
        def args = task.ext.args ?: ''
    // Size the JVM heap from task.memory so it stays within the reservation.
        def avail_mem = (task.memory.mega * 0.8).intValue()
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${interval_shard}"
        """
        gatk --java-options "-Xmx${avail_mem}M" HaplotypeCaller \
            -R $reference_fa \
            -I $bam \
            -L $interval_shard \
            -O "${prefix}.vcf.gz" \
            $args \
            --native-pair-hmm-threads $task.cpus \
            -ip $interval_padding \
            --create-output-variant-index
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${interval_shard}"
        """
        touch ${prefix}.vcf.gz
        touch ${prefix}.vcf.gz.tbi
        """
}
