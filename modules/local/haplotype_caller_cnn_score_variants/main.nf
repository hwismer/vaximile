process HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS {
    
    /*

    Use HaplotypeCaller on a single scattered interval. Post processes with CNNScoreVariants.

    */

    label 'process_low_memory'

    // GATK PINNED TO 4.3.0.0 - DO NOT BUMP TO MATCH THE 4.6.1.0 MODULES.
    // This module is part of the CNNScoreVariants germline chain
    // (HaplotypeCaller -> CNNScoreVariants -> FilterVariantTranches). CNNScoreVariants
    // was deprecated in favour of NVScoreVariants and is not available in current GATK4,
    // so the three modules in this chain must stay on a release that still ships it.
    conda "bioconda::gatk4=4.3.0.0"
    container "broadinstitute/gatk:4.3.0.0"

    tag "Scoring haplotypecaller variants from ${meta.sample_name} on ${interval_shard}"

    input:
        tuple val(meta), path(vcf), path(vcf_index), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict
        val(interval_padding)

    output:
        tuple val(meta), path("*_CNN.vcf.gz"), path("*_CNN.vcf.gz.tbi"), path(interval_shard), emit: vcf
        path "versions.yml", topic: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${interval_shard}"
        """
        gatk CNNScoreVariants \
            -V $vcf \
            -L $interval_shard \
            -ip $interval_padding \
            -R $reference_fa \
            --create-output-variant-index \
            -O "${prefix}_CNN.vcf.gz"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${interval_shard}"
        """
        touch ${prefix}_CNN.vcf.gz
        touch ${prefix}_CNN.vcf.gz.tbi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: 4.3.0.0
        END_VERSIONS
        """
}
