process MUTECT2_SCATTER {


    /*

        Call mutect 2 on a single interval shard. This process creates a new metadata block specific to
        the somatic caller of the form [somatic_name, somatic_caller, tumor_metadata, normal_metadata].

    */

    label 'process_medium'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Running Mutect2 scatter on ${somatic_meta.somatic_name} at $interval_shard"

    input:
        tuple val(somatic_meta), val(normal_sample_name), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict
        tuple path(germline_resource), path(germline_resource_index)
        tuple path(pon), path(pon_index)
        val interval_padding

    output:
        tuple val(somatic_meta), path("*_mutect.vcf.gz"),
        path("*_mutect.vcf.gz.tbi"), path(interval_shard), emit: vcf
        tuple val(somatic_meta), path("*_mutect_f1r2.tar.gz"), emit: f1r2
        tuple val(somatic_meta), path("*_mutect.vcf.gz.stats"), emit: stats

    script:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}_${interval_shard}"
        """
        gatk Mutect2 \
            -R "${reference_fa}" \
            -I ${tumor_bam} \
            -I ${normal_bam} \
            -normal ${normal_sample_name} \
            --germline-resource $germline_resource \
            --panel-of-normals $pon \
            --f1r2-tar-gz "${prefix}_mutect_f1r2.tar.gz" \
            -L $interval_shard \
            -ip $interval_padding \
            -O "${prefix}_mutect.vcf.gz" \
            --native-pair-hmm-threads $task.cpus
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}_${interval_shard}"
        """
        touch ${prefix}_mutect.vcf.gz
        touch ${prefix}_mutect.vcf.gz.tbi
        touch ${prefix}_mutect.vcf.gz.stats
        touch ${prefix}_mutect_f1r2.tar.gz
        """
}
