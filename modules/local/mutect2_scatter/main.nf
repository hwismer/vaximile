process MUTECT2_SCATTER {


    /*

        Call mutect 2 on a single interval shard. This process creates a new metadata block specific to
        the somatic caller of the form [somatic_name, somatic_caller, tumor_metadata, normal_metadata].

    */

    cpus 4
    memory "12GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Running Mutect2 scatter on ${somatic_meta.somatic_name} at $interval_shard"

    input:
        tuple val(somatic_meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict
        tuple path(germline_resource), path(germline_resource_index)
        tuple path(pon), path(pon_index)
        val interval_padding

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${interval_shard}_mutect.vcf.gz"), 
        path("${somatic_meta.somatic_name}_${interval_shard}_mutect.vcf.gz.tbi"), path(interval_shard), emit: vcf
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${interval_shard}_mutect_f1r2.tar.gz"), emit: f1r2
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${interval_shard}_mutect.vcf.gz.stats"), emit: stats

    script:

        """
        gatk Mutect2 \
            -R "${reference_fa}" \
            -I ${tumor_bam} \
            -I ${normal_bam} \
            -normal ${somatic_meta.normal_meta.sample_name} \
            --germline-resource $germline_resource \
            --panel-of-normals $pon \
            --f1r2-tar-gz "${somatic_meta.somatic_name}_${interval_shard}_mutect_f1r2.tar.gz" \
            -L $interval_shard \
            -ip $interval_padding \
            -O "${somatic_meta.somatic_name}_${interval_shard}_mutect.vcf.gz" \
            --native-pair-hmm-threads $task.cpus
        """
}
