process DEEPVARIANT {

    cpus 16
    memory "32GB"
    container "google/deepvariant:1.10.0"
    
    tag "DeepVariant on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai), path(bed)
        tuple path(reference_fa), path(reference_index)

    output:
        tuple val(meta), val("deepvariant"), path("${meta.sample_name}_deepvariant.vcf.gz"), path("${meta.sample_name}_deepvariant.vcf.gz.tbi"), emit: vcf
        tuple val(meta), val("deepvariant"), path("${meta.sample_name}_deepvariant.gvcf.gz"), path("${meta.sample_name}_deepvariant.gvcf.gz.tbi"), emit: gvcf

    script:

	def model_map = [
        exome      : 'WES',
        genome     : 'WGS',
    ]

	def model = model_map[meta.sequencing_type]

    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=$model \
        --ref=$reference_fa \
        --reads=$bam \
        --output_vcf=${meta.sample_name}_deepvariant.vcf.gz \
        --output_gvcf=${meta.sample_name}_deepvariant.gvcf.gz \
        --num_shards=$task.cpus \
        --vcf_stats_report=true \
        --regions $bed \
        --disable_small_model=true
    """
}
