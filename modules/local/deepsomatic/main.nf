process DEEPSOMATIC {

    cpus 32
    memory "64GB"
    container "google/deepsomatic:1.10.0"
    
    tag "DeepSomatic on ${somatic_meta.somatic_name}"
    input:
        tuple val(somatic_meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bed_regions)
        tuple path(reference_fa), path(reference_index)


    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_deepsomatic.vcf.gz"), path("${somatic_meta.somatic_name}_deepsomatic.vcf.gz.tbi")


    script:

	def model_map = [
        exome      : 'WES',
        genome     : 'WGS',
        exome_ffpe : 'FFPE_WES',
        genome_ffpe: 'FFPE_WGS'
    ]

	def model = model_map[somatic_meta.tumor_meta.sequencing_type]
    """
    run_deepsomatic \
        --model_type=$model \
        --ref=$reference_fa \
        --reads_normal=$normal_bam \
        --reads_tumor=$tumor_bam \
        --output_vcf=${somatic_meta.somatic_name}_deepsomatic.vcf.gz \
        --output_gvcf=${somatic_meta.somatic_name}_deepsomatic.gvcf.gz \
        --sample_name_tumor=${somatic_meta.tumor_meta.sample_name} \
        --sample_name_normal=${somatic_meta.normal_meta.sample_name}\
        --num_shards=$task.cpus \
        --logging_dir=./logs \
        --vcf_stats_report=true \
        --use_default_pon_filtering=true \
        --regions=$bed_regions

    """

}
