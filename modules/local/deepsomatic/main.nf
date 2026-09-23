process DEEPSOMATIC {

	// Run DeepSomatic on a Tumor-Normal pair.

    label 'process_max'
    container "google/deepsomatic:1.10.0"
    
    tag "DeepSomatic on ${somatic_meta.somatic_name}"
    input:
        tuple val(somatic_meta), val(tumor_sample_name), val(normal_sample_name), val(tumor_sequencing_type), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bed_regions)
        tuple path(reference_fa), path(reference_index)


    output:
        tuple val(somatic_meta), path("*_deepsomatic.vcf.gz"), path("*_deepsomatic.vcf.gz.tbi"), emit: vcf


    script:

	def model_map = [
        exome      : 'WES',
        genome     : 'WGS',
        exome_ffpe : 'FFPE_WES',
        genome_ffpe: 'FFPE_WGS'
    ]

	def model = model_map[tumor_sequencing_type]
	def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
	def args = task.ext.args ?: ''
    """
    run_deepsomatic \
        --model_type=$model \
        --ref=$reference_fa \
        --reads_normal=$normal_bam \
        --reads_tumor=$tumor_bam \
        --output_vcf=${prefix}_deepsomatic.vcf.gz \
        --output_gvcf=${prefix}_deepsomatic.gvcf.gz \
        --sample_name_tumor=${tumor_sample_name} \
        --sample_name_normal=${normal_sample_name} \
        --num_shards=$task.cpus \
        --logging_dir=./logs \
        $args \
        --regions=$bed_regions

    """

    stub:
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
    """
    touch ${prefix}_deepsomatic.vcf.gz
    touch ${prefix}_deepsomatic.vcf.gz.tbi
    """

}
