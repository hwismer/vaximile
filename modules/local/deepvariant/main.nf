process DEEPVARIANT {

    label 'process_very_high'
    container "google/deepvariant:1.10.0"
    
    tag "DeepVariant on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai), path(bed)
        tuple path(reference_fa), path(reference_index)

    output:
        tuple val(meta), val("deepvariant"), path("*_deepvariant.vcf.gz"), path("*_deepvariant.vcf.gz.tbi"), emit: vcf
        tuple val(meta), val("deepvariant"), path("*_deepvariant.gvcf.gz"), path("*_deepvariant.gvcf.gz.tbi"), emit: gvcf

    script:

	def model_map = [
        exome      : 'WES',
        genome     : 'WGS',
    ]

	def model = model_map[meta.sequencing_type]

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=$model \
        --ref=$reference_fa \
        --reads=$bam \
        --output_vcf=${prefix}_deepvariant.vcf.gz \
        --output_gvcf=${prefix}_deepvariant.gvcf.gz \
        --num_shards=$task.cpus \
        --regions $bed \
        $args
    """
}
