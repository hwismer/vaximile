process DEEPVARIANT {

    label 'process_very_high'
    container "google/deepvariant:1.10.0"
    
    tag "DeepVariant on ${meta.sample_name}"

    input:
        tuple val(meta), val(sequencing_type), path(bam), path(bai), path(bed)
        tuple path(reference_fa), path(reference_index)

    output:
        tuple val(meta), val("deepvariant"), path("*_deepvariant.vcf.gz"), path("*_deepvariant.vcf.gz.tbi"), emit: vcf

    script:

	def model_map = [
        exome      : 'WES',
        genome     : 'WGS',
        // DeepVariant has no FFPE models.
        exome_ffpe : 'WES',
        genome_ffpe: 'WGS'
    ]

	def model = model_map[sequencing_type]

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    # Intermediates go in the task directory rather than \$TMPDIR, and are removed once the VCF is written.
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=$model \
        --ref=$reference_fa \
        --reads=$bam \
        --output_vcf=${prefix}_deepvariant.vcf.gz \
        --num_shards=$task.cpus \
        --intermediate_results_dir=intermediate_results \
        --regions $bed \
        $args
    rm -rf intermediate_results
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    touch ${prefix}_deepvariant.vcf.gz
    touch ${prefix}_deepvariant.vcf.gz.tbi
    """
}
