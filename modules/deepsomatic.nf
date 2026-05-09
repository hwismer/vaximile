process DEEPSOMATIC {
    
    container "google/deepsomatic:1.10.0"
    
    //publishDir "${params.outdir}/deepsomatic/", mode: "copy"

    input:
        tuple val(somatic_meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        tuple path(reference_fa), path(reference_fai), path(reference_dict)
        path(bed)

    output:
        tuple val(somatic_meta), path("${somatic_name}_deepsomatic.vcf.gz"), path("${somatic_name}_deepsomatic.vcf.gz.tbi"), emit:deepsomatic_vcf
        tuple val(somatic_meta), path(

    script:
    def seq_type_map = [
        "exome"       : "WES",
        "genome"      : "WGS",
        "genome_FFPE" : "WGS_FFPE",
        "exome_FFPE"  : "WES_FFPE"
    ]
    def model = seq_type_map[tumor_meta.sequencing_type] ?: "Unknown"
    """
    run_deepsomatic \
        --model_type=${model} \
        --ref="${reference_fa}" \
        --reads_normal="${normal_bam}" \
        --reads_tumor="${tumor_bam}" \
        --output_vcf="${somatic_meta.somatic_name}_deepsomatic.vcf.gz" \
        --output_gvcf="${somatic_meta.somatic_name}_deepsomatic.g.vcf.gz" \
        --sample_name_tumor="${somatic_meta.tumor_meta.sample_name}" \
        --sample_name_normal="${somatic_meta.normal_meta.sample_name}" \
        --use_default_pon_filtering=true \
        --num_shards=$task.cpus \
        --logging_dir=./logs \
        --intermediate_results_dir ./intermediate_dir

    """
}
