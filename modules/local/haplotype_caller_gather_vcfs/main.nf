process HAPLOTYPE_CALLER_GATHER_VCFS {

    cpus 2
    memory "8GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Gathering haplotype caller VCFs for ${sample_meta.sample_name}"

    input:
        tuple val(sample_meta), path(vcfs)

    output:
        tuple val(sample_meta), path("${sample_meta.sample_name}_merged.vcf")

    script:
    
    def sorted_vcfs = vcfs.sort { a, b -> a.name <=> b.name }
    def vcf_as_input = sorted_vcfs.collect { vcf ->
            "--INPUT ${vcf}"
        }.join(' ')

    """
    gatk GatherVcfs \
        $vcf_as_input \
        -O "${sample_meta.sample_name}_merged.vcf"
    
    """

}
