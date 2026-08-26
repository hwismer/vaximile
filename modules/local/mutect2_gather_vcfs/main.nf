process MUTECT2_GATHER_VCFS {

    cpus 2
    memory "8GB"
    container "broadinstitute/gatk:4.6.1.0"
    
    tag "Gathering Mutect VCFs from ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(vcfs)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_merged.vcf")

    script:
    
    def sorted_vcfs = vcfs.sort { a, b -> a.name <=> b.name }
    def vcf_as_input = sorted_vcfs.collect { vcf ->
            "--INPUT ${vcf}"
        }.join(' ')

    """
    gatk GatherVcfs \
        $vcf_as_input \
        -O "${somatic_meta.somatic_name}_merged.vcf"
    
    """

}
