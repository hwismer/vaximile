process HAPLOTYPE_CALLER_GATHER_VCFS {

    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Gathering haplotype caller VCFs for ${sample_meta.sample_name}"

    input:
        tuple val(sample_meta), path(vcfs)

    output:
        tuple val(sample_meta), path("*_merged.vcf"), emit: vcf
        path "versions.yml", topic: versions

    script:

    def prefix = task.ext.prefix ?: "${sample_meta.sample_name}"
    def sorted_vcfs = vcfs.sort { a, b -> a.name <=> b.name }
    def vcf_as_input = sorted_vcfs.collect { vcf ->
            "--INPUT ${vcf}"
        }.join(' ')

    """
    gatk GatherVcfs \
        $vcf_as_input \
        -O "${prefix}_merged.vcf"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
    END_VERSIONS
    """

    stub:

    def prefix = task.ext.prefix ?: "${sample_meta.sample_name}"

    """
    touch ${prefix}_merged.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """

}
