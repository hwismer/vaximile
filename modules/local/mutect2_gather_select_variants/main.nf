process MUTECT2_GATHER_SELECT_VARIANTS {

    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Selecting Variants within interval for ${somatic_meta.somatic_name} in $interval_shard"

    input:
        tuple val(somatic_meta), val(somatic_name), path(vcf), path(vcf_index),  path(interval_shard)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${interval_shard}.vcf.gz"), emit: vcf
        path "versions.yml", topic: versions

    script:
    
    """
    gatk SelectVariants \
        -V $vcf \
        -L $interval_shard \
        -O "${somatic_name}_${interval_shard}.vcf.gz"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
    END_VERSIONS
    """

    stub:
    """
    touch "${somatic_name}_${interval_shard}.vcf.gz"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """

}
