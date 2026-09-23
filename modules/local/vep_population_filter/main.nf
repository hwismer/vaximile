process VEP_POPULATION_FILTER {

    // Filter a VEP-annotated VCF with filter_vep (the gnomAD filter is set in conf/modules.config).

    label 'process_low'
    cache "lenient"


    conda "bioconda::ensembl-vep=115"
    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(somatic_meta), path(vep_vcf)
        path vep_cache
        path vep_plugins

    output:
        tuple val(somatic_meta), path("*_vep_filter.vcf"), emit: vcf

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        filter_vep -i $vep_vcf \
            -o "${prefix}_vep_filter.vcf" \
            --format vcf \
            $args
        """
    stub:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        touch "${prefix}_vep_filter.vcf"
        """

}
