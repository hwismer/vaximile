process VEP_POPULATION_FILTER {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    cpus 2
    memory "16GB"
    cache "lenient"


    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(somatic_meta), path(vep_vcf)
        path vep_cache
        path vep_plugins

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_vep_filter.vcf"), emit: vcf

    script:
        """
        filter_vep -i $vep_vcf \
            -o "${somatic_meta.somatic_name}_vep_filter.vcf" \
            --format vcf \
            --filter "gnomADe_AF < 0.001 or not gnomADe_AF"
        """
}
