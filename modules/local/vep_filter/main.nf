process VEP_FILTER {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    label 'process_low'
    conda "bioconda::ensembl-vep=115"
    container "ensemblorg/ensembl-vep:release_115.0"

    tag "Filtering $vep_vcf with VEP based on population frequencies"

    input:
        tuple val(meta), path(vep_vcf)
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("${meta.somatic_name}_vep.vcf"), emit: vep_vcf

    script:
        def args = task.ext.args ?: ''
        """

        filter_vep -i $vep_vcf \
            -o "${meta.somatic_name}_vep.vcf" \
            --format vcf \
            $args

        """
}
