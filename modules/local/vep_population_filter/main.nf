process VEP_POPULATION_FILTER {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

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
        path "versions.yml", topic: versions

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        filter_vep -i $vep_vcf \
            -o "${prefix}_vep_filter.vcf" \
            --format vcf \
            $args
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ensembl-vep: \$(vep --help 2>&1 | grep -Eo 'ensembl-vep +: *[0-9.]+' | tr -d ' ' | cut -d: -f2)
        END_VERSIONS
        """
    stub:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        touch "${prefix}_vep_filter.vcf"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ensembl-vep: 115
        END_VERSIONS
        """

}
