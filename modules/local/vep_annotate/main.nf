process VEP_ANNOTATE {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    label 'process_medium'

    conda "bioconda::ensembl-vep=115"
    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(somatic_meta), path(vcf)
        tuple path(reference_fa), path(reference_index)
        path vep_cache
        path vep_plugins

    output:
        tuple val(somatic_meta), path("*_vep.vcf"), emit: vcf
        tuple val(somatic_meta), path("*.html"), emit: report

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        vep \
            --input_file $vcf  \
            --output_file ${prefix}_vep.vcf \
            --format vcf --vcf \
            $args \
            --fasta $reference_fa  \
            --offline --cache \
            --plugin Frameshift --plugin Wildtype \
            --dir_plugins $vep_plugins \
            --dir_cache $vep_cache
        """
}
