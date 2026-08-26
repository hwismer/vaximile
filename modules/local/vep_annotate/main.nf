process VEP_ANNOTATE {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    cpus 4
    memory "16GB"
    
    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(somatic_meta), path(vcf)
        tuple path(reference_fa), path(reference_index)
        path vep_cache
        path vep_plugins

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_vep.vcf"), emit: vcf
        tuple val(somatic_meta), path("*.html"), emit: report

    script:
        """
        vep \
            --input_file $vcf  \
            --output_file ${somatic_meta.somatic_name}_vep.vcf \
            --everything \
            --format vcf --vcf --symbol --terms SO --mane_select --canonical --tsl --biotype --hgvs \
            --fasta $reference_fa  \
            --offline --cache \
            --plugin Frameshift --plugin Wildtype \
            --pick \
            --dir_plugins $vep_plugins \
            --dir_cache $vep_cache \
            --transcript_version
        """
}
