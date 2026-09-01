process VEP_ANNOTATE {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    label 'process_high'
    conda "bioconda::ensembl-vep=115"
    container "ensemblorg/ensembl-vep:release_115.0"

    tag "VEP on ${sample_meta.sample_name}"
    input:
        tuple val(sample_meta), path(vcf)
        tuple path(reference_fa), path(reference_index)
        path vep_cache
        path vep_plugins

    output:
        tuple val(sample_meta), path("*_vep.vcf"), emit: vcf
        tuple val(sample_meta), path("*.html"), emit: report
        path "versions.yml", topic: versions

    script:
        def prefix = task.ext.prefix ?: "${sample_meta.sample_name}"
        def args = task.ext.args ?: ''
        """
        vep \
            --fork $task.cpus \
            --input_file $vcf  \
            --output_file ${prefix}_vep.vcf \
            --format vcf --vcf \
            --fasta $reference_fa  \
            --offline --cache \
            --plugin Frameshift --plugin Wildtype \
            --dir_plugins $vep_plugins \
            --dir_cache $vep_cache \
            $args
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ensembl-vep: \$(vep --help 2>&1 | grep -Eo 'ensembl-vep +: *[0-9.]+' | tr -d ' ' | cut -d: -f2)
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${sample_meta.sample_name}"
        """
        touch ${prefix}_vep.vcf
        touch ${prefix}_vep.vcf_summary.html
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ensembl-vep: 115
        END_VERSIONS
        """
}
