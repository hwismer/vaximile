process VEP_ANNOTATE {

    // Annotate a VCF with VEP, using a local cache and the pVACtools plugins.

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
        """

    stub:
        def prefix = task.ext.prefix ?: "${sample_meta.sample_name}"
        """
        touch ${prefix}_vep.vcf
        touch ${prefix}_vep.vcf_summary.html
        """
}
