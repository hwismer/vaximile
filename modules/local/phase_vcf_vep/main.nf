process PHASE_VCF_VEP {

    /*

        Use VEP to annotated the phased vcf.

    */
    label 'process_high'
    conda "bioconda::ensembl-vep=115"
    container "ensemblorg/ensembl-vep:release_115.0"

    tag "VEP on phased vcf $phased_vcf"

    input:
        tuple val(meta), path(phased_vcf)
        path reference_fa
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("phased_vcf_vep.vcf"), emit: vcf
        path "versions.yml", topic: versions

    script:
        def args = task.ext.args ?: ''
        """
        vep \
            --input_file $phased_vcf  \
            --output_file "phased_vcf_vep.vcf" \
            --format vcf --vcf \
            $args \
            --fasta $reference_fa  \
            --offline --cache $vep_cache \
            --plugin Frameshift --plugin Wildtype --plugin Downstream \
            --fork ${task.cpus} \
            --dir_plugins $vep_plugins

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ensembl-vep: \$(vep --help 2>&1 | grep -Eo 'ensembl-vep +: *[0-9.]+' | tr -d ' ' | cut -d: -f2)
        END_VERSIONS
        """

    stub:
        """
        touch phased_vcf_vep.vcf
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ensembl-vep: 115
        END_VERSIONS
        """

}
