process PHASE_VCF_INDEX {

    /*

        Index the final phased vcf.

    */

    label 'process_low'
    conda "bioconda::tabix=0.2.6"

    tag "Indexing phased vcf $phased_vcf"

    input:
        tuple val(meta), val(somatic_name), path(phased_vcf)

    output:
        tuple val(meta), path("*_phased_annotated.vcf.gz"), path("*_phased_annotated.vcf.gz.tbi"), emit: phased_vcf
        path "versions.yml", topic: versions


    script:
        def prefix = task.ext.prefix ?: "${somatic_name}"
        """
        echo ${somatic_name}
        bgzip -c $phased_vcf > ${prefix}_phased_annotated.vcf.gz

        tabix -p vcf ${prefix}_phased_annotated.vcf.gz
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tabix: 0.2.6
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_name}"
        """
        touch ${prefix}_phased_annotated.vcf.gz
        touch ${prefix}_phased_annotated.vcf.gz.tbi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tabix: 0.2.6
        END_VERSIONS
        """

}
