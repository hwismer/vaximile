process INDEX_FINAL_VCF {
    /*

    tabix index a vcf file

    */

    label 'process_low'
    conda "bioconda::tabix=0.2.6"

    tag "Indexing final vcf $vcf"

    input:
        tuple val(meta), path(vcf)
        

    output:
        tuple val(meta), path("${vcf}.gz"), path("${vcf}.gz.tbi"), emit: out
        path "versions.yml", topic: versions

    script:
        """
        bgzip $vcf
        tabix -p vcf "${vcf}.gz"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tabix: 0.2.6
        END_VERSIONS
        """

    stub:
        """
        touch "${vcf}.gz"
        touch "${vcf}.gz.tbi"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tabix: 0.2.6
        END_VERSIONS
        """

}
