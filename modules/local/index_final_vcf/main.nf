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
        tuple val(meta), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
        """
        bgzip $vcf
        tabix -p vcf "${vcf}.gz"
        """

    stub:
        """
        touch "${vcf}.gz"
        touch "${vcf}.gz.tbi"
        """

}
