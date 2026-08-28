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


    script:
        def prefix = task.ext.prefix ?: "${somatic_name}"
        """
        echo ${somatic_name}
        bgzip -c $phased_vcf > ${prefix}_phased_annotated.vcf.gz

        tabix -p vcf ${prefix}_phased_annotated.vcf.gz
        """

}
