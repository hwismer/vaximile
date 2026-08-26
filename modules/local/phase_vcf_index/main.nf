process PHASE_VCF_INDEX {

    /*

        Index the final phased vcf.

    */

    cpus 1
    memory "16GB"
    conda "bioconda::tabix=0.2.6"

    tag "Indexing phased vcf $phased_vcf"

    input:
        tuple val(meta), path(phased_vcf)

    output:
        tuple val(meta), path("${meta.somatic_name}_phased_annotated.vcf.gz"), path("${meta.somatic_name}_phased_annotated.vcf.gz.tbi"), emit: phased_vcf


    script:
        """
        echo ${meta.somatic_name}
        bgzip -c $phased_vcf > ${meta.somatic_name}_phased_annotated.vcf.gz

        tabix -p vcf ${meta.somatic_name}_phased_annotated.vcf.gz
        """

}
