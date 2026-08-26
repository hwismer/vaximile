process ADD_VCF_GT_FIELD {

    /*

    Adds the GT fields to a vcf lacking it. This is currently used for Strelka/Manta vcfs and the field is populated by 0/1 by default.

    */

    cpus 2
    memory "16GB"
    container "griffithlab/vatools:5.2.0"

    tag "Adding 0/1 default GT to variants $somatic_vcf"

    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index)
         
    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_gt.vcf"), path("${somatic_meta.somatic_name}_gt.vcf.tbi"),emit: vcf

    script:
        """
        cp $somatic_vcf_index ${somatic_meta.somatic_name}_gt.vcf.tbi
        vcf-genotype-annotator $somatic_vcf \
            "${somatic_meta.tumor_meta.sample_name}" \
            0/1 \
            -o "${somatic_meta.somatic_name}_gt.vcf"
        """
}
