process ADD_VCF_GT_FIELD {

    /*

    Adds the GT fields to a vcf lacking it. This is currently used for Strelka/Manta vcfs and the field is populated by 0/1 by default.

    */

    label 'process_low'
    container "griffithlab/vatools:5.2.0"

    tag "Adding 0/1 default GT to variants $somatic_vcf"

    input:
        tuple val(somatic_meta), val(tumor_sample_name), path(somatic_vcf), path(somatic_vcf_index)
         
    output:
        tuple val(somatic_meta), path("*_gt.vcf"), path("*_gt.vcf.tbi"),emit: vcf

    script:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        cp $somatic_vcf_index ${prefix}_gt.vcf.tbi
        vcf-genotype-annotator $somatic_vcf \
            "${tumor_sample_name}" \
            0/1 \
            -o "${prefix}_gt.vcf"
        """
}
