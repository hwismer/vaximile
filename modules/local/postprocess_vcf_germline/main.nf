process POSTPROCESS_VCF {

    /*

        Postprocess somatic variants from a vcf. This is used to normalize variants and decompose biallelics
        into different entries. The vcf is then sorted and duplicate entries removed. The final vcf is then indexed.

    */

    label 'process_medium'
    conda "bioconda::bcftools=1.23"
    container "staphb/bcftools:1.23"

    tag "Normalizing $vcf"

    input:
        tuple val(meta), val(sample_name), val(caller), path(vcf), path(vcf_index)
        tuple path(reference_fa), path(reference_index)

    output:
        tuple val(meta), val(caller),
            path("${meta.sample_name}_${caller}_variants.vcf.gz"),
            path("${meta.sample_name}_${caller}_variants.vcf.gz.tbi"), emit: postproc_vcf

    script:

        """
        bcftools norm -m -any -d exact -f $reference_fa $vcf -Oz -o norm_vcf.vcf.gz
        bcftools sort norm_vcf.vcf.gz -Oz -o "${sample_name}_${caller}_variants.vcf.gz"
        bcftools index -t "${sample_name}_${caller}_variants.vcf.gz"
        """
}
