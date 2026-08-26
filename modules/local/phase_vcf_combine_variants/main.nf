process PHASE_VCF_COMBINE_VARIANTS {

    /*

    Part of creating a phased germline vcf.
    Combines the variants from the tumor-only vcf and the germline_vcf (which has been renamed).


    */

    cpus 2
    memory "16GB"
    container "broadinstitute/gatk3:3.6-0"

    tag "Combining somatic $tumor_only_vcf and germline $germline_vcf variants"

    input:
        tuple val(somatic_meta), path(tumor_only_vcf), path(tumor_only_index), path(germline_vcf), path(germline_index)  
        tuple path(reference_fa), path(reference_fai)
        path reference_dict


    output:
        tuple val(somatic_meta), path("${somatic_meta.tumor_meta.sample_name}_combined_somatic_plus_germline.vcf")

    script:
        """
        java -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
                -R $reference_fa \
                --variant $germline_vcf \
                --variant $tumor_only_vcf \
                -o ${somatic_meta.tumor_meta.sample_name}_combined_somatic_plus_germline.vcf \
                --assumeIdenticalSamples
        """

}
