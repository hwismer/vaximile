process PHASE_VCF_COMBINE_VARIANTS {

    /*

    Part of creating a phased germline vcf.
    Combines the variants from the tumor-only vcf and the germline_vcf (which has been renamed).


    */

    label 'process_low'
    // GATK3 ONLY - deliberately has no conda spec, and must not be given a gatk4 one.
    // This tool has no GATK4 equivalent: CombineVariants and ReadBackedPhasing were both
    // dropped in GATK4. bioconda's `gatk` 3.x is only a wrapper that needs the licensed
    // jar registered by hand, so this module stays container-only.
    container "broadinstitute/gatk3:3.6-0"

    tag "Combining somatic $tumor_only_vcf and germline $germline_vcf variants"

    input:
        tuple val(somatic_meta), path(tumor_only_vcf), path(tumor_only_index), path(germline_vcf), path(germline_index)  
        tuple path(reference_fa), path(reference_fai)
        path reference_dict


    output:
        tuple val(somatic_meta), path("*_combined_somatic_plus_germline.vcf")

    script:
        def prefix = task.ext.prefix ?: "${somatic_meta.tumor_meta.sample_name}"
        """
        java -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
                -R $reference_fa \
                --variant $germline_vcf \
                --variant $tumor_only_vcf \
                -o ${prefix}_combined_somatic_plus_germline.vcf \
                --assumeIdenticalSamples
        """

}
