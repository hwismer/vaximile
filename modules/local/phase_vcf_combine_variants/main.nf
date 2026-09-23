process PHASE_VCF_COMBINE_VARIANTS {

    // Combine the tumour-only and renamed germline VCFs for phasing.

    label 'process_low'
    // GATK3 only: CombineVariants/ReadBackedPhasing have no GATK4 equivalent, and bioconda's gatk3 needs a licensed jar.
    container "broadinstitute/gatk3:3.6-0"

    tag "Combining somatic $tumor_only_vcf and germline $germline_vcf variants"

    input:
        tuple val(somatic_meta), path(tumor_only_vcf), path(tumor_only_index), path(germline_vcf), path(germline_index)  
        tuple path(reference_fa), path(reference_fai)
        path reference_dict


    output:
        tuple val(somatic_meta), path("*_combined_somatic_plus_germline.vcf"), emit: vcf

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

    stub:
        def prefix = task.ext.prefix ?: "${somatic_meta.tumor_meta.sample_name}"
        """
        touch ${prefix}_combined_somatic_plus_germline.vcf
        """

}
