process PHASE_VCF_SORT_VCF {

    /*

    Part of creating a phased vcf.
    Sort the somatic + germline combined vcf.

    */

    label 'process_medium'
    container 'broadinstitute/picard:3.4.0'

    tag "Sorting VCF ${combined_vcf}"

    input:
        tuple val(somatic_meta), path(combined_vcf)
        tuple path(reference_fa), path(reference_fai)
        path(reference_dict)

    output:
        tuple val(somatic_meta), path("*_combined.sorted.vcf"), emit: sorted_vcf

    script:
        def prefix = task.ext.prefix ?: "${somatic_meta.tumor_meta.sample_name}"
        """
        java -jar /usr/picard/picard.jar \
            SortVcf \
                -I $combined_vcf \
                -O ${prefix}_combined.sorted.vcf \
                -SD $reference_dict
        """

}
