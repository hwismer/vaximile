process PHASE_VCF_SELECT_VARIANTS {

    /*

    Part of creating a phased germline vcf
    Takes a somatic vcf and extracts just the tumor sample.

    */

    cpus 2
    memory "32GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Extracting tumor sample from $somatic_vcf"

    input:

        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index)
        tuple path(reference_fa), path(reference_fai)
        path reference_dict

    output:
        tuple val(somatic_meta), path("${somatic_meta.tumor_meta.sample_name}_tumor_only.vcf.gz"), path("${somatic_meta.tumor_meta.sample_name}_tumor_only.vcf.gz.tbi")

    script:
        """
        gatk SelectVariants \
            -V $somatic_vcf \
            -R "${reference_fa}" \
            --sample-name ${somatic_meta.tumor_meta.sample_name} \
            --create-output-variant-index \
            -O ${somatic_meta.tumor_meta.sample_name}_tumor_only.vcf.gz
        """

}
