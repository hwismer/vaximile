process PHASE_VCF_SELECT_VARIANTS {

    /*

    Part of creating a phased germline vcf
    Takes a somatic vcf and extracts just the tumor sample.

    */

    label 'process_medium'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Extracting tumor sample from $somatic_vcf"

    input:

        tuple val(somatic_meta), val(tumor_sample_name), path(somatic_vcf), path(somatic_vcf_index)
        tuple path(reference_fa), path(reference_fai)
        path reference_dict

    output:
        tuple val(somatic_meta), path("*_tumor_only.vcf.gz"), path("*_tumor_only.vcf.gz.tbi")

    script:
        def prefix = task.ext.prefix ?: "${tumor_sample_name}"
        """
        gatk SelectVariants \
            -V $somatic_vcf \
            -R "${reference_fa}" \
            --sample-name ${tumor_sample_name} \
            --create-output-variant-index \
            -O ${prefix}_tumor_only.vcf.gz
        """

    stub:
        def prefix = task.ext.prefix ?: "${tumor_sample_name}"
        """
        touch ${prefix}_tumor_only.vcf.gz
        touch ${prefix}_tumor_only.vcf.gz.tbi
        """

}
