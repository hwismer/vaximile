process FILTER_VCF {

    /*

    Filter VCF, retaining only PASS variants.

    */

    label 'process_medium'
    conda "bioconda::bcftools=1.23"
    container "staphb/bcftools:1.23"
    
    tag "Filtering VCF $vcf"

    input:
        tuple val(meta), val(caller), path(vcf), path(tbi)

    output:
        tuple val(meta), val(caller), path("*_filtered_variants.vcf.gz"), path("*_filtered_variants.vcf.gz.tbi"), emit: filtered_vcf


    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${caller}"
        """
        bcftools view -f PASS -Oz -o "${prefix}_filtered_variants.vcf.gz" $vcf
        bcftools index -t "${prefix}_filtered_variants.vcf.gz"
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${caller}"
        """
        touch ${prefix}_filtered_variants.vcf.gz
        touch ${prefix}_filtered_variants.vcf.gz.tbi
        """



}
