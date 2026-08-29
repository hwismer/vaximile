process FILTER_VCF {

    /*

    Filter VCF, retaining only PASS variants.

    */
    
    label 'process_medium'
    conda "bioconda::bcftools=1.23"
    container "staphb/bcftools:1.23"

    tag "Filtering non-passing variants from $caller $somatic_vcf"
    
    
    input:
        tuple val(somatic_meta), val(caller), path(somatic_vcf), path(tbi)
        
    output:
        tuple val(somatic_meta), val(caller), path("*_filtered_variants.vcf.gz"), path("*_filtered_variants.vcf.gz.tbi"), emit: filtered_vcf
        path "versions.yml", topic: versions


    script:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}_${caller}"
        """
        bcftools view -f PASS -Oz -o "${prefix}_filtered_variants.vcf.gz" $somatic_vcf
        bcftools index -t "${prefix}_filtered_variants.vcf.gz"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}_${caller}"
        """
        touch ${prefix}_filtered_variants.vcf.gz
        touch ${prefix}_filtered_variants.vcf.gz.tbi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: 1.23
        END_VERSIONS
        """



}
