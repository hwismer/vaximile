process FILTER_VCF {

    // Keep only PASS variants.
    
    label 'process_medium'
    conda "bioconda::bcftools=1.23.1"
    container "staphb/bcftools:1.23.1"

    tag "Filtering non-passing variants from $caller $somatic_vcf"
    
    
    input:
        // No index input: `bcftools view -f PASS` reads the VCF directly, and one
        // caller branch (Strelka via ADD_VCF_GT_FIELD) has no valid index to give.
        tuple val(somatic_meta), val(caller), path(somatic_vcf)
        
    output:
        tuple val(somatic_meta), val(caller), path("*_filtered_variants.vcf.gz"), path("*_filtered_variants.vcf.gz.tbi"), emit: filtered_vcf


    script:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}_${caller}"
        """
        bcftools view --threads $task.cpus -f PASS -Oz -o "${prefix}_filtered_variants.vcf.gz" $somatic_vcf
        bcftools index -t "${prefix}_filtered_variants.vcf.gz"
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}_${caller}"
        """
        touch ${prefix}_filtered_variants.vcf.gz
        touch ${prefix}_filtered_variants.vcf.gz.tbi
        """



}
