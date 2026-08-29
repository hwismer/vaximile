process INDEX_VCF {
    
    // TBI index a vcf file


    label 'process_low'
    conda "bioconda::bcftools=1.23"
    container "staphb/bcftools:1.23"

    tag "Indexing $vcf"

    input:
        tuple val(vcf_name), val(meta), path(vcf) // vcf name should be a string that corresponds to the file name ie. Sample1
        val(filename_suffix) // suffix corresponds to name after vcf_name ie somatic_variants -> Sample1_somatic_variants.vcf.gz

    output:
        tuple val(meta), path("${vcf_name}_${filename_suffix}.vcf.gz"), path("${vcf_name}_${filename_suffix}.vcf.gz.tbi"), emit: vcf
        path "versions.yml", topic: versions

    script:
        """
        bcftools view $vcf -Oz -o "${vcf_name}_${filename_suffix}.vcf.gz"
        bcftools index -t ${vcf_name}_${filename_suffix}.vcf.gz
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
        END_VERSIONS
        """

    stub:
        """
        touch ${vcf_name}_${filename_suffix}.vcf.gz
        touch ${vcf_name}_${filename_suffix}.vcf.gz.tbi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: 1.23
        END_VERSIONS
        """

}
