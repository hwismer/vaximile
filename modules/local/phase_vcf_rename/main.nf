process PHASE_VCF_RENAME {

    // Rename the germline VCF's sample to the tumour's, so the two can be combined.

    label 'process_low'
    conda "bioconda::bcftools=1.23.1"
    container "staphb/bcftools:1.23.1"

    tag "Renaming germline sample ${normal_meta.sample_name} to ${somatic_meta.tumor_meta.sample_name}"

    input:
        tuple val(normal_meta), val(somatic_meta), val(normal_sample_name), val(tumor_sample_name), path(germline_vcf), path(germline_vcf_index)

    output:
        tuple val(somatic_meta), path("*_germline_rename.vcf.gz"), path("*_germline_rename.vcf.gz.tbi"), emit: vcf

    script:
        def prefix = task.ext.prefix ?: "${tumor_sample_name}"
        """
        cat > sample_map.txt <<EOF
        ${normal_sample_name} ${tumor_sample_name}
        EOF
        
        bcftools reheader \
            -N sample_map.txt \
            --threads $task.cpus \
            -o ${prefix}_germline_rename.vcf.gz \
            $germline_vcf

        bcftools index -t --threads $task.cpus ${prefix}_germline_rename.vcf.gz
        """

    stub:
        def prefix = task.ext.prefix ?: "${tumor_sample_name}"
        """
        touch ${prefix}_germline_rename.vcf.gz
        touch ${prefix}_germline_rename.vcf.gz.tbi
        """


}
