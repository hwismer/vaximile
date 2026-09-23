process POSTPROCESS_VCF {

    // Normalise, split multiallelics, sort, deduplicate and index a VCF.

    label 'process_medium'
    conda "bioconda::bcftools=1.23.1"
    container "staphb/bcftools:1.23.1"

    tag "Normalizing $somatic_vcf"

    input:
        tuple val(meta), val(somatic_name), val(caller), path(somatic_vcf), path(somatic_vcf_index)
        tuple path(reference_fa), path(reference_index_dir)

    output:
        tuple val(meta), val(caller),
            path("${meta.somatic_name}_${caller}_variants.vcf.gz"),
            path("${meta.somatic_name}_${caller}_variants.vcf.gz.tbi"), emit: vt_vcf

    script:
        // ext.prefix intentionally omitted: the output stem "*_variants.vcf.gz" would also
        // match the staged input "*_filtered_variants.vcf.gz", so no safe output glob exists.
        def args = task.ext.args ?: ''
        """
        bcftools norm --threads $task.cpus $args -f $reference_fa $somatic_vcf -Oz -o norm_vcf.vcf.gz
        bcftools sort norm_vcf.vcf.gz -Oz -o "${somatic_name}_${caller}_variants.vcf.gz"
        bcftools index -t "${somatic_name}_${caller}_variants.vcf.gz"
        """

    stub:
        // The output paths are literal names built from meta.somatic_name, so the stub
        // must use the same expression rather than the somatic_name val the script uses.
        """
        touch ${meta.somatic_name}_${caller}_variants.vcf.gz
        touch ${meta.somatic_name}_${caller}_variants.vcf.gz.tbi
        """
}
