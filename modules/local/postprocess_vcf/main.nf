process POSTPROCESS_VCF {

    /*

        Postprocess somatic variants from a vcf. This is used to normalize variants and decompose biallelics
        into different entries. The vcf is then sorted and duplicate entries removed. The final vcf is then indexed.

    */

    label 'process_medium'
    container "staphb/bcftools:1.23"

    tag "Normalizing $somatic_vcf"

    input:
        tuple val(meta), val(caller), path(somatic_vcf), path(somatic_vcf_index)
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
        bcftools norm $args -f $reference_fa $somatic_vcf -Oz -o norm_vcf.vcf.gz
        bcftools sort norm_vcf.vcf.gz -Oz -o "${meta.somatic_name}_${caller}_variants.vcf.gz"
        bcftools index -t "${meta.somatic_name}_${caller}_variants.vcf.gz"
        """
}
