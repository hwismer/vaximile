process ANNOTATE_VCF_COVERAGE {

    /*

        Use the output of bamreadcount to annotate coverage given a particular sample.
    */

    cpus 4
    memory "16GB"

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_meta), val(sample_meta), path(indels), path(snvs), path(vcf)
    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${sample_meta.sample_name}_${sample_meta.molecule}_coverage.vcf")
        
    script:
        """
        vcf-readcount-annotator \
            $vcf \
            $indels \
            ${sample_meta.molecule} \
            -s ${sample_meta.sample_name} \
            -t indel \
            -o vcf1.vcf

        vcf-readcount-annotator \
            vcf1.vcf \
            $snvs \
            ${sample_meta.molecule} \
            -s ${sample_meta.sample_name} \
            -t snv \
            -o ${somatic_meta.somatic_name}_${sample_meta.sample_name}_${sample_meta.molecule}_coverage.vcf
        """

}
