process ANNOTATE_VCF_COVERAGE {

    /*

        Use the output of bamreadcount to annotate coverage given a particular sample.
    */

    label 'process_low_memory'

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_meta), val(sample_meta), val(sample_name), val(molecule), val(somatic_name), path(indels), path(snvs), path(vcf)
    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${sample_meta.sample_name}_${sample_meta.molecule}_coverage.vcf"), emit: vcf
        
    script:
        """
        vcf-readcount-annotator \
            $vcf \
            $indels \
            ${molecule} \
            -s ${sample_name} \
            -t indel \
            -o vcf1.vcf

        vcf-readcount-annotator \
            vcf1.vcf \
            $snvs \
            ${molecule} \
            -s ${sample_name} \
            -t snv \
            -o ${somatic_name}_${sample_name}_${molecule}_coverage.vcf
        """

    stub:
        """
        touch ${somatic_meta.somatic_name}_${sample_meta.sample_name}_${sample_meta.molecule}_coverage.vcf
        """

}
