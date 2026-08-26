process STRELKA_GERMLINE {

    /*

    Run Strelka in germline mode

    */

    cpus 8
    memory "24GB"
    container "mgibio/strelka:2.9.9"

    tag "Running Strelka in germline mode on ${meta.sample_name}"

   input:
        tuple val(meta), path(bam), path(bai), path(bed), path(bed_index)
        tuple path(reference_fa), path(reference_fai)

    output:
        tuple val(meta), val("strelka"), path("./strelka/results/variants/variants.vcf.gz"), path("./strelka/results/variants/variants.vcf.gz.tbi")

    script:

    def exome_flag = (meta.sequencing_type == "exome" || meta.sequencing_type == "exome_ffpe") ? "--exome" : ""

    """
    /opt/strelka/bin/configureStrelkaGermlineWorkflow.py \
        --bam $bam \
        ${exome_flag} \
        --callRegions $bed \
        --referenceFasta $reference_fa \
        --runDir "./strelka"

    ./strelka/runWorkflow.py -m local -j $task.cpus
    """

}
