include { STAR_FUSION } from "../../../modules/local/star_fusion/main"
include { ARRIBA_FUSION } from "../../../modules/local/arriba_fusion/main"

workflow FUSION_CALLING {

    take:
        star_bam
        star_chimeric_out
        reference_genome
        gtf
        ctat_resource_lib
        arriba_resources

    main:
        star_fusion = STAR_FUSION(star_chimeric_out, ctat_resource_lib)
        arriba = ARRIBA_FUSION(star_bam, reference_genome, gtf, arriba_resources)
    emit:
        arriba_fusion = arriba.arriba_fusions
        star_fusion = star_fusion.fusion_preds

}
