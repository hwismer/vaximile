include { MANTA; STRELKA; POSTPROCESS_STRELKA } from "../modules/strelka.nf"
include { BED_BGZIP_INDEX } from "../modules/utilities.nf"

workflow STRELKA_WORKFLOW {

    take:
        somatic_pairs // (somatic metamap, tumor_bam, tumor_bai, normal_bam, normal_bai)
        reference_genome // (fasta, fasta.fai, dict)
        interval_bed

    main:
        strelka_bed = BED_BGZIP_INDEX(interval_bed)
        manta = MANTA(somatic_pairs, reference_genome, strelka_bed)

        strelka_input = somatic_pairs.join(manta)
        strelka = STRELKA(strelka_input, reference_genome, strelka_bed)
        strelka_postprocess = POSTPROCESS_STRELKA(strelka.strelka_vcfs)

    emit:
        strelka_vcf = strelka_postprocess.vcf
}

