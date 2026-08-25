include { MANTA; STRELKA; POSTPROCESS_STRELKA } from "../modules/local/strelka_somatic.nf"
include { CAPTURE_KIT_BED_PROCESS } from "../modules/local/utilities.nf"

workflow STRELKA_WORKFLOW {

    take:
        somatic_pairs // (somatic metamap, tumor_bam, tumor_bai, normal_bam, normal_bai)
        reference_genome // (fasta, fasta.fai, dict)
        strelka_bed

    main:
        
        somatic_pairs_kit = somatic_pairs.map{ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            tuple(meta.capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai)
        }
        .combine(strelka_bed, by:0)
        .map{ capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bed, bed_tbi ->
            tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bed, bed_tbi)
        }
        
        manta = MANTA(somatic_pairs_kit, reference_genome)

        strelka_input = somatic_pairs.join(manta).map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir ->
            tuple(meta.capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir)
        }
        .combine(strelka_bed, by:0)
        .map { capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir, bed, bed_tbi ->
            tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir, bed, bed_tbi)
        }

        strelka = STRELKA(strelka_input, reference_genome)
        strelka_postprocess = POSTPROCESS_STRELKA(strelka.strelka_vcfs)

    emit:
        strelka_vcf = strelka_postprocess.vcf
}

