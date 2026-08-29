include { MANTA } from "../../../modules/local/manta/main"
include { STRELKA } from "../../../modules/local/strelka/main"
include { POSTPROCESS_STRELKA } from "../../../modules/local/postprocess_strelka/main"

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
            tuple(meta, meta.tumor_meta.sequencing_type, tumor_bam, tumor_bai, normal_bam, normal_bai, bed, bed_tbi)
        }
        
        manta = MANTA(somatic_pairs_kit, reference_genome).dir

        strelka_input = somatic_pairs.join(manta).map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir ->
            tuple(meta.capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir)
        }
        .combine(strelka_bed, by:0)
        .map { capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir, bed, bed_tbi ->
            tuple(meta, meta.tumor_meta.sequencing_type, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir, bed, bed_tbi)
        }

        strelka = STRELKA(strelka_input, reference_genome)

        postprocess_strelka_input = strelka.strelka_vcfs
            .map { meta, strelka_snvs, strelka_snvs_index, strelka_indels, strelka_indels_index ->
                tuple(meta, meta.tumor_meta.sample_name, meta.normal_meta.sample_name, strelka_snvs, strelka_snvs_index, strelka_indels, strelka_indels_index)
            }

        strelka_postprocess = POSTPROCESS_STRELKA(postprocess_strelka_input)

    emit:
        strelka_vcf = strelka_postprocess.vcf
}
