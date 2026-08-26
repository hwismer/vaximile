include { DEEPSOMATIC } from "../../../modules/local/deepsomatic/main"

workflow DEEPSOMATIC_WORKFLOW {

    take:
        somatic_pairs // (somatic metamap, tumor_bam, tumor_bai, normal_bam, normal_bai)
        reference_genome // (fasta, fasta.fai, dict)
        capture_kits // (kit_name, bed file)

    main:

        somatic_pairs_kit = somatic_pairs.map{ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            tuple(meta.capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai)
        }
        .combine(capture_kits, by:0)
        .map { kit, meta, tb, tbai, nbam, nbai, bed -> tuple(meta, tb, tbai, nbam, nbai, bed)
        }

        deepsomatic = DEEPSOMATIC(somatic_pairs_kit, reference_genome)

    emit:
        deepsomatic_vcf = deepsomatic
}
