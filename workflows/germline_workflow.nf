include { HAPLOTYPE_CALLER_SCATTER; HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS; HAPLOTYPE_CALLER_GATHER_SELECT_VARIANTS; HAPLOTYPE_CALLER_GATHER_VCFS; HAPLOTYPE_CALLER_FILTER_VARIANTS; POSTPROCESS_GERMLINE; } from "../modules/germline.nf"

workflow HAPLOTYPE_CALLER {

    take:
        sample_bams // (sample metamap, bam, bai)
        intervals // (shard name, interval_shard)
        num_intervals // int
        interval_padding
        reference_genome // (fasta, fasta.fai, dict)
        hapmap
        mills

    main:
        
        bams_intervals = sample_bams.combine(intervals)
        haplotype_scatter = HAPLOTYPE_CALLER_SCATTER(bams_intervals, reference_genome, interval_padding)
        cnn_score = HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS(haplotype_scatter, reference_genome, interval_padding)
        select_variants = HAPLOTYPE_CALLER_GATHER_SELECT_VARIANTS(cnn_score)
        gathered_vcf = HAPLOTYPE_CALLER_GATHER_VCFS(select_variants.groupTuple(size: num_intervals))
        filtered_vcf = HAPLOTYPE_CALLER_FILTER_VARIANTS(gathered_vcf, hapmap, mills)
        final_vcf = POSTPROCESS_GERMLINE(filtered_vcf, reference_genome)
    
    emit:
        germline_vcf = final_vcf

}


