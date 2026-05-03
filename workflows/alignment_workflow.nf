include { BWA_MAP;  CREATE_BWA_INDEX; } from "../modules/alignment.nf"
include { MARK_DUPLICATES_SPARK; BASE_RECALIBRATOR_SCATTER; BASE_RECALIBRATOR_GATHER; APPLY_BQSR_SCATTER; APPLY_BQSR_GATHER; GET_PILEUP_SUMMARIES } from "../modules/alignment.nf"

/*
process MARK_DUPLICATES_SPARK {

        cpus 24
            memory "48GB"
                container "broadinstitute/gatk:4.6.1.0"

                    input:
                            tuple val(meta), path(aligned_sam)

                                output:
                                        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_markdup.bam")
                                            
                                                script:
                                                    """
                                                        gatk MarkDuplicatesSpark \
                                                                -I $aligned_sam \
                                                                        -O "${meta.sample_name}_${meta.molecule}_markdup.bam" \
                                                                                --tmp-dir "\${PWD}"
                                                                                    """

}

process BASE_RECALIBRATOR {

        cpus 8
            memory "24GB"
                container "broadinstitute/gatk:4.6.1.0"

                    input:
                            tuple val(meta), path(dedup_bam)
                                    tuple path(reference_fa), path(reference_index), path(reference_dict)
                                            tuple path(known_sites_dbsnp), path(known_sites_dbsnp_index)
                                                    tuple path(known_sites_1000g_snps), path(known_sites_1000g_snps_index)
                                                            tuple path(known_indels), path(known_indels_index)
                                                                    tuple path(mills), path(mills_index)

                                                                        output:
                                                                                tuple val(meta), path("${meta.sample_name}_${meta.molecule}_recal_table.table")
                                                                                    
                                                                                        script:
                                                                                            """
                                                                                                gatk BaseRecalibrator \
                                                                                                        -I "${meta.sample_name}_${meta.molecule}_dedup.bam" \
                                                                                                                -O "${meta.sample_name}_${meta.molecule}_recal_table.table" \
                                                                                                                        -R $reference_fa \
                                                                                                                                --known-sites $known_sites_dbsnp \
                                                                                                                                        --known-sites $known_sites_1000g_snps \
                                                                                                                                                --known-sites $known_indels \
                                                                                                                                                        --known-sites $mills

                                                                                                                                                            """

}

process APPLY_BQSR {
        
            cpus 8
                memory "24GB"
                    container "broadinstitute/gatk:4.6.1.0"

                        input:
                                tuple val(meta), path(dedup_bam), path(recal_table)
                                        tuple path(reference_fa), path(reference_index), path(reference_dict)

                                            output:
                                                    tuple val(meta), path("${meta.sample_name}_${meta.molecule}_bqsr.bam"), path("${meta.sample_name}_${meta.molecule}_bqsr.bam.bai")
                                                        
                                                            script:
                                                                """
                                                                    gatk ApplyBQSR \
                                                                            -R $reference_fa \
                                                                                    -I "${meta.sample_name}_${meta.molecule}_dedup.bam" \
                                                                                            --bqsr-recal-file $recal_table \
                                                                                                    -O "${meta.sample_name}_${meta.molecule}_bqsr.bam" \
                                                                                                            --create-output-bam-index
                                                                                                                """

}

process GET_PILEUP_SUMMARIES {
        
            cpus 8
                memory "24GB"
                    container "broadinstitute/gatk:4.6.1.0"

                        input:
                                tuple val(meta), path(bqsr_bam), path(bqsr_bai)
                                        tuple path(common_germline), path(common_germline_index)

                                            output:
                                                    tuple val(meta), path("${meta.sample_name}_${meta.molecule}_pileups.table")
                                                        
                                                            script:
                                                                """
                                                                    gatk GetPileupSummaries \
                                                                            -I $bqsr_bam \
                                                                                    -V "${common_germline}" \
                                                                                            -L "${common_germline}" \
                                                                                                    -O "${meta.sample_name}_${meta.molecule}_pileups.table"

                                                                                                        """



}

*/


workflow DNA_ALIGNMENT_WORKFLOW {

    take:
        fastqs
        reference_genome
        bwa_index
        known_sites_dbsnp
        known_sites_1000g_snps
        known_indels
        mills
        common_germline
        intervals
        num_intervals
        
    main:
        
        if ( bwa_index ) {
           bwa_index_ch = Channel.fromPath(bwa_index).collect()
        } else {
            bwa_index_ch = CREATE_BWA_INDEX(reference_genome)
        }
        
        bwa_sam = BWA_MAP(fastqs, reference_genome, bwa_index_ch)
        mark_dup = MARK_DUPLICATES_SPARK(bwa_sam)
        
        base_recal_input = mark_dup.combine(intervals)
        base_recal = BASE_RECALIBRATOR_SCATTER(
            base_recal_input,
            reference_genome,
            known_sites_dbsnp,
            known_sites_1000g_snps,
            known_indels,
            mills,
        )

        base_recal_gather_input = base_recal.groupTuple(size: num_intervals)

        base_recal_gathered = BASE_RECALIBRATOR_GATHER(base_recal_gather_input)

        bqsr_input = mark_dup.join(base_recal_gathered).combine(intervals)

        bqsr = APPLY_BQSR_SCATTER(bqsr_input, reference_genome)

        bqsr_scattered = bqsr.groupTuple(size: num_intervals)

        bqsr_gather = APPLY_BQSR_GATHER(bqsr_scattered)

    emit:
        preproc_bams = bqsr_gather
        markdup_bams = mark_dup
        base_recal = base_recal_gathered

}

