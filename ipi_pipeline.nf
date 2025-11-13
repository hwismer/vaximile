// Default parameter input

params.normal_fastq_dir
params.tumor_fastq_dir
params.outdir = "./"

/*
process DEEPSOMATIC {

    container "/c4/home/hwismer/pvac_pipeline/containers/deepsomatic_1.9.0.sif"
    
    publishDir "${params.outdir}/deepsomatic/", mode: "copy"

    input:
        tuple val(tumor_sample_id), val(tumor_type),  path(tumor_reads), path(tumor_reads_index)
        tuple val(normal_sample_id), val(normal_type), path(normal_reads), path(normal_reads_index)
        val(somatic_name)
        path reference_fa
        path reference_index_dir

    output:
        tuple val(somatic_name), val("DeepSomatic"), path("${somatic_name}_deepsomatic_filtered.vcf.gz")

    script:
        """
        run_deepsomatic \
            --model_type=WES \
            --ref="${reference_index_dir}/${reference_fa}" \
            --reads_normal="${normal_reads}" \
            --reads_tumor="${tumor_reads}" \
            --output_vcf="${somatic_name}_deepsomatic.vcf.gz" \
            --output_gvcf="${somatic_name}_deepsomatic.g.vcf.gz" \
            --sample_name_tumor="${tumor_sample_id}" \
            --sample_name_normal="${normal_sample_id}" \
            --use_default_pon_filtering=true \
            --num_shards=$task.cpus \
            --logging_dir="${somatic_name}_logs.txt"
        """
}
*/

process VT_DECOMPOSE_GERMLINE {
    conda "bioconda::vt bioconda::tabix=0.2.6"

    publishDir "${params.outdir}/germline/HaplotypeCaller", mode: "copy"

    input:
        tuple val(sample_id), val(sample_type)
        tuple path(vcf), path(vcf_index)

    output:
        tuple path("${sample_id}_germline_decomp.vcf.gz"), path("${sample_id}_germline_decomp.vcf.gz.tbi"), emit: vcf
        tuple val(sample_id), val(sample_type), emit: sample_info
    script:
        """
        vt decompose -s $vcf -o "${sample_id}_germline_decomp.vcf.gz"

        tabix -p vcf "${sample_id}_germline_decomp.vcf.gz"
        """
}



process PHASE_VCF_INDEX {
    
    conda "bioconda::tabix=0.2.6"

    publishDir "${params.outdir}/variants/phased_variants/", mode: "copy"

    input:
        path(phased_vcf)
        tuple val(somatic_name), val(somatic_caller), val(tumor_name), val(normal_name)

    output:
        tuple path("${somatic_name}_phased_annotated.vcf.gz"), path("${somatic_name}_phased_annotated.vcf.gz.tbi"), emit: phased_vcf


    script:
        """
        bgzip -c $phased_vcf > ${somatic_name}_phased_annotated.vcf.gz
        
        tabix -p vcf ${somatic_name}_phased_annotated.vcf.gz
        """

}

process PHASE_VCF_VEP {

    //container "/c4/home/hwismer/0_execs/vep/vep.sif"
    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        path phased_vcf
        path reference_fa
        path vep_cache
        path(vep_plugins)

    output:
        path "phased_vcf_vep.vcf", emit: vcf

    script:
        """
        vep \
            --input_file $phased_vcf  \
            --output_file "phased_vcf_vep.vcf" \
            --format vcf --vcf --symbol --terms SO --tsl --biotype \
            --hgvs --fasta $reference_fa  \
            --offline --cache $vep_cache \
            --plugin Frameshift --plugin Wildtype \
            --pick \
            --fork ${task.cpus} \
            --dir_plugins $vep_plugins
            #[--transcript_version]
    
        """

}

process PHASE_VCF_RBPHASING {

    container "broadinstitute/gatk3:3.6-0" 

    input:
        path(combined_sorted_vcf)
        path(reference_fa)
        path(reference_index_dir)
        tuple val(tumor_sample), val(sample_type), path(tumor_reads), path(tumor_reads_index)

    output:
        path("phased.vcf"), emit:vcf
        

    script:

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T ReadBackedPhasing \
                -R $reference_index_dir/$reference_fa \
                -I $tumor_reads \
                --variant $combined_sorted_vcf \
                -L $combined_sorted_vcf \
                -o phased.vcf

        """




}

process PHASE_VCF_SORT_VCF {

    container 'broadinstitute/picard:3.4.0'
    
    input:
        path(somatic_germline_vcf_combined)
        path(reference_dict)

    output:
        path("combined_somatic_plus_germline.sorted.vcf"), emit: vcf

    script:
        """
        java -jar /usr/picard/picard.jar \
            SortVcf \
                -I $somatic_germline_vcf_combined \
                -O combined_somatic_plus_germline.sorted.vcf \
                -SD $reference_dict


        """



}


process PHASE_VCF_COMBINE_VARIANTS {

    //container "/c4/home/hwismer/pvac_pipeline/containers/broadinstitute-gatk3-3.6-0.sif"
    container "broadinstitute/gatk3:3.6-0"

    //publishDir "${params.outdir}/phase_vcf/", mode: "copy"

    input:
        path(tumor_only_vcf)
        path(tumor_only_vcf_index)
        tuple path(germline_vcf), path(germline_vcf_index)
        path(reference_fa)
        path(reference_index_dir)

    output:
        path("combined_somatic_plus_germline.vcf"), emit: combined_vcf

    script:
        """
        
        java -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
                -R $reference_index_dir/$reference_fa \
                --variant $germline_vcf \
                --variant $tumor_only_vcf \
                -o combined_somatic_plus_germline.vcf \
                --assumeIdenticalSamples

        """

}


process PHASE_VCF_SELECT_VARIANTS {

    //container "/c4/home/hwismer/pvac_pipeline/containers/gatk_4.6.1.0.sif"
    container "broadinstitute/gatk:4.6.1.0"
    
    //publishDir "${params.outdir}/phase_vcf/", mode: "copy"

    input:
        tuple val(sample_name), val(sample_type), path(sample_reads), path(sample_reads_index)
        tuple path(somatic_vcf), path(somatic_vcf_index)
        path(reference_fa)
        path(reference_index_dir)


    output:
        path("tumor_only.vcf.gz"), emit: vcf
        path("tumor_only.vcf.gz.tbi"), emit: vcf_index

    script:
        """
        gatk SelectVariants \
            -V $somatic_vcf \
            -R "${reference_index_dir}/${reference_fa}" \
            --sample-name $sample_name \
            -O tumor_only.vcf.gz

        gatk IndexFeatureFile \
            -I tumor_only.vcf.gz
        
        """

}

process POSTPROCESS_HAPLOTYPE_SCATTER {

    //container "/c4/home/hwismer/pvac_pipeline/containers/gatk_4.3.0.0.sif"

    container "broadinstitute/gatk:4.3.0.0"

    //publishDir "${params.outdir}/germline/HaplotypeCaller", mode: "copy"

    input:
        tuple val(sample_id), val(sample_type)
        path(vcfs)
        path(reference_fa)
        path(reference_index_dir)
        path(hapmap)
        path hapmap_index
        path(mills)
        path mills_index


    output:
        tuple val(sample_id), val(sample_type), emit: germline_sample_info
        tuple path("${sample_id}_germline_filtered.vcf.gz"), path("${sample_id}_germline_filtered.vcf.gz.tbi"), emit: germline_vcf

    script:

        def sorted_vcfs = vcfs.sort { vcf ->
            def matcher = vcf.name =~ /(\d+)-scattered/
                matcher.find() ? matcher.group(1).toInteger() : 0
        }

        def vcf_as_input = sorted_vcfs.collect { vcf ->
            "--INPUT ${vcf}"
        }.join(' ')

        """
        gatk GatherVcfs \
            $vcf_as_input \
            -O "${sample_id}_germline_merged.vcf.gz"
    
        gatk IndexFeatureFile \
            -I "${sample_id}_germline_merged.vcf.gz"

        gatk FilterVariantTranches \
            -V "${sample_id}_germline_merged.vcf.gz" \
            --resource $hapmap \
            --resource $mills \
            --info-key CNN_1D \
            --snp-tranche 99.95 \
            --indel-tranche 99.4 \
            -O "${sample_id}_germline_filtered.vcf.gz"

        gatk IndexFeatureFile \
            -I "${sample_id}_germline_filtered.vcf.gz"

        """
}





process HAPLOTYPE_CALLER_SCATTER {
    
    maxForks 50
    cpus 2
    memory "16GB"

    //container "/c4/home/hwismer/pvac_pipeline/containers/gatk_4.3.0.0.sif"

    container "broadinstitute/gatk:4.3.0.0"


    // publishDir "${params.outdir}/germline/scatter/", mode: "copy"


    input:
        tuple val(sample_id), val(sample_type), path(sample_reads), path(sample_reads_index)
        each path(interval_shard)
        path reference_fa
        path reference_index_dir

    output:
        tuple val(sample_id), val(sample_type), emit: sample_info
        path("${sample_id}_${interval_shard}_CNN.vcf.gz"), emit: vcf

    script:

        """
        gatk HaplotypeCaller \
            -R $reference_index_dir/$reference_fa \
            -I $sample_reads \
            -L $interval_shard \
            -O "${sample_id}_${interval_shard}.vcf.gz" \
            -ERC NONE

        gatk CNNScoreVariants \
            -V "${sample_id}_${interval_shard}.vcf.gz" \
            -L $interval_shard \
            -R $reference_index_dir/$reference_fa \
            -O "${sample_id}_${interval_shard}_CNN.vcf.gz" \


        """
     

}


process ANNOTATE_VCF_EXPRESSION {

    container "griffithlab/vatools:5.2.0"

    publishDir "${params.outdir}/coverage/", mode: "copy"

    input:
        path(coverage_annotated_vcf)
        tuple val(sample_id), val(sample_type), path(kallisto_quant_dir)
        tuple val(somatic_name), val(somatic_caller), val(tumor_name), val(normal_name)

    output:
        path("${somatic_name}_cov_expr_annotated.vcf")

    script:
        """
        vcf-expression-annotator \
            $coverage_annotated_vcf \
            -s $sample_id \
            "${kallisto_quant_dir}/abundance.tsv" \
            kallisto transcript \
            -o "${somatic_name}_cov_expr_annotated.vcf"

        """

}

process ANNOTATE_VCF_COVERAGE {
    
    container "griffithlab/vatools:5.2.0"

    publishDir "${params.outdir}/coverage/", mode: "copy"

    input:
        path(vcf)
        tuple val(sample_name), val(sample_type), path(brc_indels), path(brc_snvs)
        tuple val(somatic_name), val(somatic_caller), val(tumor_name), val(normal_name)

    output:
        path("${somatic_name}_annotated.vcf")

    script:
        """
        vcf-readcount-annotator \
            $vcf \
            $brc_snvs \
            RNA \
            -s $sample_name \
            -t snv \
            -o "${somatic_name}_snv_annotated.vcf"

        vcf-readcount-annotator \
            "${somatic_name}_snv_annotated.vcf" \
            $brc_indels \
            RNA \
            -s $sample_name \
            -t indel \
            -o ${somatic_name}_annotated.vcf

        """

}

process BAMREADCOUNT {
    //container "/c4/home/hwismer/pvac_pipeline/containers/bam_readcount_helper-cwl_1.2.1.sif"

    container "mgibio/bam_readcount_helper-cwl:1.2.1"

    publishDir "${params.outdir}/coverage/", mode: "copy"

    input:
        path(vt_vcf)
        tuple val(somatic_name), val(somatic_caller), val(tumor_name), val(normal_name)
        path reference_fa
        tuple val(sample_name), val(sample_type), path(sample_bam), path(sample_bam_index)

    output:
        tuple val(sample_name), val(sample_type), path("${sample_name}_bamrc_helper/${sample_name}_bam_readcount_indel.tsv"), path("${sample_name}_bamrc_helper/${sample_name}_bam_readcount_snv.tsv"), emit: brc_files
        tuple val(somatic_name), val(somatic_caller), val(tumor_name), val(normal_name), emit: somatic_info
    script:
        """
        mkdir ${sample_name}_bamrc_helper
        bam_readcount_helper.py \
            $vt_vcf \
            $sample_name \
            $reference_fa \
            $sample_bam \
            NOPREFIX \
            ${sample_name}_bamrc_helper
        """

}

process VT_DECOMPOSE {
    conda "bioconda::vt"

    //publishDir "${params.outdir}/vcf_decompose/", mode: "copy"

    input:
        tuple path(somatic_vcf), path(somatic_vcf_index)
        tuple val(somatic_name), val(somatic_caller), val(tumor_name), val(normal_name)

    output:
        path("${somatic_name}_vt_decomp.vcf.gz"), emit: vcf
        tuple val(somatic_name), val(somatic_caller), val(tumor_name), val(normal_name), emit: caller_info

    script:
        """
        vt decompose -s $somatic_vcf -o "${somatic_name}_vt_decomp.vcf.gz"
        """
}


process VEP_ANNOTATE {
    
    container "ensemblorg/ensembl-vep:release_115.0"

    //publishDir "${params.outdir}/vep_annotated/", mode: "copy"

    input: 
        path(somatic_vcf)
        tuple val(somatic_name), val(somatic_caller), val(tumor_name), val(normal_name)
        path reference_fa
        path vep_cache
        path vep_plugins

    output:
        path "${somatic_name}_vep.vcf", emit: vcf
        tuple val(somatic_name), val(somatic_caller), val(tumor_name), val(normal_name), emit: caller_info

    script:
        """
        vep \
            --input_file $somatic_vcf  \
            --output_file "${somatic_name}_vep.vcf" \
            --format vcf --vcf --symbol --terms SO --tsl --biotype \
            --hgvs --fasta $reference_fa  \
            --offline --cache $vep_cache \
            --plugin Frameshift --plugin Wildtype \
            --pick \
            --dir_plugins $vep_plugins
            #[--transcript_version]
        """
}


process FILTER_VCF {

    conda "bioconda::bcftools=1.22"

    input:
        path(vcf)
        tuple val(somatic_name)

    output:
        path filtered_vcf, emit:filtered_vcf

    script:
    """
        bcftools view -f PASS $vcf -o passing_variants.vcf

    """

}

process INDEX_FINAL_VCF {

    conda "bioconda::tabix=0.2.6"

    publishDir "${params.outdir}/variants", mode: "copy"

    input:
        path(vcf)

    output:
        tuple path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
        """
        bgzip $vcf
        tabix -p vcf "${vcf}.gz"
        """

}


/*
process POSTPROCESS_STRELKA {

    publishDir "${params.outdir}/strelka", mode: "copy"

    input:
        tuple val(somatic_name), val(strelka_caller), val(tumor_sample_id), val(normal_sample_id), path(strelka_results_dir)
        tuple val(somatic_name), val(manta_caller), val(tumor_sample_id), val(normal_sample_id), path(manta_results_dir)


}
*/


process STRELKA {

   //container "/c4/home/hwismer/pvac_pipeline/containers/strelka_manta/strelka2-manta_latest.sif"
   container 'quay.io/wtsicgp/strelka2-manta'

   publishDir "${params.outdir}/somatic/", mode: "copy"

   input:
        tuple val(tumor_sample_id), val(tumor_type),  path(tumor_reads), path(tumor_reads_index)
        tuple val(normal_sample_id), val(normal_type), path(normal_reads), path(tumor_reads_index)
        val(somatic_name)
        path reference_fa
        path reference_index_dir

    output:
        tuple val(somatic_name), val("Strelka"), val(tumor_sample_id), val(normal_sample_id), path("./strelka/results/"), emit: strelka
        tuple val(somatic_name), val("Manta"), val(tumor_sample_id), val(normal_sample_id), path("./manta/results/"), emit: manta

    script:
        """
        configManta.py \
            --normalBam $normal_reads \
            --tumorBam $tumor_reads \
            --referenceFasta "${reference_index_dir}/${reference_fa}" \
            --runDir ./manta/ \
            --exome

        ./manta/runWorkflow.py -j $task.cpus

        configureStrelkaSomaticWorkflow.py \
            --normalBam $normal_reads \
            --tumorBam $tumor_reads \
            --referenceFasta "${reference_index_dir}/${reference_fa}" \
            --indelCandidates ./manta/results/variants/candidateSmallIndels.vcf.gz \
            --runDir "./strelka" \
            --exome

        ./strelka/runWorkflow.py -m local -j $task.cpus

        """

}

process POSTPROCESS_MUTECT2_SCATTER {

    //container "/c4/home/hwismer/pvac_pipeline/containers/gatk_4.6.1.0.sif"
    container "broadinstitute/gatk:4.6.1.0"
    

    publishDir "${params.outdir}/somatic/mutect2", mode: "copy"


    input:
        tuple val(somatic_name), val(somatic_caller), val(tumor_sample), val(normal_sample)
        path(vcfs)
        path(f1r2s)
        path(stats)
        tuple val(sample_id), val(sample_type), path (tumor_pileups)
        tuple val(sample_id), val(sample_type), path (normal_pileups)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple path("${somatic_name}_mutect_postprocessed.vcf.gz"), path("${somatic_name}_mutect_postprocessed.vcf.gz.tbi"), emit: vcf
        tuple val(somatic_name), val(somatic_caller), val(tumor_sample), val(normal_sample), emit: caller_info

    script:

        def sorted_vcfs = vcfs.sort { vcf ->
            def matcher = vcf.name =~ /(\d+)-scattered/
                matcher.find() ? matcher.group(1).toInteger() : 0
        }

        def vcf_as_input = sorted_vcfs.collect { vcf ->
            "--INPUT ${vcf}"
        }.join(' ')

        def f1r2_as_input = f1r2s.collect { f1r2 ->
            "-I ${f1r2}"
        }.join(' ')

        def stat_as_input = stats.collect {stat ->
            "-stats ${stat}"
        }.join(' ')

        """
        
        gatk CalculateContamination \
            -I $tumor_pileups \
            -matched $normal_pileups \
            -O contamination.table


        gatk GatherVcfs \
            $vcf_as_input \
            -O "${somatic_name}_merged.vcf"


        gatk LearnReadOrientationModel \
            $f1r2_as_input \
            -O "${somatic_name}_orientmodel.tar.gz"

        gatk MergeMutectStats \
            $stat_as_input \
            -O "${somatic_name}_merged.stats"

        gatk FilterMutectCalls \
            -R "${reference_index_dir}/${reference_fa}" \
            -V "${somatic_name}_merged.vcf" \
            --orientation-bias-artifact-priors "${somatic_name}_orientmodel.tar.gz" \
            --contamination-table contamination.table \
            -stats "${somatic_name}_merged.stats" \
            -O "${somatic_name}_mutect_postprocessed.vcf.gz"

        """
}


process MUTECT2_SCATTER {

    maxForks 50
    cpus 2
    memory "16GB"

    //container "/c4/home/hwismer/pvac_pipeline/containers/gatk_4.6.1.0.sif"
    container "broadinstitute/gatk:4.6.1.0"
    

    // publishDir "${params.outdir}/somatic/mutect2", mode: "copy"


    input:
        tuple val(tumor_sample_id), val(tumor_type),  path(tumor_reads), path(tumor_reads_index)
        tuple val(normal_sample_id), val(normal_type), path(normal_reads), path(normal_reads_index)

        val(somatic_name)
        each path(interval_shard)
        path reference_fa
        path reference_index_dir
        path known_sites
        path known_sites_dir
        path pon
        path pon_index_dir

    output:
        tuple val(somatic_name), val("Mutect2"), val(tumor_sample_id), val(normal_sample_id), emit: caller_info
        path("*.vcf.gz"), emit: vcf
        path("*vcf.gz.tbi"), emit: vcf_index
        path("*_mutect_f1r2.tar.gz"), emit: f1r2
        path("*.stats"), emit: stat

    script:

        """
        gatk Mutect2 \
            -R "${reference_index_dir}/${reference_fa}" \
            -I ${tumor_reads} \
            -I ${normal_reads} \
            -normal ${normal_sample_id} \
            --germline-resource "${known_sites}" \
            --panel-of-normals "${pon}" \
            --f1r2-tar-gz "${somatic_name}_${interval_shard}_mutect_f1r2.tar.gz" \
            -L $interval_shard \
            -O "${somatic_name}_${interval_shard}_mutect.vcf.gz"

        """
}


process SPLIT_INTERVALS {

    tag "Split intervals for ${params.scatter_count} shards"
    
    // Use the GATK docker image for consistency
    //container "/c4/home/hwismer/pvac_pipeline/containers/gatk_4.6.1.0.sif"    
    container "broadinstitute/gatk:4.6.1.0"
    

    input:
        path reference_fa
        path reference_index_dir
        path intervals_file
        val scatter_count
    
    output:
        path "*-scattered.interval_list", emit: interval_shards

    script:
    """
    gatk SplitIntervals \
        -R $reference_index_dir/$reference_fa \
        -L ${intervals_file} \
        --scatter-count $scatter_count \
        -O . 
        
    """
}

process KALLISTO_QUANT {
    conda "bioconda::kallisto=0.51.1"

    publishDir "${params.outdir}/rnaseq/kallisto/", mode: "copy"

    input:
        tuple val(sample_id), path(fastq_1), path(fastq_2), val(sample_type)
        path kallisto_index

    output:
        tuple val(sample_id), val(sample_type), path("${sample_id}_kallisto")

    script:
        """
        kallisto quant -i $kallisto_index -o ${sample_id}_kallisto -t ${task.cpus} $fastq_1 $fastq_2  
        """
}

process PREPROCESS_BAM {
    
    //container "/c4/home/hwismer/pvac_pipeline/containers/gatk_4.6.1.0.sif"
    
    container "broadinstitute/gatk:4.6.1.0"

    publishDir "${params.outdir}/preprocess_bam/", mode: "copy"

    input:
        tuple val(sample_id), val(sample_type), path(reads), path(reads_index)
        path reference_fa
        path reference_index_dir
        path known_sites
        path known_sites_index

    output:
        tuple val(sample_id), val(sample_type), path("${sample_id}/${sample_id}_bqsr.bam"), path("${sample_id}/${sample_id}_bqsr.bai"), emit: data
        tuple val(sample_id), val(sample_type), path ("${sample_id}/${sample_id}_pileups.table"), emit: pileups
        path("${sample_id}/${sample_id}_recal_table.table"), emit: recal_table
        path("${sample_id}/${sample_id}_dedup.bam"), emit: dedup_bam

    script:
        """
        gatk MarkDuplicatesSpark \
            -I $reads \
            -O "${sample_id}/${sample_id}_dedup.bam" \
            --tmp-dir "\${PWD}"

        gatk BaseRecalibrator \
            -I "${sample_id}/${sample_id}_dedup.bam" \
            -O "${sample_id}/${sample_id}_recal_table.table" \
            -R "${reference_index_dir}/${reference_fa}" \
            --known-sites "$known_sites"
        
        gatk ApplyBQSR \
            -R $reference_index_dir/$reference_fa \
            -I "${sample_id}/${sample_id}_dedup.bam" \
            --bqsr-recal-file "${sample_id}/${sample_id}_recal_table.table" \
            -O "${sample_id}/${sample_id}_bqsr.bam" \
            --create-output-bam-index

        gatk GetPileupSummaries \
            -I "${sample_id}/${sample_id}_bqsr.bam" \
            -V "${known_sites}" \
            -L "${known_sites}" \
            -O "${sample_id}/${sample_id}_pileups.table"
        """
}

process BWA_MAP {

    //container "/c4/home/hwismer/pvac_pipeline/containers/custom_ipi_utils.sif"

    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.22.1"

    publishDir "${params.outdir}/alignment/bwa", mode: "copy"

    input:
        tuple val(sample_id), val(sample_type), path(fastq1), path(fastq2)
        path reference_fa
        path reference_index_dir

    output:
        tuple val(sample_id), val(sample_type), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: data
    
    script:
    """
    NEW_RG="@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA"
    bwa mem -M -t $task.cpus -R \$NEW_RG ${reference_index_dir}/${reference_fa} ${fastq1} ${fastq2} \
        | samtools sort --threads $task.cpus -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    """
}


process FASTP {
    
    //container "/c4/home/hwismer/pvac_pipeline/containers/custom_ipi_utils.sif"

    conda "bioconda::fastp=1.0.1"

    input:
        tuple val(sample_id), path(read1), path(read2), val(sample_type)

    output:
        tuple val(sample_id), val(sample_type), path("${sample_id}_R1_*fastq.gz*"), path("${sample_id}_R2_*fastq.gz*"), emit: fastqs
        path("${sample_id}_*{html,json}*"), emit: reports
    
    tag "FastP on ${sample_id}"

    publishDir "${params.outdir}/qc/${sample_id}_fastp/", mode: 'copy'

    script:
        """
        fastp --thread $task.cpus \
              -i $read1 \
              -I $read2 \
              -o "${sample_id}_R1_fastp.fastq.gz" \
              -O "${sample_id}_R2_fastp.fastq.gz" \
              -R "${sample_id}_fastp_report" \
              -h "${sample_id}_fastp_report.html" \
              -j "${sample_id}_fastp_report.json" \

        """
}


// Workflow block

workflow {


    // INPUT PARSING

    // Read in input data consiting of:
    // 1) Tumor FastQ WES/WGS Pairs
    // 2) Normal FastQ WES/WGS Pairs
    // 3) Tumor FastQ RNAseq Pairs

    tumor_fastqs = channel.fromFilePairs( params.tumor_fastq_dir, 
                                           checkIfExists: true,
                                           flat: true)
                            .merge(Channel.of("Tumor"))

    normal_fastqs = channel.fromFilePairs( params.normal_fastq_dir,
                                            checkIfExists: true,
                                            flat: true)
                                .merge(Channel.of("Normal"))

    tumor_rna_fastqs = channel.fromFilePairs(params.tumor_fastq_rna_dir, 
                                              checkIfExists: true, 
                                              flat: true)
                            .merge(Channel.of("Tumor RNA"))


    // Merge all fastq pair channels to input into fastp qc and BWA alignment
    fastqs_ch = tumor_fastqs.concat(normal_fastqs).concat(tumor_rna_fastqs)
    fastp = FASTP(fastqs_ch)
    mapped = BWA_MAP(fastp.fastqs, params.reference_fa, params.reference_index_dir)
    
    // Split BWA aligned channel into WES/WGS and RNA channels
    mapped.data.branch { output ->
            tumor_normal: output[1] == "Normal" || output[1] == "Tumor"
            tumor_rna: output[1] == "Tumor RNA"
    }.set { mapped_branch }
    tumor_normal_mapped = mapped_branch.tumor_normal
    tumor_rna = mapped_branch.tumor_rna

    

    // PREPROCESSING



    // Preprocess WES/WGS samples according to GATK standard (BQSR, pileups etc.)
    preprc = PREPROCESS_BAM(tumor_normal_mapped,
                            params.reference_fa, params.reference_index_dir, 
                            params.common_germline, params.common_germline_index)

    // Branch preprocessed data into tumor and normal channels
    preprc.data.branch { output ->
            normal: output[1] == "Normal"
            tumor: output[1] == "Tumor"
    }.set { tumor_normal }

    
    // Branch preprocessed pileup data into tumor and normal channels for mutect later
    preprc.pileups.branch { output ->
        normal: output[1] == "Normal"
        tumor: output[1] == "Tumor"
    }.set { pileups }

    


    // TRANSCRIPT ABUNDANCE ESTIMATION

    // Run kallisto to get transcript abundance estimates to later annotate VCF with
    kallisto = KALLISTO_QUANT(tumor_rna_fastqs,
                              params.kallisto_index)


    

    // INTERVAL CREATION FOR SOMATIC CALLERS

    intervals = SPLIT_INTERVALS(params.reference_fa, params.reference_index_dir,
                                params.intervals_file, params.scatter_count)




    // GERMLINE VARIANT CALLING

    // Run HaplotypeCaller with scatter/gather approach and preprocess with CNNScoreVariants and FilterVariantTranches
    // For later use in phasing the somatic VCF with proximal variants

    germline = HAPLOTYPE_CALLER_SCATTER(tumor_normal.tumor,
                                        intervals.flatten(),
                                        params.reference_fa, params.reference_index_dir)

    germline_postprocess = POSTPROCESS_HAPLOTYPE_SCATTER(germline.sample_info.first(), 
                                                        germline.vcf.collect(),
                                                        params.reference_fa,
                                                        params.reference_index_dir,
                                                        params.hapmap,
                                                        params.hapmap_index,
                                                        params.mills,
                                                        params.mills_index)

    // Use vt decompose on germline calls
    vt_germline = VT_DECOMPOSE_GERMLINE(germline_postprocess.germline_sample_info, germline_postprocess.germline_vcf)




    // SOMATIC VARIANT CALLING
    
    // Run Mutect2 with scatter/gather approach
    mutect_scattered = MUTECT2_SCATTER(tumor_normal.tumor, tumor_normal.normal, params.somatic_name,
                    intervals.flatten(),
                    params.reference_fa, params.reference_index_dir,
                    params.known_sites, params.known_sites_index,
                    params.pon, params.pon_index)

    // Postprocess Mutect2 Output using pileups, stats, etc.
    mutect_postprocess = POSTPROCESS_MUTECT2_SCATTER(mutect_scattered.caller_info.first(),
                                                    mutect_scattered.vcf.collect(),
                                                    mutect_scattered.f1r2.collect(),
                                                    mutect_scattered.stat.collect(),
                                                    pileups.tumor,
                                                    pileups.normal,
                                                    params.reference_fa,
                                                    params.reference_index_dir)
    
    // Run Strelka2
    strelka = STRELKA(tumor_normal.tumor, tumor_normal.normal, 
            params.somatic_name,
            params.reference_fa, params.reference_index_dir)


    
    

    // SOMATIC VCF PREPROCESSING


    //filtered_vcf = FILTER_SOMATIC_NONPASSING(mutect_postprocess.vcf, mutect_postprocess.caller_info)

    // Use vt decompose to get simplest/standardized variant representation for vcf compatibility
    vt = VT_DECOMPOSE(mutect_postprocess.vcf, mutect_postprocess.caller_info)
    
    // Annotate VCF with VEP
    vep = VEP_ANNOTATE(vt.vcf,
                       vt.caller_info,
                       params.reference_fa,
                       params.vep_cache,
                       params.vep_plugins)

    // Add Read Coverage to VCF
    // Currently only for Tumor RNA as Mutect2 handles the DNA already
    bamreadcount = BAMREADCOUNT(vep.vcf,
                                vep.caller_info,
                                params.reference_fa,
                                tumor_rna)

    vcf_annotated_coverage = ANNOTATE_VCF_COVERAGE(vep.vcf, bamreadcount.brc_files, vep.caller_info)

    vcf_annotated_expression = ANNOTATE_VCF_EXPRESSION(vcf_annotated_coverage, kallisto, vep.caller_info)

    vcf_final = INDEX_FINAL_VCF(vcf_annotated_expression)
    



    // PERFORM VCF PHASING USING GERMLINE CALLS

    // Create Tumor-Only VCF From Final Somatic VCF (vcf_final)
    vcf_phase_select_variants = PHASE_VCF_SELECT_VARIANTS(tumor_normal.tumor,
                           vcf_final,
                           params.reference_fa,
                           params.reference_index_dir)

    // Combine Tumor-Only VCF with germline variants
    vcf_phase_combine = PHASE_VCF_COMBINE_VARIANTS(vcf_phase_select_variants.vcf,
                                                   vcf_phase_select_variants.vcf_index,
                                                   vt_germline.vcf,
                                                   params.reference_fa,
                                                   params.reference_index_dir)

    // Sort combined VCF
    vcf_phase_sort = PHASE_VCF_SORT_VCF(vcf_phase_combine.combined_vcf, params.reference_dict)

    // Call GATK ReadBackedPhasing (required gatk 3.6.0) to phase VCF
    vcf_phase_rbphase = PHASE_VCF_RBPHASING(vcf_phase_sort.vcf, 
                                            params.reference_fa, params.reference_index_dir,
                                            tumor_normal.tumor)
    /// VEP Annotation phased vcf
    vcf_phase_vep = PHASE_VCF_VEP(vcf_phase_rbphase.vcf,
                                            params.reference_fa, params.vep_cache, params.vep_plugins)        
    
    // Zip and index phased vcf
    vcf_phased = PHASE_VCF_INDEX(vcf_phase_vep.vcf, vep.caller_info)

        
        
    vcf_final.view()

    vcf_phased.view()



}
