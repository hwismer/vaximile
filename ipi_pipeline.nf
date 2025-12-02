// Default parameter input

params.outdir = "./vax_pipeline_out/"

process DEEPSOMATIC {

    container "/c4/home/hwismer/pvac_pipeline/containers/deepsomatic_1.9.0.sif"
    
    publishDir "${params.outdir}/deepsomatic/", mode: "copy"

    input:
        tuple val(tumor_meta), path(tumor_bam), path(tumor_bam_index)
        tuple val(normal_meta), path(normal_bam), path(normal_bam_index)
        val somatic_name
        path reference_fa
        path reference_index_dir

    output:
        tuple val(somatic_meta), path("${somatic_name}_deepsomatic.vcf.gz"), path("${somatic_name}_deepsomatic.vcf.gz.tbi"), emit:deepsom_vcf

    script:

        somatic_meta = [
            somatic_name: somatic_name,
            somatic_caller: 'deepsomatic',
            tumor_metamap: tumor_meta,
            normal_metamap: normal_meta
        ]

        """
        
        mkdir -p \$TMPDIR

        run_deepsomatic \
            --model_type=${tumor_meta.sequencing_type} \
            --ref="${reference_fa}" \
            --reads_normal="${normal_bam}" \
            --reads_tumor="${tumor_bam}" \
            --output_vcf="${somatic_name}_deepsomatic.vcf.gz" \
            --output_gvcf="${somatic_name}_deepsomatic.g.vcf.gz" \
            --sample_name_tumor="${tumor_meta.sample_name}" \
            --sample_name_normal="${normal_meta.sample_name}" \
            --use_default_pon_filtering=true \
            --num_shards=$task.cpus \
            --logging_dir=./logs \
            --intermediate_results_dir ./intermediate_dir

        """
}

process VT_DECOMPOSE_GERMLINE {
    


    conda "bioconda::vt bioconda::tabix=0.2.6"

    publishDir "${params.outdir}/germline/HaplotypeCaller", mode: "copy"

    input:
        tuple val(meta), path(germline_vcf), path(germline_vcf_index)

    output:
        tuple val(meta), path("${meta.sample_name}_germline.vcf.gz"), path("${meta.sample_name}_germline.vcf.gz.tbi"), emit: germline_vcf
    script:
        """
        vt decompose -s $germline_vcf -o "${meta.sample_name}_germline.vcf.gz"

        tabix -p vcf "${meta.sample_name}_germline.vcf.gz"
        """
}




process PHASE_VCF_INDEX {
    
    conda "bioconda::tabix=0.2.6"

    publishDir "${params.outdir}/variants/phased_variants/", mode: "copy"

    input:
        tuple val(meta), path(phased_vcf)

    output:
        tuple val(meta), path("${meta.somatic_name}_phased_annotated.vcf.gz"), path("${meta.somatic_name}_phased_annotated.vcf.gz.tbi"), emit: phased_vcf


    script:
        """
        echo ${meta.somatic_name}
        bgzip -c $phased_vcf > ${meta.somatic_name}_phased_annotated.vcf.gz
        
        tabix -p vcf ${meta.somatic_name}_phased_annotated.vcf.gz
        """

}

process PHASE_VCF_VEP {

    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(meta), path(phased_vcf)
        path reference_fa
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("phased_vcf_vep.vcf"), emit: vcf

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
        tuple val(meta), path(combined_sorted_vcf)
        path(reference_fa)
        path(reference_index_dir)
        tuple val(meta), path(normal_reads), path(normal_reads_index)

    output:
        tuple val(meta), path("phased.vcf")

    script:

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T ReadBackedPhasing \
                -R $reference_fa \
                -I $normal_reads \
                --variant $combined_sorted_vcf \
                -L $combined_sorted_vcf \
                -o phased.vcf

        """
}

process PHASE_VCF_RENAME {
    container "biocontainers/bcftools:v1.9-1-deb_cv1"

    input:
        tuple val(meta), path(germline_vcf), path(germline_vcf_index)
        tuple val(tumor_meta), path(tumor_reads), path(tumor_reads_index)

    output:
        tuple val(meta), path("vt_germline_rename.vcf.gz" ), path("vt_germline_rename.vcf.gz.tbi")

    script:
        """
        echo ${tumor_meta.sample_name} > new_names.txt
        
        bcftools reheader \
            --samples new_names.txt \
            --output vt_germline_rename.vcf.gz \
            $germline_vcf 

        bcftools index -t vt_germline_rename.vcf.gz
        """


}

process PHASE_VCF_SORT_VCF {

    container 'broadinstitute/picard:3.4.0'
    
    input:
        tuple val(meta), path(combined_vcf)
        path(reference_dict)

    output:
        tuple val(meta), path("combined.sorted.vcf"), emit: vcf

    script:
        """
        java -jar /usr/picard/picard.jar \
            SortVcf \
                -I $combined_vcf \
                -O combined.sorted.vcf \
                -SD $reference_dict

        """


}


process PHASE_VCF_COMBINE_VARIANTS {

    container "broadinstitute/gatk3:3.6-0"

    input:
        tuple val(meta), path(subset_vcf), path(subset_vcf_index)
        tuple val(germline_meta), path(germline_vcf), path(germline_vcf_index)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(meta), path("combined_somatic_plus_germline.vcf"), emit: combined_vcf

    script:
        """
        
        java -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
                -R $reference_fa \
                --variant $germline_vcf \
                --variant $subset_vcf \
                -o combined_somatic_plus_germline.vcf \
                --assumeIdenticalSamples

        """

}


process PHASE_VCF_SELECT_VARIANTS {


    container "broadinstitute/gatk:4.6.1.0"
    
    input:
        tuple val(sample_meta), path(sample_reads), path(sample_reads_index)
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(somatic_meta), path("tumor_only.vcf.gz"), path("tumor_only.vcf.gz.tbi"), emit:vcf

    script:
        """
        gatk SelectVariants \
            -V $somatic_vcf \
            -R "${reference_fa}" \
            --sample-name ${sample_meta.sample_name} \
            -O tumor_only.vcf.gz

        gatk IndexFeatureFile \
            -I tumor_only.vcf.gz
        
        """

}

process POSTPROCESS_HAPLOTYPE_SCATTER {

    cpus 4
    memory "64GB"

    container "broadinstitute/gatk:4.3.0.0"

    input:
        tuple val(meta), path(vcfs)
        path reference_fa
        path reference_index_dir
        path hapmap
        path hapmap_index
        path mills
        path mills_index


    output:
        tuple val(meta), path("${meta.sample_name}_germline_filtered.vcf.gz"), path("${meta.sample_name}_germline_filtered.vcf.gz.tbi"), emit: germline_vcf

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
            -O "${meta.sample_name}_germline_merged.vcf.gz"
    
        gatk IndexFeatureFile \
            -I "${meta.sample_name}_germline_merged.vcf.gz"

        gatk FilterVariantTranches \
            -V "${meta.sample_name}_germline_merged.vcf.gz" \
            --resource $hapmap \
            --resource $mills \
            --info-key CNN_1D \
            --snp-tranche 99.95 \
            --indel-tranche 99.4 \
            -O "${meta.sample_name}_germline_filtered.vcf.gz"

        gatk IndexFeatureFile \
            -I "${meta.sample_name}_germline_filtered.vcf.gz"

        """
}





process HAPLOTYPE_CALLER_SCATTER {
    
    cpus 2
    memory "16GB"

    container "broadinstitute/gatk:4.3.0.0"

    input:
        tuple val(meta), path(sample_reads), path(sample_reads_index)
        each path(interval_shard)
        path reference_fa
        path reference_index_dir

    output:
        tuple val(meta), path("${meta.sample_name}_${interval_shard}_CNN.vcf.gz"), emit: vcfs

    script:

        """
        gatk --java-options "-Xmx12G" HaplotypeCaller \
            -R $reference_fa \
            -I $sample_reads \
            -L $interval_shard \
            -O "${meta.sample_name}_${interval_shard}.vcf.gz" \
            -ERC NONE

        gatk --java-options "-Xmx12G" CNNScoreVariants \
            -V "${meta.sample_name}_${interval_shard}.vcf.gz" \
            -L $interval_shard \
            -R $reference_fa \
            -O "${meta.sample_name}_${interval_shard}_CNN.vcf.gz" \


        """
     

}


process ANNOTATE_VCF_EXPRESSION {

    container "griffithlab/vatools:5.2.0"

    publishDir "${params.outdir}/coverage/", mode: "copy"

    input:
        tuple val(somatic_meta), path(vcf)
        tuple val(sample_meta), path(kallisto_dir)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_cov_expr_annotated.vcf")

    script:
        """
        vcf-expression-annotator \
            $vcf \
            -s ${sample_meta.sample_name} \
            "${kallisto_dir}/abundance.tsv" \
            kallisto transcript \
            -o "${somatic_meta.somatic_name}_cov_expr_annotated.vcf"

        """

}

process ANNOTATE_VCF_COVERAGE {
    
    container "griffithlab/vatools:5.2.0"

    publishDir "${params.outdir}/coverage/", mode: "copy"

    input:
        tuple val(somatic_meta), path(vcf), val(sample_meta),
        path(brc_indels), path(brc_snvs)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_annotated.vcf")

    script:
        """
        vcf-readcount-annotator \
            $vcf \
            $brc_snvs \
            RNA \
            -s ${sample_meta.sample_name} \
            -t snv \
            -o "${somatic_meta.somatic_name}_snv_annotated.vcf"

        vcf-readcount-annotator \
            "${somatic_meta.somatic_name}_snv_annotated.vcf" \
            $brc_indels \
            RNA \
            -s ${sample_meta.sample_name} \
            -t indel \
            -o ${somatic_meta.somatic_name}_annotated.vcf

        """

}

process BAMREADCOUNT {

    container "mgibio/bam_readcount_helper-cwl:1.2.1"

    publishDir "${params.outdir}/coverage/", mode: "copy"

    input:
        tuple val(somatic_meta), path(vt_vcf)
        path reference_fa
        tuple val(sample_meta), path(sample_bam), path(sample_bam_index)

    output:
        tuple val(somatic_meta), path(vt_vcf), val(sample_meta),
        path("${sample_meta.sample_name}_bamrc_helper/${sample_meta.sample_name}_bam_readcount_indel.tsv"), 
        path("${sample_meta.sample_name}_bamrc_helper/${sample_meta.sample_name}_bam_readcount_snv.tsv"), emit: brc_files
    script:
        """
        mkdir ${sample_meta.sample_name}_bamrc_helper
        bam_readcount_helper.py \
            $vt_vcf \
            ${sample_meta.sample_name} \
            $reference_fa \
            $sample_bam \
            NOPREFIX \
            ${sample_meta.sample_name}_bamrc_helper
        """

}

process VT_SOMATIC_POSTPROCESS {
    
    conda "bioconda::vt bioconda::samtools"
    publishDir "${params.outdir}/somatic/${meta.somatic_caller}/"
    input:
        tuple val(meta), path(somatic_vcf), path(somatic_vcf_index)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(meta), path("${meta.somatic_name}_${meta.somatic_caller}_variants.vcf.gz"), path("${meta.somatic_name}_${meta.somatic_caller}_variants.vcf.gz.tbi"), emit: vt_vcf

    script:
        """
        vt decompose -s $somatic_vcf -o "${meta.somatic_caller}_decomp.vcf.gz"
        vt index "${meta.somatic_caller}_decomp.vcf.gz"
        vt normalize "${meta.somatic_caller}_decomp.vcf.gz" -r $reference_fa -o decomp_norm.vcf.gz
        vt index decomp_norm.vcf.gz
        vt uniq decomp_norm.vcf.gz -o "${meta.somatic_name}_${meta.somatic_caller}_variants.vcf.gz"
        vt index "${meta.somatic_name}_${meta.somatic_caller}_variants.vcf.gz"
        """
}


process VEP_ANNOTATE {

    
    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(meta), path(somatic_vcf)
        path reference_fa
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("${meta.somatic_name}_vep.vcf"), emit: vep_vcf

    script:
        """
        vep \
            --input_file $somatic_vcf  \
            --output_file "${meta.somatic_name}_vep.vcf" \
            --format vcf --vcf --symbol --terms SO --tsl --biotype \
            --hgvs --fasta $reference_fa  \
            --offline --cache $vep_cache \
            --plugin Frameshift --plugin Wildtype \
            --pick \
            --dir_plugins $vep_plugins
            #[--transcript_version]
        """
}



process INDEX_FINAL_VCF {

    cpus 2
    memory "32GB"

    conda "bioconda::tabix=0.2.6"

    publishDir "${params.outdir}/variants", mode: "copy"

    input:
        tuple val(meta), path(vcf)
        

    output:
        tuple val(meta), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
        """
        bgzip $vcf
        tabix -p vcf "${vcf}.gz"
        """

}


process MERGE_SOMATIC_VCFS {

    cpus 4
    memory "64GB"
    
    container "broadinstitute/gatk3:3.6-0"

    publishDir "${params.outdir}/somatic/"

    input:
        tuple val(mutect_meta), path(mutect_vcf), path(mutect_vcf_index)
        tuple val(strelka_meta), path(strelka_vcf), path(strelka_vcf_index)
        path(reference_fa)
        path(reference_index_dir)
    
    output:
        tuple val(merged_meta), path("${merged_meta.somatic_name}_variants.vcf.gz")


    script:
        merged_meta = [
            somatic_name: "${mutect_meta.somatic_name}",
            somatic_caller: 'mutect_strelka',
            tumor_metamap: "${mutect_meta.tumor_meta}",
            normal_metamap: "${mutect_meta.normal_meta}"
        ]

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
            -R $reference_fa \
            -genotypeMergeOptions PRIORITIZE \
            --rod_priority_list mutect,strelka \
            -V:mutect $mutect_vcf \
            -V:strelka $strelka_vcf \
            -o "${merged_meta.somatic_name}_variants.vcf.gz"

        """
    
}

process FILTER_VCF {

   cpus 4
   memory "32GB"
   
   container "biocontainers/bcftools:v1.9-1-deb_cv1"

   input:
       tuple val(somatic_meta), path(somatic_vcf)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_variants.vcf.gz"), path("${somatic_meta.somatic_name}_variants.vcf.gz.tbi"), emit: vcf

    script:
        """
        bcftools index $somatic_vcf
        bcftools view -f PASS -Oz -o "${somatic_meta.somatic_name}_variants.vcf.gz" $somatic_vcf
        bcftools index -t "${somatic_meta.somatic_name}_variants.vcf.gz"
        """
        


}


process ADD_VCF_GT_FIELD {

    cpus 4
    memory "32GB"
    
    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index)
         
    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${somatic_meta.somatic_caller}_gt.vcf.gz"), emit: gt_vcf

    script:
        """
        vcf-genotype-annotator $somatic_vcf \
            "${somatic_meta.tumor_metamap.sample_name}" \
            0/1 \
            -o "${somatic_meta.somatic_name}_${somatic_meta.somatic_caller}_gt.vcf.gz"

        """

}

process POSTPROCESS_STRELKA {

    cpus 4
    memory "32GB"
    
    container "biocontainers/bcftools:v1.9-1-deb_cv1"
    
    input:
        tuple val(somatic_meta),
              path(strelka_snvs), path(strelka_snvs_index), 
              path(strelka_indels), path(strelka_indels_index)

    output:
        tuple val(somatic_meta),
              path("${somatic_meta.somatic_name}_strelka.vcf.gz"), path("${somatic_meta.somatic_name}_strelka.vcf.gz.tbi"), 
              emit: strelka_vcf
    
    script:
        """
        bcftools concat \
            --allow-overlaps \
            --remove-duplicates \
            -Oz \
            -o "${somatic_meta.somatic_name}_strelka_merged.vcf.gz" \
            --threads $task.cpus \
            $strelka_snvs $strelka_indels

        bcftools view -s NORMAL,TUMOR "${somatic_meta.somatic_name}_strelka_merged.vcf.gz" \
            -Oz -o "${somatic_meta.somatic_name}_strelka_merged_order_samples.vcf.gz"

        echo ${somatic_meta.normal_metamap.sample_name} > new_names.txt
        echo ${somatic_meta.tumor_metamap.sample_name} >> new_names.txt
        
        bcftools reheader \
            --samples new_names.txt \
            --output "${somatic_meta.somatic_name}_strelka.vcf.gz" \
            "${somatic_meta.somatic_name}_strelka_merged_order_samples.vcf.gz"

        bcftools index -t "${somatic_meta.somatic_name}_strelka.vcf.gz"
        """
}


process STRELKA {

   container 'quay.io/wtsicgp/strelka2-manta'

   input:
        tuple val(tumor_meta), path(tumor_bam), path(tumor_bam_index)
        tuple val(normal_meta), val(normal_bam), path(normal_bam_index)
        val(somatic_name)
        path reference_fa
        path reference_index_dir

    output:

        tuple val(somatic_meta),
              path("./strelka/results/variants/somatic.snvs.vcf.gz"), path("./strelka/results/variants/somatic.snvs.vcf.gz.tbi"),
              path("./strelka/results/variants/somatic.indels.vcf.gz"), path("./strelka/results/variants/somatic.indels.vcf.gz.tbi"),
              emit: strelka_vcf

    script:
        somatic_meta = [
            somatic_name: somatic_name,
            somatic_caller: 'strelka',
            tumor_metamap: tumor_meta,
            normal_metamap: normal_meta
        ]
        """
        configManta.py \
            --normalBam $normal_bam \
            --tumorBam $tumor_bam \
            --referenceFasta "${reference_fa}" \
            --runDir ./manta/ \
            --exome

        ./manta/runWorkflow.py -j $task.cpus

        configureStrelkaSomaticWorkflow.py \
            --normalBam $normal_bam \
            --tumorBam $tumor_bam \
            --referenceFasta "${reference_fa}" \
            --indelCandidates ./manta/results/variants/candidateSmallIndels.vcf.gz \
            --runDir "./strelka" \
            --exome

        ./strelka/runWorkflow.py -m local -j $task.cpus
        
        """

}

process POSTPROCESS_MUTECT2_SCATTER {

    cpus 4
    memory "64GB"

    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple val(meta), path(vcfs), path(vcf_indices), path(f1r2s), path(stats)
        tuple val(tumor_meta), path(tumor_bam), path(tumor_bam_index), path(tumor_pileups)
        tuple val(normal_meta), path(normal_bam), path(normal_bam_index), path(normal_pileups)
        //tuple path(tumor_pileups), path(normal_pileups)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(meta), path("${meta.somatic_name}_mutect_postproc.vcf.gz"), emit: mutect_vcf
        //path("${meta.somatic_name}_mutect_postproc.vcf.gz.tbi"), emit: mutect_vcf

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
            -O "${meta.somatic_name}_merged.vcf"


        gatk LearnReadOrientationModel \
            $f1r2_as_input \
            -O "${meta.somatic_name}_orientmodel.tar.gz"

        gatk MergeMutectStats \
            $stat_as_input \
            -O "${meta.somatic_name}_merged.stats"

        gatk FilterMutectCalls \
            -R $reference_fa \
            -V "${meta.somatic_name}_merged.vcf" \
            --orientation-bias-artifact-priors "${meta.somatic_name}_orientmodel.tar.gz" \
            --contamination-table contamination.table \
            -stats "${meta.somatic_name}_merged.stats" \
            -O "${meta.somatic_name}_mutect_postproc.vcf.gz"

        """
}


process MUTECT2_SCATTER {

    cpus 2
    memory "24GB"

    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple val(tumor_meta), path(tumor_bam), path(tumor_bam_index), path(tumor_pileups)
        tuple val(normal_meta), path(normal_bam), path(normal_bam_index), path(normal_pileups)
        val somatic_name
        each path(interval_shard)
        path reference_fa
        path reference_index_dir
        path known_sites
        path known_sites_dir
        path pon
        path pon_index_dir

    output:
        tuple val(somatic_meta),
              path("${somatic_name}_${interval_shard}_mutect.vcf.gz"), path("${somatic_name}_${interval_shard}_mutect.vcf.gz.tbi"), 
              path("${somatic_name}_${interval_shard}_mutect_f1r2.tar.gz"), path("*.stats"), 
              emit: mutect_scatter_vcf

        //tuple path(tumor_pileups), path(normal_pileups), emit: tumor_normal_pileups


    script:

        somatic_meta = [
            somatic_name: somatic_name,
            somatic_caller: 'mutect',
            tumor_metamap: tumor_meta,
            normal_metamap: normal_meta
        ]

        """
        gatk Mutect2 \
            -R "${reference_fa}" \
            -I ${tumor_bam} \
            -I ${normal_bam} \
            -normal ${normal_meta.sample_name} \
            --germline-resource "${known_sites}" \
            --panel-of-normals "${pon}" \
            --f1r2-tar-gz "${somatic_name}_${interval_shard}_mutect_f1r2.tar.gz" \
            -L $interval_shard \
            -O "${somatic_name}_${interval_shard}_mutect.vcf.gz"

        """
}

process SPLIT_INTERVALS {

    cpus 4
    memory "16GB"

    tag "Split intervals for ${params.scatter_count} shards"
    
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
        -R $reference_fa \
        -L $intervals_file \
        --scatter-count $scatter_count \
        --interval-padding 100 \
        -O .
        
    """
}

process KALLISTO_QUANT {
    conda "bioconda::kallisto=0.51.1"

    publishDir "${params.outdir}/rnaseq/kallisto/", mode: "copy"

    input:
        tuple val(meta), path(read1), path(read2)
        path kallisto_index

    output:
        tuple val(meta), path("${meta.sample_name}_kallisto")

    script:
        """
        kallisto quant -i $kallisto_index -o ${meta.sample_name}_kallisto -t ${task.cpus} $read1 $read2 
        """
}


process KALLISTO_INDEX {

    conda "bioconda::kallisto=0.51.1"

    input:
        path(kallisto_reference)

    output:
        path("kallisto_index.idx")

    script:
        """
        kallisto index -i kallisto_index.idx $kallisto_reference

        """

}

process PREPROCESS_BAM {

    cpus 32

    memory "128GB"
    
    container "broadinstitute/gatk:4.6.1.0"

    publishDir "${params.outdir}/preprocess_bam/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), path(reads), path(reads_index)
        path reference_fa
        path reference_fa_index
        path known_sites
        path known_sites_index

    output:
        tuple val(meta), 
            path("${meta.sample_name}_${meta.molecule}_bqsr.bam"),
            path("${meta.sample_name}_${meta.molecule}_bqsr.bai"), 
        emit: preproc_bams
        
        tuple val(meta), 
            path("${meta.sample_name}_${meta.molecule}_bqsr.bam"),
            path("${meta.sample_name}_${meta.molecule}_bqsr.bai"),
            path("${meta.sample_name}_${meta.molecule}_pileups.table"), 
        emit: preproc_bams_pileups
        
        //path("${sample_id}/${sample_id}_recal_table.table"), emit: recal_table
        //path("${sample_id}/${sample_id}_dedup.bam"), emit: dedup_bam

    script:
        """
        gatk MarkDuplicatesSpark \
            -I $reads \
            -O "${meta.sample_name}_${meta.molecule}_dedup.bam" \
            --tmp-dir "\${PWD}"

        gatk BaseRecalibrator \
            -I "${meta.sample_name}_${meta.molecule}_dedup.bam" \
            -O "${meta.sample_name}_${meta.molecule}_recal_table.table" \
            -R $reference_fa \
            --known-sites $known_sites
        
        gatk ApplyBQSR \
            -R $reference_fa \
            -I "${meta.sample_name}_${meta.molecule}_dedup.bam" \
            --bqsr-recal-file "${meta.sample_name}_${meta.molecule}_recal_table.table" \
            -O "${meta.sample_name}_${meta.molecule}_bqsr.bam" \
            --create-output-bam-index

        gatk GetPileupSummaries \
            -I "${meta.sample_name}_${meta.molecule}_bqsr.bam" \
            -V "${known_sites}" \
            -L "${known_sites}" \
            -O "${meta.sample_name}_${meta.molecule}_pileups.table"
        """
}

process BWA_MAP {

    cpus 32

    memory "64GB"

    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.22.1"

    publishDir "${params.outdir}/alignment/bwa/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path reference_fa
        path reference_index_dir
        path(bwa_index)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_sorted.bam.bai"), emit: mapped_bam
    
    script:
    """
    NEW_RG="@RG\\tID:${meta.sample_name}\\tSM:${meta.sample_name}\\tLB:${meta.sample_name}\\tPL:ILLUMINA"

    bwa mem -M -t $task.cpus -R \$NEW_RG $reference_fa ${fastq1} ${fastq2} \
        | samtools sort --threads $task.cpus -o ${meta.sample_name}_${meta.molecule}_sorted.bam
    samtools index ${meta.sample_name}_${meta.molecule}_sorted.bam
    """
}


process CREATE_BWA_INDEX {
    
    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.22.1"


    input:
        path reference_fa
        path reference_fa_index

    output:
        path "*{.bwt,.sa,.pac,.amb,.ann}", emit:bwa_index

    script:
        """
        bwa index $reference_fa
        """




}

process STAR_FUSION {

    container "trinityctat/starfusion:1.15.0"

    publishDir "${params.outdir}/fusions/star_fusion/"

    input:
        tuple val(meta), path(chimeric_out)
        path ctat_resource_lib

    output:

        tuple val(meta), path("./${meta.sample_name}_starfusion/*.fusion_predictions.tsv"), emit: fusion_preds
        tuple val(meta), path("./{meta.sample_name}_starfusion/*.fusion_predictions.abridged.tsv"), emit: abridged_preds


    script:
        """

        STAR-Fusion --genome_lib_dir $ctat_resource_lib \
             -J $chimeric_out \
             --output_dir "./${meta.sample_name}_starfusion"


        """



}

process STAR_INDEX_BAM {

    container "biocontainers/samtools:v1.9-4-deb_cv1"
    
    publishDir "${params.outdir}/alignment/star/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), path(bam)
        tuple val(meta), path(final_log)
        tuple val(meta), path(sj_out)
        tuple val(meta), path(chimeric_out)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_STAR_sorted.bam.bai"), emit: star_bam
        path final_log
        path sj_out
        path chimeric_out

    script:
        """
        samtools sort --threads $task.cpus  $bam -o "${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"
        samtools index -@ $task.cpus  "${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"

        """


}

process STAR_ALIGN {

    cpus 64

    memory "256GB"

    container "alexdobin/star:2.7.10a_alpha_220506"

    //publishDir "${params.outdir}/alignment/star_raw/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path(star_index_dir)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Aligned.out.bam"), emit: star_bam 
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Log.final.out"), emit:final_log
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_SJ.out.tab"), emit: sj_out
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Chimeric.out.junction"), emit: chimeric_out
    script:
        """
        STAR \
            --runThreadN $task.cpus \
            --genomeDir $star_index_dir \
            --readFilesIn $fastq1 $fastq2 \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped None \
            --twopassMode Basic \
            --outSAMstrandField intronMotif \
            --outSAMunmapped Within \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 8 \
            --chimOutJunctionFormat 1 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --alignIntronMax 100000 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outSAMattrRGline ID:"${meta.sample_name}" SM:"${meta.sample_name}" \
            --chimMultimapScoreRange 3 \
            --chimScoreJunctionNonGTAG -4 \
            --chimMultimapNmax 20 \
            --chimNonchimScoreDropMin 10 \
            --peOverlapNbasesMin 12 \
            --peOverlapMMp 0.1 \
            --alignInsertionFlush Right \
            --alignSplicedMateMapLminOverLmate 0 \
            --alignSplicedMateMapLmin 30 \
            --outFileNamePrefix ./${meta.sample_name}_${meta.molecule}_

        """

}


process CREATE_STAR_INDEX {

    container "alexdobin/star:2.7.10a_alpha_220506"

    input:
        path(reference_fa)
        path(reference_index_files)
        path(gtf)

    output:
        path("./STARGenomeDir"), emit: star_index

    script:
        """
        gzip -d -c $gtf > gencode.gtf

        STAR \
            --runThreadN $task.cpus \
            --runMode genomeGenerate \
            --genomeDir ./STARGenomeDir \
            --genomeFastaFiles $reference_fa \
            --sjdbGTFfile gencode.gtf

        """


}



process FASTP {
    
    cpus 32
    memory "128GB"

    conda "bioconda::fastp=1.0.1"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_R1_fastp.fastq.gz"), path("${meta.sample_name}_${meta.molecule}_R2_fastp.fastq.gz"), emit: fastqs
        path("${meta.sample_name}_*{html,json}*"), emit: reports
    
    tag "FastP on ${meta.sample_name} w/ ${meta.molecule}"

    publishDir "${params.outdir}/fastp_qc/${meta.sample_name}_${meta.molecule}_fastp/", mode: 'copy'

    script:
        """
        fastp --thread $task.cpus \
              -i ${reads[0]} \
              -I ${reads[1]} \
              -o "${meta.sample_name}_${meta.molecule}_R1_fastp.fastq.gz" \
              -O "${meta.sample_name}_${meta.molecule}_R2_fastp.fastq.gz" \
              -R "${meta.sample_name}_${meta.molecule}_fastp_report" \
              -h "${meta.sample_name}_${meta.molecule}_fastp_report.html" \
              -j "${meta.sample_name}_${meta.molecule}_fastp_report.json" \

        """
}

process OPTITYPE_HLA_CALLS {
    
    cpus 1
    memory "16GB"

    conda "python=3.10 pandas=2.1"

    publishDir "${params.outdir}/HLA/"

    input:
        tuple val(meta), path(optitype_result_tsv)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.sample_type}_${meta.molecule}_hla_pvacinput.csv")

    script:
        """
        #!/usr/bin/env python3

        import pandas as pd

        df = pd.read_csv("$optitype_result_tsv", sep = "\t")
        
        alleles = []
        for allele in ["A1", "A2", "B1", "B2", "C1", "C2"]:
            alleles.append("HLA-" + df[allele].iloc[0])
        allele_csv = ",".join(alleles)

        with open("${meta.sample_name}_${meta.sample_type}_${meta.molecule}_hla_pvacinput.csv", "w") as f:
            print(allele_csv, file=f)

        """

}
process POSTPROCESS_OPTITYPE {

    cpus 1
    memory "32GB"

    publishDir "${params.outdir}/HLA/optitype/${meta.sample_name}_optitype", mode: "copy"

    input:
        tuple val(meta), path(optitype_output_dir)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.sample_type}_${meta.molecule}_optitype.tsv"), emit: result_tsv
        tuple val(meta), path("${meta.sample_name}_${meta.sample_type}_${meta.molecule}_optitype_coverage.pdf"), emit: coverage_plot

    script:
        """
        RESULT_FILE=\$(find ${optitype_output_dir}/* -name "*_result.tsv" | head -n 1)    
        mv "\$RESULT_FILE" "${meta.sample_name}_${meta.sample_type}_${meta.molecule}_optitype.tsv"
        
        COVERAGE_FILE=\$(find ${optitype_output_dir}/* -name "*_coverage_plot.pdf" | head -n 1)    
        mv "\$COVERAGE_FILE" "${meta.sample_name}_${meta.sample_type}_${meta.molecule}_optitype_coverage.pdf"

        """
}


process OPTITYPE {
    
    container "fred2/optitype:release-v1.3.1"

    cpus 16

    memory "128GB"

    input:
        tuple val(meta), path(fastq1), path(fastq2)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_optitype")

    script:

        def molecule_flag = meta.molecule.toLowerCase()
        """
        cat << EOF > OptiType.ini
        [mapping]
        razers3=/usr/local/bin/razers3
        threads=${task.cpus}
        [ilp]
        solver=cbc
        threads=${task.cpus}
        [behavior]
        deletebam=true
        unpaired_weight=0
        use_discordant=false
        EOF
        
        python /usr/local/bin/OptiType/OptiTypePipeline.py \
            -i $fastq1 $fastq2 \
            --$molecule_flag \
            -c OptiType.ini \
            --outdir "${meta.sample_name}_${meta.molecule}_optitype"
        """
}


process HLAHD {
    
    cpus 32

    memory "32GB"


    container "griffithlab/hlahd:1.0"
    
    publishDir "${params.outdir}/HLA/hlahd/${meta.sample_name}_hlahd", mode: "copy"


    input:
        tuple val(meta), path(fastq1), path(fastq2)
    
    output:
        tuple val(meta), path("./${meta.sample_name}/result/${meta.sample_name}_final.result.txt"), emit: hla_calls
        path("./${meta.sample_name}/result/*")

    script:
        """
        /opt/hlahd/bin/hlahd.1.6.1.sh \
            -f /opt/hlahd/freq_data \
            -t $task.cpus \
            $fastq1 $fastq2 \
            /opt/hlahd/HLA_gene.split.txt \
            /opt/hlahd/dictionary \
            "${meta.sample_name}" \
            ./
        """


}

process PVACSEQ {

    cpus 32
    memory "256GB"

    container "griffithlab/pvactools:6.0.1"

    publishDir "${params.outdir}/pvactools/"

    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index)
        tuple val(phased_meta), path(phased_vcf), path(phased_vcf_index)
        tuple val(tumor_meta), path(tumor_reads), path(tumor_reads_index)
        tuple val(normal_meta), path(normal_reads), path(normal_reads_index)
        path(hla_pvac_input)
    output:
        path("${somatic_meta.somatic_name}_pvacseq"), emit: pvacseq_dir

    script:
        """
        pvacseq run \
            $somatic_vcf \
            ${tumor_meta.sample_name} \
            \$(head $hla_pvac_input -n 1) \
            all \
            "${somatic_meta.somatic_name}_pvacseq" \
            -e1 8,9,10,11 \
            -e2 12,13,14,15,16,17,18 \
            --normal-sample-name ${normal_meta.sample_name} \
            --phased-proximal-variants-vcf $phased_vcf \
            --iedb-install-directory /opt/iedb \
            --pass-only \
            -t $task.cpus
        """
}


workflow {
    
    // INPUT PARSING
    // READ IN SAMPLE DATA FROM SAMPLESHEET


    samplemap_inputs = Channel.fromPath(params.sample_sheet)
        | splitCsv( header: true )
            | map { row ->
                meta = [
                    sample_name: row.sample_name,
                    sample_type: row.sample_type,
                    molecule: row.molecule,
                    sequencing_type: row.sequencing_type
                ]
                
                reads = [
                    file(row.fastqr1, checkIfExists: true),
                    file(row.fastqr2, checkIfExists: true)
                ]
            [meta, reads]
        }
    

    samples_by_type = samplemap_inputs.branch { meta, reads ->
            dna_samples: meta.sample_type == "DNA"
            rna_samples: meta.sample_type == "RNA"
    }


            
    /* 
    // Reads in the GATK resource bucket for reference genome hg38
    
    all_reference_files = Channel
                .fromPath("${params.reference_index_dir}*", checkIfExists: true)
                .collect()

    all_reference_files
        | flatten
        | branch {
            fasta: it.name.endsWith('.fa') || it.name.endsWith('.fasta')
            index_files: true
        }
    | set { ref_files }

    reference_index_files = ref_files.index_files.collect()
    reference_fa = ref_files.fasta.collect()
    reference_dict = Channel.fromPath("${params.reference_index_dir}*.dict")
    */


    // PULL REFERENCE FASTA AND REFERENCE FASTA SUPPLEMENTAL FILES

    reference_fa = Channel.fromPath(params.reference_fa).first()
    reference_index_files = Channel.fromPath(params.reference_index_dir).collect()
    reference_dict = Channel.fromPath(params.reference_dict).first()
    common_germline = Channel.fromPath(params.common_germline).first()
    common_germline_index = Channel.fromPath(params.common_germline_index).first()
    known_sites = Channel.fromPath(params.known_sites).first()
    known_sites_index = Channel.fromPath(params.known_sites_index).first()
    pon = Channel.fromPath(params.pon).first()
    pon_index = Channel.fromPath(params.pon_index).first()
    hapmap = Channel.fromPath(params.hapmap).first()
    hapmap_index = Channel.fromPath(params.hapmap_index).first()
    mills = Channel.fromPath(params.mills).first()
    mills_index = Channel.fromPath(params.mills_index).first()
    intervals_file = Channel.fromPath(params.intervals_file).first()
    kallisto_reference = Channel.fromPath(params.kallisto_reference).first()
    

    // RUN FASTP QC ON ALL SAMPLES

    fastp = FASTP(samplemap_inputs)
    fastp_by_molec = fastp.fastqs.branch{meta, fastq1, fastq2 ->
                                        dna: meta.molecule == "DNA"
                                        rna: meta.molecule == "RNA"
                                        }



    // HLA TYPING: RUN OPTITYPE AND HLA-HD

    hla_optitype = OPTITYPE(fastp_by_molec.dna)
    hla_hlahd = HLAHD(fastp_by_molec.dna)
    hla_optitype_postprocess = POSTPROCESS_OPTITYPE(hla_optitype)
    if (params.use_clinical_hla_calls) {
        hla_calls = params.clinical_hla_calls
        //hla_calls = PARSE_CLINICAL_CALLS() // NOT IMPLEMENTED
    } else { 
        hla_calls = OPTITYPE_HLA_CALLS(hla_optitype_postprocess.result_tsv)
        hla_calls = hla_calls.branch { meta, hla_call ->
                      normal_dna: meta.sample_type == "Normal" && meta.molecule == "DNA"
                        return hla_call
                      other: true
                        return hla_call
                    }
        hla_calls = hla_calls.normal_dna
    }

    
    
    // ALIGNMENT OF DNA SEQUENCING USING BWA


    bwa_index = CREATE_BWA_INDEX(reference_fa, reference_index_files)

    bwa_mapped = BWA_MAP(fastp_by_molec.dna, reference_fa, reference_index_files, bwa_index.bwa_index)
    bwa_mapped.mapped_bam.branch { meta, bam, bam_index ->
                                 mapped_dna: meta.molecule == "DNA"
                                 mapped_rna: meta.molecule == "RNA"
                                 other: true
                                 }
                             | set { mapped_dna_rna_branch }


    // ALIGNMENT OF RNA SEQUENCING USING STAR
    // STAR alignment process includes chimeric reads for directly going from BAM -> Star-Fusion
    gencode_gtf = Channel.fromPath(file(params.gencode_gtf))
    star_index = CREATE_STAR_INDEX(reference_fa, reference_index_files, gencode_gtf).star_index
    star_rna_align = STAR_ALIGN(fastp_by_molec.rna, star_index.collect())
    star_rna = STAR_INDEX_BAM(star_rna_align)
    

    // BAM PREPROCESSING OF DNA: GATK BEST PRACTICES

    preprc = PREPROCESS_BAM(mapped_dna_rna_branch.mapped_dna,
                            reference_fa, reference_index_files, 
                            common_germline, common_germline_index)

    
    preproc_bams_type_branched = preprc.preproc_bams
                                    | branch {meta, bam, bai ->
                                        normal: meta.sample_type == "Normal"
                                        tumor: meta.sample_type == "Tumor"
    }
    preproc_bams_pileups_type_branched = preprc.preproc_bams_pileups
                                            | branch {meta, bam, bai, pileups ->
                                                normal: meta.sample_type == "Normal"
                                                tumor: meta.sample_type == "Tumor"
                                            }
    

    //deepsomatic = DEEPSOMATIC(preproc_bams_type_branched.tumor, preproc_bams_type_branched.normal,
    //                          params.somatic_name, reference_fa, reference_index_files)
    


    // RNA: TRANSCRIPT ABUNDANCE ESTIMATION
    

    kallisto_index = KALLISTO_INDEX(kallisto_reference)
    kallisto = KALLISTO_QUANT(fastp_by_molec.rna,
                              kallisto_index)



    // INTERVAL CREATION FOR SOMATIC CALLERS

    intervals = SPLIT_INTERVALS(reference_fa, reference_index_files,
                                intervals_file, params.scatter_count)



    // GERMLINE VARIANT CALLING

    // Run HaplotypeCaller with scatter/gather approach and preprocess with CNNScoreVariants and FilterVariantTranches
    // For later use in phasing the somatic VCF with proximal variants
    germline = HAPLOTYPE_CALLER_SCATTER(preproc_bams_type_branched.normal,
                intervals.flatten(),
                reference_fa, reference_index_files)
        | groupTuple
    germline_postprocess = POSTPROCESS_HAPLOTYPE_SCATTER(germline, 
                                                        reference_fa,
                                                        reference_index_files,
                                                        hapmap,
                                                        hapmap_index,
                                                        mills,
                                                        mills_index)
    // Use vt decompose on germline calls
    vt_germline = VT_DECOMPOSE_GERMLINE(germline_postprocess)


    //
    // SOMATIC VARIANT CALLING
    
    // MUTECT2
    // Run Mutect2 with scatter/gather approach
    
    mutect_scattered = MUTECT2_SCATTER(preproc_bams_pileups_type_branched.tumor, preproc_bams_pileups_type_branched.normal, 
                    params.somatic_name,
                    intervals.flatten(),
                    reference_fa, reference_index_files,
                    known_sites, known_sites_index,
                    pon, pon_index)

    mutect_gathered = mutect_scattered.mutect_scatter_vcf
                        | groupTuple
                       
    // Postprocess Mutect2 Output using pileups, stats, etc.
    mutect_postprocess = POSTPROCESS_MUTECT2_SCATTER(mutect_gathered,
                                                     preproc_bams_pileups_type_branched.tumor, 
                                                     preproc_bams_pileups_type_branched.normal,
                                                     reference_fa,
                                                     reference_index_files)
    
    
    // STRELKA
    strelka = STRELKA(preproc_bams_type_branched.tumor, preproc_bams_type_branched.normal, 
            params.somatic_name,
            reference_fa, reference_index_files)

    strelka_postprocess = POSTPROCESS_STRELKA(strelka)
    
    strelka = ADD_VCF_GT_FIELD(strelka_postprocess.strelka_vcf)

    somatic_vcfs = mutect_postprocess.concat(strelka)
    
    somatic_filtered_vcfs = FILTER_VCF(somatic_vcfs)
    

    // POST-PROCESS SOMATIC VCFS

    // Variant Decomposition & Normalization
    somatic_vcfs_vt = VT_SOMATIC_POSTPROCESS(somatic_filtered_vcfs, reference_fa, reference_index_files)

    somatic_vcf_branched = somatic_vcfs_vt.branch { vcf ->
        mutect: vcf[0].somatic_caller == "mutect"
        strelka: vcf[0].somatic_caller == "strelka"
    }

    merged_vcf = MERGE_SOMATIC_VCFS(somatic_vcf_branched.mutect,
                                      somatic_vcf_branched.strelka,
                                      reference_fa,
                                      reference_index_files)


    // PVACtools VCF PREPARATION
    
    // Annotate VCF with VEP
    vep = VEP_ANNOTATE(merged_vcf,
                       reference_fa,
                       params.vep_cache,
                       params.vep_plugins)
    
    
    // Add Read Coverage to VCF with Bamreadcount
    bamreadcount = BAMREADCOUNT(vep.vep_vcf,
                                reference_fa,
                                star_rna.star_bam)
    vcf_annotated_coverage = ANNOTATE_VCF_COVERAGE(bamreadcount)
    
    // Add transcript abundance estimation from kallisto
    vcf_annotated_expression = ANNOTATE_VCF_EXPRESSION(vcf_annotated_coverage, kallisto)

    // Index final somatic vcf for input to PVACseq
    vcf_final = INDEX_FINAL_VCF(vcf_annotated_expression)
    

    

    // PERFORM VCF PHASING USING GERMLINE CALLS

    // Create Tumor-Only VCF From Final Somatic VCF (vcf_final)
    vcf_phase_select_variants = PHASE_VCF_SELECT_VARIANTS(preproc_bams_type_branched.tumor,
                           vcf_final,
                           reference_fa,
                           reference_index_files)

    vcf_phase_rename_samples = PHASE_VCF_RENAME(vt_germline,
                                                preproc_bams_type_branched.tumor)

    // Combine Tumor-Only VCF with germline variants
    vcf_phase_combine = PHASE_VCF_COMBINE_VARIANTS(vcf_phase_select_variants,
                                                   vcf_phase_rename_samples,
                                                   reference_fa,
                                                   reference_index_files)

    // Sort combined VCF
    vcf_phase_sort = PHASE_VCF_SORT_VCF(vcf_phase_combine, reference_dict)
    

    // Call GATK ReadBackedPhasing (required gatk 3.6.0) to phase VCF
    vcf_phase_rbphase = PHASE_VCF_RBPHASING(vcf_phase_sort.vcf, 
                                            reference_fa, reference_index_files,
                                            preproc_bams_type_branched.tumor)

    /// VEP Annotation phased vcf
    vcf_phase_vep = PHASE_VCF_VEP(vcf_phase_rbphase,
                                    reference_fa, params.vep_cache, params.vep_plugins)        
    
    // Zip and index phased vcf
    vcf_phased = PHASE_VCF_INDEX(vcf_phase_vep)

    
    PVACSEQ = PVACSEQ(vcf_final, 
                      vcf_phased,
                      preproc_bams_type_branched.tumor,
                      preproc_bams_type_branched.normal,
                      hla_calls)
}
