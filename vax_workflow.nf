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

process VT_POSTPROCESS_GERMLINE {

    /*

        Germline variant post-processing similar to VT_SOMATIC_POSTPROCESS using
        variant normalization, decomposition, sort, and duplicate removal and indexing.

    */

   
    cpus 2
    memory "8GB"

    container "staphb/bcftools:1.23"


    input:
        tuple val(meta), path(germline_vcf), path(germline_vcf_index)
        path reference_fa
        path reference_index_files

    output:
        tuple val(meta), path("${meta.sample_name}_germline_normalized.vcf.gz"), path("${meta.sample_name}_germline_normalized.vcf.gz.tbi"), emit: germline_vcf
    script:
        """
        bcftools view -f PASS $germline_vcf -Oz -o filtered_germline.vcf.gz
        bcftools norm -m -any -f $reference_fa filtered_germline.vcf.gz -Oz -o norm_vcf.vcf.gz
        bcftools sort norm_vcf.vcf.gz -Oz -o norm_sort.vcf.gz
        bcftools norm -d exact norm_sort.vcf.gz -Oz -o "${meta.sample_name}_germline_normalized.vcf.gz"
        bcftools index -t "${meta.sample_name}_germline_normalized.vcf.gz"
     """
}

process VEP_ANNOTATE_GERMLINE {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    cpus 8
    memory "16GB"
    cache "lenient"

    
    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(meta), path(germline_vcf), path(germline_vcf_index)
        path reference_fa
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("germline_vep.vcf"), emit: vep_vcf

    script:
        """
        vep \
            --input_file $germline_vcf  \
            --output_file germline_vep.vcf \
            --everything \
            --format vcf --vcf --symbol --terms SO --mane_select --canonical --tsl --biotype --hgvs \
            --fasta $reference_fa  \
            --offline --cache \
            --plugin Frameshift --plugin Wildtype \
            --pick \
            --dir_plugins $vep_plugins \
            --dir_cache $vep_cache \
            --transcript_version
        
        """
}

process INDEX_FINAL_VCF_GERMLINE {

    /*

    tabix index a vcf file

    */

    cpus 2
    memory "8GB"

    container "staphb/bcftools:1.23"


    publishDir "${params.outdir}/${meta.somatic_sample}/germline/", mode: "copy"

    input:
        tuple val(meta), path(germline_vcf)
        

    output:
        tuple val(meta), path("${meta.sample_name}_germline_filt_vep.vcf.gz"), path("${meta.sample_name}_germline_filt_vep.vcf.gz.tbi")

    script:
        """
        bcftools view $germline_vcf -Oz -o "${meta.sample_name}_germline_filt_vep.vcf.gz"
        bcftools index -t "${meta.sample_name}_germline_filt_vep.vcf.gz"
        """

}




process PHASE_VCF_INDEX {

    /*

        Index the final phased vcf.

    */

    cpus 1
    memory "16GB"
    
    conda "bioconda::tabix=0.2.6"

    publishDir "${params.outdir}/${meta.somatic_name}/variants/phased_variants/", mode: "copy"

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

    /*

        Use VEP to annotated the phased vcf.

    */
    cpus 8
    memory "32GB"


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
            --plugin Frameshift --plugin Wildtype --plugin Downstream \
            --pick \
            --fork ${task.cpus} \
            --dir_plugins $vep_plugins \
            --transcript_version
    
        """

}

process PHASE_VCF_RBPHASING {

    /*

    Use deprecated ReadBackedPhasing from GATK 3.6.0 to phase variants 
    in the somatic + germline combined and sorted vcf

    */
    
    cpus 4
    memory "32GB"


    container "broadinstitute/gatk3:3.6-0" 

    input:
        tuple val(somatic_meta), path(combined_sorted_vcf), val(tumor_meta), path(tumor_reads), path(tumor_reads_index)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(somatic_meta), path("phased.vcf")

    script:

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T ReadBackedPhasing \
                -R $reference_fa \
                -I $tumor_reads \
                --variant $combined_sorted_vcf \
                -L $combined_sorted_vcf \
                -o phased.vcf

        """
}

process PHASE_VCF_RENAME {

    /*

    Part of creating a phase vcf file.
    Germline sample name will be the name of the NORMAL sample, but to combine variants with the tumor sample, the names must match.
    Here the sample name in the germline vcf is renamed to the name of the tumor sample.

    */

    cpus 2
    memory "16GB"

    container "biocontainers/bcftools:v1.9-1-deb_cv1"

    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index),
            path(tumor_only_vcf), path(tumor_only_vcf_index)

    output:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path("vt_germline_rename.vcf.gz" ), path("vt_germline_rename.vcf.gz.tbi"),
            path(tumor_only_vcf), path(tumor_only_vcf_index)

    script:
        """
        echo ${somatic_meta.tumor_metamap.sample_name} > new_names.txt
        
        bcftools reheader \
            --samples new_names.txt \
            --output vt_germline_rename.vcf.gz \
            $germline_vcf 

        bcftools index -t vt_germline_rename.vcf.gz
        """


}

process PHASE_VCF_SORT_VCF {

    /*

    Part of creating a phased vcf.
    Sort the somatic + germline combined vcf.

    */

    cpus 2
    memory "32GB"

    container 'broadinstitute/picard:3.4.0'
    
    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index),
            path(combined_vcf)
        path(reference_dict)

    output:
        tuple val(somatic_meta), path("combined.sorted.vcf"), emit: vcf

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

    /*

    Part of creating a phased germline vcf.
    Combines the variants from the tumor-only vcf and the germline_vcf (which has been renamed).


    */

    cpus 2
    memory "32GB"

    container "broadinstitute/gatk3:3.6-0"

    input:

        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index),
            path(tumor_only_vcf), path(tumor_only_vcf_index)

        path(reference_fa)
        path(reference_index_files)



    output:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index),
            path("combined_somatic_plus_germline.vcf")
         
    script:
        """
        
        java -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
                -R $reference_fa \
                --variant $germline_vcf \
                --variant $tumor_only_vcf \
                -o combined_somatic_plus_germline.vcf \
                --assumeIdenticalSamples

        """

}


process PHASE_VCF_SELECT_VARIANTS {

    /*

    Part of creating a phased germline vcf
    Takes a somatic vcf and extracts just the tumor sample.

    */

    cpus 2
    memory "32GB"
    cache "lenient"


    container "broadinstitute/gatk:4.6.1.0"
    
    input:

        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index),
            path("tumor_only.vcf.gz"), path("tumor_only.vcf.gz.tbi")

    script:
        """
        gatk SelectVariants \
            -V $somatic_vcf \
            -R "${reference_fa}" \
            --sample-name ${somatic_meta.tumor_metamap.sample_name} \
            -O tumor_only.vcf.gz

        gatk IndexFeatureFile \
            -I tumor_only.vcf.gz
        
        """

}

process POSTPROCESS_HAPLOTYPE_SCATTER {

    /*

    Gathers scattered HaplotypeCaller intervals and filters them using FilterVariantTranches

    */

    cpus 4
    memory "32GB"
    cache "lenient"

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

    /*

    Use HaplotypeCaller on a single scattered interval. Post processes with CNNScoreVariants.

    */
    
    cpus 4
    memory "24GB"
    cache "lenient"

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

    /*

    Use the abundance estimates from kallist to annotate transcript expression in a vcf file.

    */
    cpus 2
    memory "16GB"

    container "griffithlab/vatools:5.2.0"

    publishDir "${params.outdir}/${vcf_meta.somatic_name}/coverage/", mode: "copy"

    input:
        tuple val(vcf_meta), path(vcf), val(kallisto_meta), path(kallisto_dir)


    output:
        tuple val(vcf_meta), path("${vcf_meta.somatic_name}_cov_expr_annotated.vcf")

    script:

        """
        vcf-expression-annotator \
            $vcf \
            -s ${kallisto_meta.sample_name} \
            "${kallisto_dir}/abundance.tsv" \
            kallisto transcript \
            -o "${vcf_meta.somatic_name}_cov_expr_annotated.vcf"

        """

}

process ANNOTATE_VCF_COVERAGE {

    /*

        Use the output of bamreadcount to annotate coverage given a particular sample. 
    */

    cpus 4
    memory "32GB"
    
    container "griffithlab/vatools:5.2.0"

    publishDir "${params.outdir}/${somatic_meta.somatic_name}/coverage/", mode: "copy"

    input:
        tuple val(somatic_meta), path(vcf), path(brc_indels), path(brc_snvs)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_annotated.vcf")

    script:
        """

        echo "${somatic_meta.tumor_metamap.sample_name}"
        echo "${somatic_meta.tumor_metamap.molecule}"
        vcf-readcount-annotator \
            $vcf \
            $brc_snvs \
            RNA \
            -s ${somatic_meta.tumor_metamap.sample_name} \
            -t snv \
            -o "${somatic_meta.somatic_name}_snv_annotated.vcf"

        vcf-readcount-annotator \
            "${somatic_meta.somatic_name}_snv_annotated.vcf" \
            $brc_indels \
            RNA \
            -s ${somatic_meta.tumor_metamap.sample_name} \
            -t indel \
            -o ${somatic_meta.somatic_name}_annotated.vcf

        """

}

process BAMREADCOUNT {

    /*

    Run bamreadcount to add coverage information to a vcf using a BAM file. The sample names in the VCF
    must match the read groups and sample names in the BAM.

    */

    cpus 4
    memory "32GB"
    cache "lenient"

    container "mgibio/bam_readcount_helper-cwl:1.2.1"

    publishDir "${params.outdir}/${vcf_meta.somatic_name}/coverage/", mode: "copy"

    input:
        tuple val(vcf_meta), path(vcf), val(star_meta), path(star_bam), path(star_bai)
        path reference_fa

    output:
        tuple val(vcf_meta), path(vcf),
            path("${star_meta.sample_name}_bamrc_helper/${star_meta.sample_name}_bam_readcount_indel.tsv"), 
            path("${star_meta.sample_name}_bamrc_helper/${star_meta.sample_name}_bam_readcount_snv.tsv"), emit: brc_files
    script:
        
        """
        mkdir ${star_meta.sample_name}_bamrc_helper
        bam_readcount_helper.py \
            $vcf \
            ${star_meta.sample_name} \
            $reference_fa \
            $star_bam \
            NOPREFIX \
            ${star_meta.sample_name}_bamrc_helper
        """

}

process VT_SOMATIC_POSTPROCESS {

    /*

        Postprocess somatic variants from a vcf. This is used to normalize variants and decompose biallelics
        into different entries. The vcf is then sorted and duplicate entries removed. The final vcf is then indexed.

    */

    cpus 4
    memory "16GB"
    cache "lenient"
    
    container "staphb/bcftools:1.23"

    publishDir "${params.outdir}/${meta.somatic_name}/somatic/${meta.somatic_caller}/", mode:"copy"
    
    input:
        tuple val(meta), path(somatic_vcf), path(somatic_vcf_index)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(meta), 
            path("${meta.somatic_name}_${meta.somatic_caller}_normalized_variants.vcf.gz"), 
            path("${meta.somatic_name}_${meta.somatic_caller}_normalized_variants.vcf.gz.tbi"), emit: vt_vcf

    script:

        """
        bcftools norm -m -any -f $reference_fa $somatic_vcf -Oz -o norm_vcf.vcf.gz
        bcftools sort norm_vcf.vcf.gz -Oz -o norm_sort.vcf.gz
        bcftools norm -d exact norm_sort.vcf.gz -Oz -o "${meta.somatic_name}_${meta.somatic_caller}_normalized_variants.vcf.gz"
        bcftools index -t "${meta.somatic_name}_${meta.somatic_caller}_normalized_variants.vcf.gz"
        """
}


process VEP_ANNOTATE {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    cpus 8
    memory "16GB"
    cache "lenient"

    
    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(meta), path(somatic_vcf)
        path reference_fa
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("somatic_vep.vcf"), emit: vep_vcf

    script:
        """
        vep \
            --input_file $somatic_vcf  \
            --output_file somatic_vep.vcf \
            --everything \
            --format vcf --vcf --symbol --terms SO --mane_select --canonical --tsl --biotype --hgvs \
            --fasta $reference_fa  \
            --offline --cache \
            --plugin Frameshift --plugin Wildtype \
            --pick \
            --dir_plugins $vep_plugins \
            --dir_cache $vep_cache \
            --transcript_version
        
        """
}

process VEP_FILTER {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    cpus 2
    memory "16GB"
    cache "lenient"

    
    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(meta), path(vep_vcf)
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("${meta.somatic_name}_vep.vcf"), emit: vep_vcf

    script:
        """

        filter_vep -i $vep_vcf \
            -o "${meta.somatic_name}_vep.vcf" \
            --format vcf \
            --filter "gnomADe_AF < 0.001 or not gnomADe_AF"

        """
}




process INDEX_FINAL_VCF {
    /*

    tabix index a vcf file

    */

    cpus 2
    memory "8GB"

    conda "bioconda::tabix=0.2.6"

    publishDir "${params.outdir}/${meta.somatic_name}/variants", mode: "copy"

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

    /*

    Use the deprecated CombineVariants from GATK 3.6.0 to combine vcf files. 
    
    Currently just from mutect + strelka, but I should probably change this to
    take in a dynamic number of vcf files.

    */

    cpus 4
    memory "32GB"
    cache "lenient"
    
    container "broadinstitute/gatk3:3.6-0"

    publishDir "${params.outdir}/${merged_meta.somatic_name}/somatic/", mode: "copy"

    input:
        tuple val(somatic_name), val(call_sets) 
        path(reference_fa)
        path(reference_index_dir)
    
    output:
        tuple val(merged_meta), path("${merged_meta.somatic_name}_variants.vcf.gz")


    script:

        def call1 = call_sets[0]
        def call2 = call_sets[1]

        def vcf1_meta = call1[0]
        def vcf1  = call1[1]
        def tbi1  = call1[2]

        def vcf2_meta = call2[0]
        def vcf2  = call2[1]
        def tbi2  = call2[2]

        merged_meta = [
            somatic_name: vcf1_meta.somatic_name,
            somatic_caller: 'mutect_strelka',
            tumor_metamap: vcf1_meta.tumor_metamap,
            normal_metamap: vcf1_meta.normal_metamap
        ]

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
            -R $reference_fa \
            -genotypeMergeOptions PRIORITIZE \
            --rod_priority_list mutect,strelka \
            -V:${vcf1_meta.somatic_caller} $vcf1 \
            -V:${vcf2_meta.somatic_caller} $vcf2 \
            -o "${merged_meta.somatic_name}_variants.vcf.gz"

        """
    
}

process FILTER_VCF {

    /*

    Filter VCF, retaining only PASS variants.

    */
    
    cpus 4
    memory "16GB"
    
    container "biocontainers/bcftools:v1.9-1-deb_cv1"
    
    publishDir "${params.outdir}/${somatic_meta.somatic_name}/somatic/${somatic_meta.somatic_caller}/", mode:"copy"

    
    input:
        tuple val(somatic_meta), path(somatic_vcf)
        
    output:
        tuple val(somatic_meta), 
            path("${somatic_meta.somatic_name}_variants.vcf.gz"), 
            path("${somatic_meta.somatic_name}_variants.vcf.gz.tbi"), emit: filter_vcf
        
        tuple val(somatic_meta), 
            path("${somatic_meta.somatic_name}_variants_unfiltered.vcf.gz"), 
            path("${somatic_meta.somatic_name}_variants_unfiltered.vcf.gz.tbi"), emit: unfiltered_vcf
        
    script:
        """
        cp $somatic_vcf "${somatic_meta.somatic_name}_variants_unfiltered.vcf.gz"
        bcftools index -t "${somatic_meta.somatic_name}_variants_unfiltered.vcf.gz"
        bcftools view -f PASS -Oz -o "${somatic_meta.somatic_name}_variants.vcf.gz" $somatic_vcf
        bcftools index -t "${somatic_meta.somatic_name}_variants.vcf.gz"
        """
        


}


process ADD_VCF_GT_FIELD {

    /*

    Adds the GT fields to a vcf lacking it. This is currently used for Strelka/Manta vcfs and the field is populated by 0/1 by default.

    */

    cpus 2
    memory "16GB"
    
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

    /*

    Combine the SNVs and Indels from a Strelka and Manta run and rename the default TUMOR and NORMAL
    samples to the specified sample names.

    */

    cpus 4
    memory "16GB"
    
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

    /*
    
    Run Strelka and Manta SNV and Indel callers.


    */

   cpus 16
   memory "32GB"
   cache "lenient"

   container 'quay.io/wtsicgp/strelka2-manta'

   input:
        tuple val(somatic_name), val(tumor_meta), path(tumor_bam), path(tumor_bam_index),
                                 val(normal_meta), path(normal_bam), path(normal_bam_index)
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


               
    /*
    
    Workflow gathers all the shards of scattered mutect2 run and this process runs all
    post-processing steps.

    */

    cpus 4
    memory "32GB"
    cache "lenient"

    container "broadinstitute/gatk:4.6.1.0"
    
    // publishDir "${params.outdir}/${meta.somatic_name}/variants/${meta.somatic_name}_mutect_postproc.vcf.gz", mode: "copy"


    input:
        tuple val(somatic_name), val(meta), val(vcfs), path(vcf_indices), path(f1r2s), path(stats), path(tumor_pileups), path(normal_pileups)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(meta), path("${meta.somatic_name}_mutect_postproc.vcf.gz"), emit: mutect_vcf

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


    /*

        Call mutect 2 on a single interval shard. This process creates a new metadata block specific to
        the somatic caller of the form [somatic_name, somatic_caller, tumor_metadata, normal_metadata].

    */

    cpus 4
    memory "16GB"
    cache "lenient"

    container "broadinstitute/gatk:4.6.1.0"

    input:
    
        tuple val(somatic_name), val(tumor_meta), path(tumor_bam), path(tumor_bam_index),
                                 val(normal_meta), path(normal_bam), path(normal_bam_index)
        each path(interval_shard)
        path reference_fa
        path reference_index_dir
        path germline_resource
        path germline_resource_index
        path pon
        path pon_index_dir

    output:
        tuple val(somatic_meta),
                path("${somatic_name}_${interval_shard}_mutect.vcf.gz"), 
                path("${somatic_name}_${interval_shard}_mutect.vcf.gz.tbi"), 
                path("${somatic_name}_${interval_shard}_mutect_f1r2.tar.gz"), 
                path("*.stats"),
                emit: mutect_scatter_vcf

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
            --germline-resource $germline_resource \
            --panel-of-normals "${pon}" \
            --f1r2-tar-gz "${somatic_name}_${interval_shard}_mutect_f1r2.tar.gz" \
            -L $interval_shard \
            -O "${somatic_name}_${interval_shard}_mutect.vcf.gz" \
            --native-pair-hmm-threads $task.cpus

        """
}

process SPLIT_INTERVALS {

    /*

        Given a file of genomic intervals, split into scatter_count number of shards.

    */

    cpus 2
    memory "8GB"
    cache "lenient"

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


    /*

        Get transcript abundance estimations using kallisto quant.

    */

    cpus 8
    memory "32GB"

    conda "bioconda::kallisto=0.51.1"

    publishDir "${params.outdir}/${meta.somatic_sample}/rnaseq/kallisto/", mode: "copy"

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

    /*

    Create a kallisto index using a reference genome.

    */
    
    cpus 32
    memory "32GB"
    cache 'lenient'
    

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

    /*

    IN THE FUTURE THESE STEPS COULD BE SPLIT UP AND PARALLELIZED.

    Pre-process mapped BAM files according to GATK best practices. 
    Currently this calls:
        MarkDuplicatesSpark
        BaseRecalibrator
        ApplyBQSR
        GetPileupSummaries

    */



    cpus 16
    memory "48GB"
    cache "lenient"
    
    container "broadinstitute/gatk:4.6.1.0"

    publishDir "${params.outdir}/${meta.somatic_sample}/preprocess_bam/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), path(reads), path(reads_index)
        path reference_fa
        path reference_fa_index
        path known_sites_dbsnp
        path known_sites_dbsnp_index
        path known_sites_1000g_snps
        path known_sites_1000g_snps_index
        path known_indels
        path known_indels_index
        path mills
        path mills_index
        path common_germline
        path common_germline_index


    output:
        tuple val(meta), 
            path("${meta.sample_name}_${meta.molecule}_bqsr.bam"),
            path("${meta.sample_name}_${meta.molecule}_bqsr.bai"), 
        emit: preproc_bams
        
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_pileups.table"), emit: preproc_bams_pileups
        

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
            --known-sites $known_sites_dbsnp \
            --known-sites $known_sites_1000g_snps \
            --known-sites $known_indels \
            --known-sites $mills 
        
        gatk ApplyBQSR \
            -R $reference_fa \
            -I "${meta.sample_name}_${meta.molecule}_dedup.bam" \
            --bqsr-recal-file "${meta.sample_name}_${meta.molecule}_recal_table.table" \
            -O "${meta.sample_name}_${meta.molecule}_bqsr.bam" \
            --create-output-bam-index

        gatk GetPileupSummaries \
            -I "${meta.sample_name}_${meta.molecule}_bqsr.bam" \
            -V "${common_germline}" \
            -L "${common_germline}" \
            -O "${meta.sample_name}_${meta.molecule}_pileups.table"
        """
}

process BWA_MAP {

    /*
        Map fastq files using BWA. Outputs a sorted BAM file and its index
        Reads groups are created using metadata information and currently are basically just the same name.

    */

    cpus 16
    memory "32GB"
    cache "lenient"
    
    container "iarcbioinfo/bwa-mem2-tools:v1.0"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path reference_fa
        path reference_index_dir
        path bwa_index

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}.sam")
    
    script:
    """
    NEW_RG="@RG\\tID:${meta.sample_name}\\tSM:${meta.sample_name}\\tLB:${meta.sample_name}\\tPL:ILLUMINA"

    bwa-mem2 mem -t $task.cpus -R \$NEW_RG $reference_fa ${fastq1} ${fastq2} > "${meta.sample_name}_${meta.molecule}.sam"

    """
}

process BWA_POSTPROCESS {
    
    cpus 32
    memory "64GB"
    cache "lenient"

    conda "bioconda::samtools" 
    
    publishDir "${params.outdir}/${meta.somatic_sample}/alignment/bwa/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), val(bwa_sam)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_sorted.bam.bai"), emit: mapped_bam
    
    script:
    """

    samtools sort --threads $task.cpus $bwa_sam -o "${meta.sample_name}_${meta.molecule}_sorted.bam"
    samtools index "${meta.sample_name}_${meta.molecule}_sorted.bam"

    """
}


process CREATE_BWA_INDEX {

    /*
        SWITCH TO BWA2?
        Create a bwa index for bwa mem mapping.

    */

    cpus 16
    memory "32GB"
    
    container "iarcbioinfo/bwa-mem2-tools:v1.0"

    cache 'lenient'

    input:
        path reference_fa
        path reference_fa_index

    output:
        path "*{.bwt.2bit.64,.sa,.pac,.amb,.ann,.0123}", emit: bwa_index

    script:
        """
        bwa-mem2 index $reference_fa
        """




}

process STAR_FUSION {

    /*
        NOT IMPLEMENTED YET
        Run star-fusion to detect RNA fusion events.

    */

    cpus 4
    memory "64GB"

    container "trinityctat/starfusion:1.15.0"

    publishDir "${params.outdir}/${meta.somatic_sample}/fusions/", mode: "copy"

    input:
        tuple val(meta), path(chimeric_out), path(fastq1), path(fastq2)
        path ctat_resource_lib

    output:

        tuple val(meta), path("${meta.sample_name}_starfusion/*.fusion_predictions.tsv"), emit: fusion_preds
        tuple val(meta), path("${meta.sample_name}_starfusion/*.fusion_predictions.abridged.tsv"), emit: abridged_preds
        tuple val(meta), path("${meta.sample_name}_starfusion/*.coding_effect.tsv"), emit: coding_effect
        tuple val(meta), path("${meta.sample_name}_starfusion/"), emit: all_output


    script:
        """
        STAR-Fusion --genome_lib_dir $ctat_resource_lib \
             -J $chimeric_out \
             --examine_coding_effect \
             --FusionInspector validate \
             --left_fq $fastq1 \
             --right_fq $fastq2 \
             --denovo_reconstruct \
             --output_dir "./${meta.sample_name}_starfusion"

        """



}

process STAR_INDEX_BAM {

    /*

        Index the BAM file from a star process.

    */

    cpus 8
    memory "32GB"

    container "biocontainers/samtools:v1.9-4-deb_cv1"
    
    publishDir "${params.outdir}/${meta.somatic_sample}/alignment/star/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), path(bam)
        tuple val(meta), path(final_log)
        tuple val(meta), path(sj_out)
        tuple val(meta), path(chimeric_out), path(fastq1), path(fastq2)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_STAR_sorted.bam.bai"), emit: star_bam
        tuple val(meta), path(final_log), emit: final_log
        tuple val(meta), path(sj_out), emit: sj_out
        tuple val(meta), path(chimeric_out), path(fastq1), path(fastq2), emit: chimeric_out

    script:
        """
        samtools sort --threads $task.cpus  $bam -o "${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"
        samtools index -@ $task.cpus  "${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"

        """


}

process ARRIBA_FUSION {

    cpus 8

    memory "48GB"

    conda "bioconda::arriba=2.5.1"

    publishDir "${params.outdir}/${meta.somatic_sample}/fusions/arriba", mode: "copy"

    input:
        tuple val(meta), path(star_aligned_bam)
        path reference_fa
        path reference_fa_index
        path gtf
        path arriba_blacklist
        path arriba_known_fusions
        path arriba_protein_domains

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_arriba_fusions.tsv"), emit: arriba_fusions
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_arriba_fusions.discarded.tsv"), emit: discarded_fusions


    script:

        """

        arriba -x $star_aligned_bam \
            -g $gtf \
            -a $reference_fa \
            -b $arriba_blacklist \
            -k $arriba_known_fusions \
            -p $arriba_protein_domains \
            -o "${meta.sample_name}_${meta.molecule}_arriba_fusions.tsv" \
            -O "${meta.sample_name}_${meta.molecule}_arriba_fusions.discarded.tsv"
        """


}

process STAR_ALIGN {

    /*

    Align RNA reads with STAR. Parameters included from star-fusion to be able to use the output of this process
    in a downstream star-fusion or arriba process without having to re-map.

    */

    cpus 32

    memory "80GB"

    container "alexdobin/star:2.7.10a_alpha_220506"

    //publishDir "${params.outdir}/alignment/star_raw/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path(star_index_dir)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Aligned.out.bam"), emit: star_bam 
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Log.final.out"), emit:final_log
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_SJ.out.tab"), emit: sj_out
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Chimeric.out.junction"),path(fastq1), path(fastq2), emit: chimeric_out
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
            --chimSegmentMin 10 \
            --chimJunctionOverhangMin 10 \
            --outFilterMultimapNmax 50 \
            --chimOutJunctionFormat 1 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --alignIntronMax 100000 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outSAMattrRGline ID:"${meta.sample_name}" SM:"${meta.sample_name}" \
            --chimMultimapScoreRange 3 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimMultimapNmax 50 \
            --chimNonchimScoreDropMin 10 \
            --chimOutType Junctions WithinBAM HardClip \
            --chimScoreDropMax 30 \
            --peOverlapNbasesMin 10 \
            --peOverlapMMp 0.1 \
            --alignInsertionFlush Right \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --alignSplicedMateMapLmin 30 \
            --outFileNamePrefix ./${meta.sample_name}_${meta.molecule}_

        """

}


process CREATE_STAR_INDEX {

    /*
    
    Use a reference fasta and gtf to create a star index for star 2.7.10

    */

    cpus 32
    memory "64GB"
    cache 'lenient'

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

    /*

    Run fastp to do initial qc on fastq files. 

    Potentially add automated adapter trimming here.

    */
    
    cpus 8
    memory "16GB"

    conda "bioconda::fastp=1.0.1"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_R1_fastp.fastq.gz"), path("${meta.sample_name}_${meta.molecule}_R2_fastp.fastq.gz"), emit: fastqs
        path("${meta.sample_name}_*{html,json}*"), emit: reports
    
    tag "FastP on ${meta.sample_name} w/ ${meta.molecule}"

    publishDir "${params.outdir}/${meta.somatic_sample}/fastp_qc/${meta.sample_name}_${meta.molecule}_fastp/", mode: 'copy'

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


process HLAHD_HLA_CALLS {

    /*

    Parse the hla-hd final_result file to a .csv compatible with pvac (all hla alleles on single line)

    */

    cpus 1
    memory "2GB"

    conda "python=3.10 pandas=2.1"

    publishDir "${params.outdir}/${meta.somatic_sample}/HLA/pvac_input", mode:"copy"

    input:
        tuple val(meta), path(hla_result)

    output:
        tuple val(meta), path("${meta.sample_name}_hla_calls.csv")

    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import csv

    hlahd = pd.read_csv("$hla_result", sep = "\t", names = ["HLA", "Allele 1", "Allele 2"], nrows=21)
    hlahd = hlahd[(hlahd["Allele 1"] != "Not typed") & (hlahd["Allele 2"] != "Not typed")]

    allele_2_new = []
    allele_1_new = []
    for allele_1, allele_2 in zip(hlahd["Allele 1"], hlahd["Allele 2"]):
        allele_1_split = allele_1.split(":")
        allele_1 = allele_1_split[0] + ":" + allele_1_split[1]
        allele_1_new.append(allele_1)
                        
        if allele_2 == "-":
            allele_2_new.append(allele_1)
        else:
            allele_2_split = allele_2.split(":")
            allele_2 = allele_2_split[0] + ":" + allele_2_split[1]
            allele_2_new.append(allele_2)

    hlahd["Allele 1"] = allele_1_new
    hlahd["Allele 2"] = allele_2_new
    #hlahd = hlahd[hlahd["HLA"].isin(["A","B","C","DRB1","DQA1","DQB1"])]

    alleles = set()
    for allele, allele_1,allele_2 in zip(hlahd["HLA"],hlahd["Allele 1"], hlahd["Allele 2"]):
            if allele != "A" and allele != "B" and allele != "C":
                alleles.add(allele_1.split("-")[1])
                alleles.add(allele_2.split("-")[1])
            else:
                alleles.add(allele_1)
                alleles.add(allele_2)
                                                                    
    alleles = list(alleles)

    with open("${meta.sample_name}_hla_calls.csv", "w", newline="",encoding="utf-8") as f:
        writer = csv.writer(f,lineterminator="\\n")
        writer.writerow(alleles)

    """
        

}



process OPTITYPE_HLA_CALLS {

    /*

        Use the optitype tsv file to create pvac-compatible .csv (all hla alleles listed on single line)

    */
    
    cpus 1
    memory "4GB"

    conda "python=3.10 pandas=2.1"

    publishDir "${params.outdir}/${meta.somatic_sample}/HLA/", mode: "copy"

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

    /*

    Parse optitype output folder to get just the tsv and pdf output

    */

    cpus 1
    memory "4GB"

    publishDir "${params.outdir}/${meta.somatic_sample}/HLA/optitype/${meta.sample_name}_optitype", mode: "copy"

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

    /*
        HLA type patient's fastqs using OPTITYPE
        This process can be fairly memory intensive. 
        One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.

    */
    
    container "fred2/optitype:release-v1.3.1"

    cpus 8

    memory "150GB"

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

    /*

    Run HLA-HD to perform HLA typing on patient fastqs

    This has given me some trouble before with temp directories and issues where bowtie never initializes.
    This should be fixed by unzipping fastqs at the beginning of the script.
    One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.


    */
    
    cpus 6
    memory "24GB"

    container "griffithlab/hlahd:1.0"
    
    publishDir "${params.outdir}/${meta.somatic_sample}/HLA/hlahd/${meta.sample_name}_hlahd", mode: "copy"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
    
    output:
        tuple val(meta), path("./${meta.sample_name}/result/${meta.sample_name}_final.result.txt"), emit: hla_calls
        path("./${meta.sample_name}/result/")

    script:
        """
        ulimit -n 1024
        echo "\$TMPDIR"
        #mkdir -p tmp
        #mkdir -p /tmp/

        #export TMPDIR=/tmp/

        gzip -dc $fastq1 > fastq_r1.fastq
        gzip -dc $fastq2 > fastq_r2.fastq

        /opt/hlahd/bin/hlahd.1.6.1.sh \
            -f /opt/hlahd/freq_data \
            -t $task.cpus \
            fastq_r1.fastq fastq_r2.fastq \
            /opt/hlahd/HLA_gene.split.txt \
            /opt/hlahd/dictionary \
            "${meta.sample_name}" \
            ./
        """


}

process PVACSEQ {

    /*

    Run PVACseq on a single sample.

    Inputs:
        *(See PVACseq input preparation documents for more info.)
        Somatic VCF
        Phased Germline VCF
        Tumor Sample Metadata
        Normal Sample Metadata

    Output:
        Pvacseq folder containing:
            Combined neoantigen predictions
            MHC I neoantigen predictions
            MHC II neoantigen predictions

    */

    cpus 8
    memory "90GB"

    container "griffithlab/pvactools:6.0.3"

    publishDir "${params.outdir}/${somatic_meta.somatic_name}/pvactools/", mode: "copy"

    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            path(phased_vcf), path(phased_vcf_index),
            path(hla_pvac_input)
    output:
        path("${somatic_meta.somatic_name}_pvacseq"), emit: pvacseq_dir

    script:
        """
        pvacseq run \
            $somatic_vcf \
            ${somatic_meta.tumor_metamap.sample_name} \
            \$(head $hla_pvac_input -n 1) \
            all \
            "${somatic_meta.somatic_name}_pvacseq" \
            -e1 8,9,10,11 \
            -e2 12,13,14,15,16,17,18 \
            --phased-proximal-variants-vcf $phased_vcf \
            --normal-sample-name ${somatic_meta.normal_metamap.sample_name} \
            --iedb-install-directory /opt/iedb \
            --pass-only \
            -t $task.cpus
        """
}

process PVACFUSE {
    
    cpus 8
    memory "90GB"

    container "griffithlab/pvactools:6.0.3"

    publishDir "${params.outdir}/${sample_meta.somatic_sample}/pvactools/", mode: "copy"

    input:
        tuple val(sample_meta), path(arriba_fusions), path(hla_pvac_input), path(starfusion_calls)

    output:
        path("${sample_meta.somatic_name}_pvacseq"), emit: pvacseq_dir
    
    script:
        """
        pvacfuse run \
            $arriba_fusions \
            "${sample_meta.somatic_sample}" \
            \$(head $hla_pvac_input -n 1) \
            all \
            "${sample_meta.somatic_sample}_pvacfuse" \
            --starfusion-file $starfusion_calls \
            -e1 8,9,10,11 \
            -e2 12,13,14,15,16,17,18 \
            --iedb-install-directory /opt/iedb \
            -t $task.cpus
        """

}

process PVACSEQ_SELECT_ALGOS {

    /*

    Run PVACseq on a single sample.

    Inputs:
        *(See PVACseq input preparation documents for more info.)
        Somatic VCF
        Phased Germline VCF
        Tumor Sample Metadata
        Normal Sample Metadata

    Output:
        Pvacseq folder containing:
            Combined neoantigen predictions
            MHC I neoantigen predictions
            MHC II neoantigen predictions

    */

    cpus 8
    memory "90GB"

    container "griffithlab/pvactools:6.0.3"

    publishDir "${params.outdir}/${somatic_meta.somatic_name}/pvactools/", mode: "copy"

    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            path(phased_vcf), path(phased_vcf_index),
            path(hla_pvac_input)
    output:
        path("${somatic_meta.somatic_name}_pvacseq_select_algos"), emit: pvacseq_dir

    script:
        """
        pvacseq run \
            $somatic_vcf \
            ${somatic_meta.tumor_metamap.sample_name} \
            \$(head $hla_pvac_input -n 1) \
            MHCflurry MHCflurryEL MHCnuggetsI MHCnuggetsII  NetMHCIIpan NetMHCIIpanEL NetMHCpan NetMHCpanEL \
            "${somatic_meta.somatic_name}_pvacseq_select_algos" \
            -e1 8,9,10,11 \
            -e2 12,13,14,15,16,17,18 \
            --phased-proximal-variants-vcf $phased_vcf \
            --normal-sample-name ${somatic_meta.normal_metamap.sample_name} \
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
                    somatic_sample: row.somatic_sample,
                    sample_name: row.sample_name,
                    sample_type: row.sample_type,
                    molecule: row.molecule,
                    sequencing_type: row.sequencing_type,
                    hla: row.hla
                ]
                
                reads = [
                    file(row.fastqr1, checkIfExists: true),
                    file(row.fastqr2, checkIfExists: true)
                ]
            return [meta, reads]
        }
    

    // PULL REFERENCE FASTA AND REFERENCE FASTA SUPPLEMENTAL FILES
    
    reference_fa = Channel.fromPath(params.reference_fa).first()
    reference_index_files = Channel.fromPath(params.reference_index_dir).collect()
    reference_dict = Channel.fromPath(params.reference_dict).first()
    
    common_germline = Channel.fromPath(params.common_germline).first()
    common_germline_index = Channel.fromPath(params.common_germline_index).first()
    
    known_sites_dbsnp = Channel.fromPath(params.known_sites_dbsnp).first()
    known_sites_dbsnp_index = Channel.fromPath(params.known_sites_dbsnp_index).first()

    known_sites_1000g_snps = Channel.fromPath(params.known_sites_1000g_snps).first()
    known_sites_1000g_snps_index = Channel.fromPath(params.known_sites_1000g_snps_index).first()

    known_indels = Channel.fromPath(params.known_indels).first()
    known_indels_index = Channel.fromPath(params.known_indels_index).first()
    
    mills = Channel.fromPath(params.mills).first()
    mills_index = Channel.fromPath(params.mills_index).first()

    gnomad = Channel.fromPath(params.gnomad).first()
    gnomad_index = Channel.fromPath(params.gnomad_index).first()
    
    pon = Channel.fromPath(params.pon).first()
    pon_index = Channel.fromPath(params.pon_index).first()
    
    hapmap = Channel.fromPath(params.hapmap).first()
    hapmap_index = Channel.fromPath(params.hapmap_index).first()
    
    intervals_file = Channel.fromPath(params.intervals_file).first()
    
    kallisto_reference = Channel.fromPath(params.kallisto_reference).first()
    gencode_gtf = Channel.fromPath(file(params.gencode_gtf)).first()
    
    ctat_resource_dir = Channel.fromPath(file(params.ctat_resource_dir)).first()
   
    arriba_blacklist = Channel.fromPath(file(params.arriba_blacklist)).first()
    arriba_known_fusions = Channel.fromPath(file(params.arriba_known_fusions)).first()
    arriba_protein_domains = Channel.fromPath(file(params.arriba_protein_domains)).first()

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

    hla_calls_hlahd = HLAHD_HLA_CALLS(hla_hlahd.hla_calls)
    
    hla_calls = hla_calls_hlahd.map { meta, hla_call ->
                                      def hla_value = meta.hla != "CALL" ? meta.hla : hla_call
                                      return [meta,hla_value]
                                    }.filter { meta, hla_call ->
                                               meta.sample_type == "Normal"
                                        }


    // ALIGNMENT OF DNA SEQUENCING USING BWA

    if (params.bwa_index) {
        bwa_index = Channel.fromPath(params.bwa_index).collect()
    } else {
        bwa_index = CREATE_BWA_INDEX(reference_fa, reference_index_files)
    }

    bwa_sam = BWA_MAP(fastp_by_molec.dna, reference_fa, reference_index_files, bwa_index)
    bwa_mapped = BWA_POSTPROCESS(bwa_sam)


    // ALIGNMENT OF RNA SEQUENCING USING STAR
    // STAR alignment process includes chimeric reads for directly going from BAM -> Star-Fusion
    
    if (params.star_index) {
        star_index = Channel.fromPath(params.star_index)
    }
    else {
        star_index = CREATE_STAR_INDEX(reference_fa, reference_index_files, gencode_gtf).star_index
    }
    
    star_rna_align = STAR_ALIGN(fastp_by_molec.rna, params.star_index)
    star_rna = STAR_INDEX_BAM(star_rna_align)

    star_fusion = STAR_FUSION(star_rna_align.chimeric_out, ctat_resource_dir)
    
    arriba_fusions = ARRIBA_FUSION(star_rna_align.star_bam,
                                    reference_fa,
                                    reference_index_files,
                                    gencode_gtf,
                                    arriba_blacklist,
                                    arriba_known_fusions,
                                    arriba_protein_domains)

    
    
    //PVACFUSE(

    // BAM PREPROCESSING OF DNA: GATK BEST PRACTICES

    preprc = PREPROCESS_BAM(bwa_mapped,
                            reference_fa, reference_index_files,
                            known_sites_dbsnp, known_sites_dbsnp_index,
                            known_sites_1000g_snps, known_sites_1000g_snps_index,
                            known_indels, known_indels_index,
                            mills, mills_index,
                            common_germline, common_germline_index)


    
    // RNA: TRANSCRIPT ABUNDANCE ESTIMATION
    
    if (params.kallisto_index) {
        kallisto_index = params.kallisto_index
    } else {
        kallisto_index = KALLISTO_INDEX(kallisto_reference)

    }
    kallisto = KALLISTO_QUANT(fastp_by_molec.rna,
                              kallisto_index)



    // INTERVAL CREATION FOR SOMATIC CALLERS

    intervals = SPLIT_INTERVALS(reference_fa, reference_index_files,
                                intervals_file, params.scatter_count)


    
    // PROBLEMATIC ?

    preproc_bams_type_branched = preprc.preproc_bams
                                   | branch {meta, bam, bai ->
                                        normal: meta.sample_type == "Normal"
                                        tumor: meta.sample_type == "Tumor"
    }
    
    //preproc_bams_pileups_type_branched = preprc.preproc_bams_pileups
    //                                        | branch {meta, bam, bai, pileups ->
    //                                            normal: meta.sample_type == "Normal"
    //                                            tumor: meta.sample_type == "Tumor"
    //                                        }
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
    vt_germline = VT_POSTPROCESS_GERMLINE(germline_postprocess,
                                          reference_fa,
                                          reference_index_files)
    
    vep_germline = VEP_ANNOTATE_GERMLINE(vt_germline,
                                         reference_fa,
                                         params.vep_cache,
                                         params.vep_plugins)


    final_germline = INDEX_FINAL_VCF_GERMLINE(vep_germline)

    preproc_bams_by_sample = preprc.preproc_bams.map { meta, bam, bai -> tuple(meta.somatic_sample, [meta, bam, bai]) } 
                                                     | groupTuple
                                                     | map { somatic_id, samples ->
                                                                def tumor  = samples.find { it[0].sample_type == 'Tumor' }
                                                                def normal = samples.find { it[0].sample_type == 'Normal' }
                                                                                     
                                                                return [ somatic_id, 
                                                                         tumor[0], tumor[1], tumor[2],
                                                                         normal[0], normal[1], normal[2]]
                                                     }
    
    // SOMATIC VARIANT CALLING

    // MUTECT2
    // Run Mutect2 with scatter/gather approach
   
    mutect_scattered = MUTECT2_SCATTER(preproc_bams_by_sample,
                    intervals.flatten(),
                    reference_fa, reference_index_files,
                    gnomad, gnomad_index,
                    pon, pon_index) 
   
    mutect_gathered = mutect_scattered 
                        | groupTuple 
                        | map { meta, vcf, vcf_index, f1r2, stats -> tuple(meta.somatic_name, [meta, vcf, vcf_index, f1r2, stats]) }
    
    pileups = preprc.preproc_bams_pileups.map { meta, pileups -> tuple(meta.somatic_sample, [meta, pileups]) } 
                                         | groupTuple
                                         | map { somatic_name, samples ->
                                                 def tumor = samples.find { it[0].sample_type == "Tumor" }
                                                 def normal = samples.find { it[0].sample_type == "Normal" }

                                             return [somatic_name, tumor[1], normal[1] ]
                                     }



    mutect_scattered_pileups = mutect_gathered
        | join(pileups)
        | map { somatic_name, inner, tumor_pileups, normal_pileups ->
                def (meta, vcfs, tbis, f1r2s, stats) = inner
                tuple(somatic_name, meta, vcfs, tbis, f1r2s, stats, tumor_pileups, normal_pileups)
        }

    
    // Postprocess Mutect2 Output using pileups, stats, etc.
    mutect_postprocess = POSTPROCESS_MUTECT2_SCATTER(mutect_scattered_pileups,
                                                     reference_fa,
                                                     reference_index_files)

    
    // STRELKA
    strelka = STRELKA(preproc_bams_by_sample, reference_fa, reference_index_files)

    strelka_postprocess = POSTPROCESS_STRELKA(strelka)
    
    strelka = ADD_VCF_GT_FIELD(strelka_postprocess.strelka_vcf)

    somatic_vcfs = mutect_postprocess.concat(strelka)
    
    somatic_filtered_vcfs = FILTER_VCF(somatic_vcfs).filter_vcf
   
    // POST-PROCESS SOMATIC VCFS

    // Variant Decomposition & Normalization
    somatic_vcfs_vt = VT_SOMATIC_POSTPROCESS(somatic_filtered_vcfs, reference_fa, reference_index_files)
        | map { meta, vcf, tbi  -> tuple(meta.somatic_name, [meta, vcf, tbi]) }
        | groupTuple
        

    merged_vcf = MERGE_SOMATIC_VCFS(somatic_vcfs_vt,
                                      reference_fa,
                                      reference_index_files)

    

    // PVACtools VCF PREPARATION
    
    // Annotate VCF with VEP
    vep_annot = VEP_ANNOTATE(merged_vcf,
                       reference_fa,
                       params.vep_cache,
                       params.vep_plugins)

    vep = VEP_FILTER(vep_annot, params.vep_cache, params.vep_plugins)


    
    vep = vep | map { meta, vcf -> tuple(meta.tumor_metamap.sample_name, [meta, vcf]) } | groupTuple

    star_rna = star_rna.star_bam | map { meta, bam, bai -> tuple(meta.sample_name, [meta, bam, bai]) } | groupTuple
   

    vep_star = vep
        | join(star_rna)
        | map { id, vcf_info, bam_info ->
            def (vcf_meta, vcf) = vcf_info[0]
            def (star_meta, star_bam, star_bai) = bam_info[0]
            return [vcf_meta, vcf, star_meta, star_bam, star_bai ]
        }
    

    // Add Read Coverage to VCF with Bamreadcount
    bamreadcount = BAMREADCOUNT(vep_star,
                                reference_fa)
    
    vcf_annotated_coverage = ANNOTATE_VCF_COVERAGE(bamreadcount)

    
    // Add transcript abundance estimation from kallisto
    

    kallisto = kallisto | map { meta, abundance -> tuple(meta.sample_name, [meta, abundance]) } | groupTuple
    vcf_annot = vcf_annotated_coverage | map { meta, vcf -> tuple(meta.tumor_metamap.sample_name, [meta, vcf]) } | groupTuple
    
    vcf_annotated = vcf_annot
        | join(kallisto)
        | map { id, vcf_info, kallisto_info ->
                def (vcf_meta, vcf) = vcf_info[0]
                def (kallisto_meta, abundance) = kallisto_info[0]
                return [vcf_meta, vcf, kallisto_meta, abundance ]
        }

    vcf_annotated_expression = ANNOTATE_VCF_EXPRESSION(vcf_annotated)

    
    // Index final somatic vcf for input to PVACseq
    vcf_final = INDEX_FINAL_VCF(vcf_annotated_expression)
    
    
    vcf_final_grouped = vcf_final
        | map { meta, vcf, vcf_index -> tuple(meta.normal_metamap.sample_name, [meta, vcf, vcf_index]) }
        | groupTuple
    germline_grouped = vt_germline
        | map { meta, vcf, vcf_index -> tuple(meta.sample_name, [meta, vcf, vcf_index]) }
        | groupTuple


    vcf_final_somatic_germline_grouped = vcf_final_grouped
        | join(germline_grouped)
        | map { id, somatic_vcf_tuple, germline_vcf_tuple ->
                def (somatic_meta, somatic_vcf, somatic_vcf_index) = somatic_vcf_tuple[0]
                def (germline_meta, germline_vcf, germline_vcf_index) = germline_vcf_tuple[0]
                return [somatic_meta, somatic_vcf, somatic_vcf_index, germline_meta, germline_vcf, germline_vcf_index ]
        }

    



    // PERFORM VCF PHASING USING GERMLINE CALLS

    // Create Tumor-Only VCF From Final Somatic VCF (vcf_final)
    vcf_phase_select_variants = PHASE_VCF_SELECT_VARIANTS(vcf_final_somatic_germline_grouped,
                           reference_fa,
                           reference_index_files)
    
    vcf_phase_rename_samples = PHASE_VCF_RENAME(vcf_phase_select_variants)
    

    // Combine Tumor-Only VCF with germline variants
    vcf_phase_combine = PHASE_VCF_COMBINE_VARIANTS(vcf_phase_rename_samples,
                                                   reference_fa,
                                                   reference_index_files)
    // Sort combined VCF
    vcf_phase_sort = PHASE_VCF_SORT_VCF(vcf_phase_combine, reference_dict)
    
    
    // Call GATK ReadBackedPhasing (required gatk 3.6.0) to phase VCF
    
    vcf_phase_sort_group = vcf_phase_sort
        | map {meta, vcf -> tuple(meta.tumor_metamap.sample_name, [meta, vcf]) }
        | groupTuple

    preproc_bams_type_branched_group = preproc_bams_type_branched.tumor
        | map { meta, bam, bai -> tuple(meta.sample_name, [meta, bam, bai]) }
        | groupTuple

    vcf_phase_sort_group_bams = vcf_phase_sort_group 
        | join(preproc_bams_type_branched_group)
        | map { id, phased_vcf_tuple, bam_tuple ->
                def (vcf_meta, vcf) = phased_vcf_tuple[0]
                def (bam_meta, bam, bai) = bam_tuple[0]
                return [vcf_meta, vcf, bam_meta, bam, bai]
        }
    
    vcf_phase_rbphase = PHASE_VCF_RBPHASING(vcf_phase_sort_group_bams,
                                            reference_fa, reference_index_files)
    
    /// VEP Annotation phased vcf
    vcf_phase_vep = PHASE_VCF_VEP(vcf_phase_rbphase,
                                    reference_fa, 
                                    params.vep_cache, params.vep_plugins)        
    // Zip and index phased vcf
    vcf_phased = PHASE_VCF_INDEX(vcf_phase_vep)
    
    final_phased = vcf_phased
        | map { meta, vcf, vcf_index -> tuple(meta.somatic_name, [meta, vcf, vcf_index]) }

    final_somatic = vcf_final
        | map { meta, vcf, vcf_index -> tuple(meta.somatic_name, [meta, vcf, vcf_index]) }

    final_hla = hla_calls
        | map { meta, hla_calls -> tuple(meta.somatic_sample, [meta, hla_calls]) }
    
    pvacseq_input = final_somatic
        | join(final_phased)
        | join(final_hla)
        | map {id, somatic_tuple, phased_tuple, hla_tuple ->
            def (somatic_meta, somatic_vcf, somatic_vcf_index) = somatic_tuple
            def (phased_meta, phased_vcf, phased_vcf_index) = phased_tuple
            def (hla_meta, hla_calls) = hla_tuple
            return [somatic_meta, somatic_vcf, somatic_vcf_index, phased_vcf, phased_vcf_index, hla_calls ]
        }

            

    PVACSEQ = PVACSEQ(pvacseq_input)


    final_fusions = arriba_fusions.arriba_fusions
        | map { meta, fusion_tsv -> tuple(meta.somatic_sample, [meta, fusion_tsv]) }
    final_starfusions = star_fusion.fusion_preds
        | map { meta, star_fusion_pred -> tuple(meta.somatic_sample, [meta, star_fusion_pred]) }

    pvacfuse_input = final_fusions | join(final_hla) | join(final_starfusions)
        | map { id, arriba_fusion_tuple, hla_tuple, starfusion_tuple ->
                def (arriba_meta, arriba_fusions) = arriba_fusion_tuple
                def (hla_meta, hla_calls) = hla_tuple
                def (starfusion_meta, starfusion_calls) = starfusion_tuple
                return [ arriba_meta, arriba_fusions, hla_calls, starfusion_calls ]
        }


    PVACFUSE = PVACFUSE(pvacfuse_input)

}
