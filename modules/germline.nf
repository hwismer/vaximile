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
