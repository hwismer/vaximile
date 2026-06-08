
process HAPLOTYPE_CALLER_FILTER_VARIANTS {
    
    cpus 4
    memory "32GB"
    container "broadinstitute/gatk:4.3.0.0"

    tag "Filtering germline variants in ${meta.sample_name}"

    input:
        tuple val(meta), path(vcf)
        tuple path(hapmap), path(hapmap_index)
        tuple path(mills), path(mills_index)

    output:
        tuple val(meta), val("haplotypecaller"), path("${meta.sample_name}_germline_filtered.vcf.gz"), path("${meta.sample_name}_germline_filtered.vcf.gz.tbi"), emit: germline_vcf
    
    script:
    """
    gatk FilterVariantTranches \
        -V $vcf \
        --resource $hapmap \
        --resource $mills \
        --info-key CNN_1D \
        --snp-tranche 99.95 \
        --indel-tranche 99.4 \
        -O ${meta.sample_name}_germline_filtered.vcf.gz \
        --create-output-variant-index
    """
}


process HAPLOTYPE_CALLER_GATHER_VCFS {

    cpus 2
    memory "8GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Gathering haplotype caller VCFs for ${sample_meta.sample_name}"

    input:
        tuple val(sample_meta), path(vcfs)

    output:
        tuple val(sample_meta), path("${sample_meta.sample_name}_merged.vcf")

    script:
    
    def sorted_vcfs = vcfs.sort { a, b -> a.name <=> b.name }
    def vcf_as_input = sorted_vcfs.collect { vcf ->
            "--INPUT ${vcf}"
        }.join(' ')

    """
    gatk GatherVcfs \
        $vcf_as_input \
        -O "${sample_meta.sample_name}_merged.vcf"
    
    """

}

process HAPLOTYPE_CALLER_GATHER_SELECT_VARIANTS {

    cpus 2
    memory "8GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Select variants from ${sample_meta.sample_name} within ${interval_shard}"

    input:
        tuple val(sample_meta), path(vcf), path(vcf_index), path(interval_shard)

    output:
        tuple val(sample_meta), path("${sample_meta.sample_name}_${interval_shard}.vcf.gz")

    script:
    
    """
    gatk SelectVariants \
        -V $vcf \
        -L $interval_shard \
        -O "${sample_meta.sample_name}_${interval_shard}.vcf.gz"
    """
}


process HAPLOTYPE_CALLER_SCATTER {

    /*

    Use HaplotypeCaller on a single scattered interval. Post processes with CNNScoreVariants.

    */

    cpus 4
    memory "16GB"
    container "broadinstitute/gatk:4.3.0.0"

    tag "Haplotype Caller on ${meta.sample_name} on ${interval_shard}"

    input:
        tuple val(meta), path(bam), path(bai), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict
        val(interval_padding)

    output:
        tuple val(meta), path("${meta.sample_name}_${interval_shard}.vcf.gz"), path("${meta.sample_name}_${interval_shard}.vcf.gz.tbi"),  path(interval_shard)

    script:
        """
        gatk HaplotypeCaller \
            -R $reference_fa \
            -I $bam \
            -L $interval_shard \
            -O "${meta.sample_name}_${interval_shard}.vcf.gz" \
            -ERC NONE \
            --native-pair-hmm-threads $task.cpus \
            -ip $interval_padding \
            --create-output-variant-index
        """
}

process HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS {
    
    /*

    Use HaplotypeCaller on a single scattered interval. Post processes with CNNScoreVariants.

    */

    cpus 4
    memory "16GB"

    container "broadinstitute/gatk:4.3.0.0"

    tag "Scoring haplotypecaller variants from ${meta.sample_name} on ${interval_shard}"

    input:
        tuple val(meta), path(vcf), path(vcf_index), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict
        val(interval_padding)

    output:
        tuple val(meta), path("${meta.sample_name}_${interval_shard}_CNN.vcf.gz"), path("${meta.sample_name}_${interval_shard}_CNN.vcf.gz.tbi"), path(interval_shard)

    script:
        """
        gatk CNNScoreVariants \
            -V $vcf \
            -L $interval_shard \
            -ip $interval_padding \
            -R $reference_fa \
            --create-output-variant-index \
            -O "${meta.sample_name}_${interval_shard}_CNN.vcf.gz"
        """
}


process DEEPVARIANT {

    cpus 16
    memory "32GB"
    container "google/deepvariant:1.10.0"
    
    tag "DeepVariant on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai), path(bed)
        tuple path(reference_fa), path(reference_index)

    output:
        tuple val(meta), val("deepvariant"), path("${meta.sample_name}_deepvariant.vcf.gz"), path("${meta.sample_name}_deepvariant.vcf.gz.tbi"), emit: vcf
        tuple val(meta), val("deepvariant"), path("${meta.sample_name}_deepvariant.gvcf.gz"), path("${meta.sample_name}_deepvariant.gvcf.gz.tbi"), emit: gvcf

    script:

	def model_map = [
        exome      : 'WES',
        genome     : 'WGS',
    ]

	def model = model_map[meta.sequencing_type]

    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=$model \
        --ref=$reference_fa \
        --reads=$bam \
        --output_vcf=${meta.sample_name}_deepvariant.vcf.gz \
        --output_gvcf=${meta.sample_name}_deepvariant.gvcf.gz \
        --num_shards=$task.cpus \
        --vcf_stats_report=true \
        --regions $bed \
        --disable_small_model=true
    """
}


process STRELKA_GERMLINE {

    /*

    Run Strelka in germline mode

    */

    cpus 8
    memory "24GB"
    container "mgibio/strelka:2.9.9"

    tag "Running Strelka in germline mode on ${meta.sample_name}"

   input:
        tuple val(meta), path(bam), path(bai), path(bed), path(bed_index)
        tuple path(reference_fa), path(reference_fai)

    output:
        tuple val(meta), val("strelka"), path("./strelka/results/variants/variants.vcf.gz"), path("./strelka/results/variants/variants.vcf.gz.tbi")

    script:

    def exome_flag = (meta.sequencing_type == "exome" || meta.sequencing_type == "exome_ffpe") ? "--exome" : ""

    """
    /opt/strelka/bin/configureStrelkaGermlineWorkflow.py \
        --bam $bam \
        ${exome_flag} \
        --callRegions $bed \
        --referenceFasta $reference_fa \
        --runDir "./strelka"

    ./strelka/runWorkflow.py -m local -j $task.cpus
    """

}

process FILTER_VCF {

    /*

    Filter VCF, retaining only PASS variants.

    */

    cpus 4
    memory "16GB"
    container "staphb/bcftools:1.23"
    
    tag "Filtering VCF $vcf"

    input:
        tuple val(meta), val(caller), path(vcf), path(tbi)

    output:
        tuple val(meta), val(caller), path("${meta.sample_name}_${caller}_filtered_variants.vcf.gz"), path("${meta.sample_name}_${caller}_filtered_variants.vcf.gz.tbi"), emit: filtered_vcf


    script:
        """
        bcftools view -f PASS -Oz -o "${meta.sample_name}_${caller}_filtered_variants.vcf.gz" $vcf
        bcftools index -t "${meta.sample_name}_${caller}_filtered_variants.vcf.gz"
        """



}

process POSTPROCESS_VCF {

    /*

        Postprocess somatic variants from a vcf. This is used to normalize variants and decompose biallelics
        into different entries. The vcf is then sorted and duplicate entries removed. The final vcf is then indexed.

    */

    cpus 4
    memory "16GB"
    container "staphb/bcftools:1.23"

    tag "Normalizing $vcf"

    input:
        tuple val(meta), val(caller), path(vcf), path(vcf_index)
        tuple path(reference_fa), path(reference_index)

    output:
        tuple val(meta), val(caller),
            path("${meta.sample_name}_${caller}_variants.vcf.gz"),
            path("${meta.sample_name}_${caller}_variants.vcf.gz.tbi"), emit: postproc_vcf

    script:

        """
        bcftools norm -m -any -d exact -f $reference_fa $vcf -Oz -o norm_vcf.vcf.gz
        bcftools sort norm_vcf.vcf.gz -Oz -o "${meta.sample_name}_${caller}_variants.vcf.gz"
        bcftools index -t "${meta.sample_name}_${caller}_variants.vcf.gz"
        """
}

process MERGE_GERMLINE_VCFS {

    /*

    Use the deprecated CombineVariants from GATK 3.6.0 to combine vcf files.

    */

    cpus 4
    memory "16GB"
    container "broadinstitute/gatk3:3.6-0"

    tag "Combining VCFs: $vcf1 $vcf2 $vcf3"

    input:
        tuple val(sample_meta),
            val(vcf1_caller), path(vcf1), path(vcf1_index),
            val(vcf2_caller), path(vcf2), path(vcf2_index),
            val(vcf3_caller), path(vcf3), path(vcf3_index)
        tuple path(reference_fa), path(reference_index)
        path(reference_dict)

    output:
        tuple val(sample_meta), path("${sample_meta.sample_name}_germline_variants.vcf.gz")


    script:

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
            -R $reference_fa \
            -genotypeMergeOptions PRIORITIZE \
            --rod_priority_list $vcf1_caller,$vcf2_caller,$vcf3_caller \
            -V:$vcf1_caller $vcf1 \
            -V:$vcf2_caller $vcf2 \
            -V:$vcf3_caller $vcf3 \
            --minimumN 2 \
            -o "${sample_meta.sample_name}_germline_variants.vcf.gz"
        """

}

process VEP_ANNOTATE {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    cpus 6
    memory "16GB"
    container "ensemblorg/ensembl-vep:release_115.0"

    tag "VEP on ${sample_meta.sample_name}"
    input:
        tuple val(sample_meta), path(vcf)
        tuple path(reference_fa), path(reference_index)
        path vep_cache
        path vep_plugins

    output:
        tuple val(sample_meta), path("${sample_meta.sample_name}_vep.vcf"), emit: vcf
        tuple val(sample_meta), path("*.html"), emit: report

    script:
        """
        vep \
            --input_file $vcf  \
            --output_file ${sample_meta.sample_name}_vep.vcf \
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
