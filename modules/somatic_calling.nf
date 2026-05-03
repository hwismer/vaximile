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
        path intervals

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

