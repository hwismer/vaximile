// Helper processes used in the main workflow
include { SPLIT_INTERVALS } from "../../modules/local/split_intervals/main"
include { CAPTURE_KIT_BED_PROCESS } from "../../modules/local/capture_kit_bed_process/main"
include { PULL_VEP_PVAC_PLUGINS } from "../../modules/local/pull_vep_pvac_plugins/main"
include { PULL_VEP_CACHE } from "../../modules/local/pull_vep_cache/main"
include { PULL_CTAT_RESOURCE_BUNDLE } from "../../modules/local/pull_ctat_resource_bundle/main"
include { PULL_ARRIBA_RESOURCES } from "../../modules/local/pull_arriba_resources/main"
include { PULL_ASCAT_RESOURCES } from "../../modules/local/pull_ascat_resources/main"
include { MULTIQC } from "../../modules/local/multiqc/main"
include { ASCAT } from "../../modules/local/ascat/main"
include { dedupe_libraries; fan_out_pairs; pair_tumor_normal } from "../../subworkflows/pipeline_init.nf"

// Subworkflows
include { DNA_QC_WORKFLOW } from "../../subworkflows/dna_preprocessing.nf"
include { RNA_QC_WORKFLOW } from "../../subworkflows/rnaseq.nf"
include { SOMALIER } from "../../subworkflows/dna_preprocessing.nf"
include { DNA_ALIGN_AND_PREPROC } from "../../subworkflows/dna_preprocessing.nf"
include { BWA_INDEX } from "../../subworkflows/dna_preprocessing.nf"
include { PREPARE_REFERENCE_FASTA } from "../../subworkflows/reference.nf"
include { HLA_TYPING_WORKFLOW } from "../../subworkflows/hla.nf"
include { HLA_LOH_WORKFLOW } from "../../subworkflows/hla.nf"
include { ASCAT_MQC } from "../../modules/local/ascat_mqc/main"
include { RNASEQ_WORKFLOW } from "../../subworkflows/rnaseq.nf"
include { MUTECT2 } from "../../subworkflows/somatic_variant_calling.nf"
include { SOMATIC_CONSENSUS } from "../../subworkflows/somatic_variant_calling.nf"
include { STRELKA_WORKFLOW } from "../../subworkflows/somatic_variant_calling.nf"
include { GERMLINE_WORKFLOW } from "../../subworkflows/germline_variant_calling.nf"
include { FUSION_CALLING } from "../../subworkflows/rnaseq.nf"
include { PVAC_INPUT_PREP_WORKFLOW } from "../../subworkflows/pvactools.nf"
include { PVACTOOLS_WORKFLOW } from "../../subworkflows/pvactools.nf"
include { DEEPSOMATIC_WORKFLOW } from "../../subworkflows/somatic_variant_calling.nf"

//**************************************************************************************************************************************

// Make a (vcf, vcf index) channel automatically given a vcf path. Tbi must also exist

def make_vcf_channel(vcf_param) {
    return vcf_param
        ? Channel.value(tuple(file(vcf_param), file("${vcf_param}.tbi")))
        : Channel.empty()
}

//**********************************************************************************************************************

workflow VAXIMILE {

    take:
    ch_samplesheet   // channel: [ meta, [ fastq_1, fastq_2 ] ]
    ch_capture_kits  // channel: [ kit_name, bed ]

    main:

    //**************************************************************************************************************************************

    reference_fa = Channel.fromPath(params.reference_fa)
    prepared_reference = PREPARE_REFERENCE_FASTA(reference_fa)
    reference_genome = prepared_reference.fa_fai_pair.first()
    reference_dict = prepared_reference.dict.first()
    

    // BWA_INDEX handles both cases: it validates and loads a prebuilt index when
    // --bwa_index is given, and builds one otherwise. Either way it emits the same shape.
    bwa_index = BWA_INDEX(reference_genome, params.bwa_index).bwa_index
    
    transcriptome_reference = Channel.fromPath(params.transcriptome_reference).first()
    gtf = Channel.fromPath(file(params.gtf)).first()
    human_ref_peptides = Channel.fromPath(file(params.human_ref_peptides)).first()

    common_germline = make_vcf_channel(params.common_germline)
    known_sites_dbsnp = make_vcf_channel(params.known_sites_dbsnp)
    known_sites_1000g_snps = make_vcf_channel(params.known_sites_1000g_snps)
    known_indels = make_vcf_channel(params.known_indels)
    mills = make_vcf_channel(params.mills)
    gnomad = make_vcf_channel(params.gnomad)
    pon = make_vcf_channel(params.pon)
    hapmap = make_vcf_channel(params.hapmap)
    
    somalier_sites = Channel.fromPath(params.somalier_sites).first()
    
    hla_fasta = Channel.fromPath(params.hla_fasta).first()
    hla_fasta_fai = Channel.fromPath(params.hla_fasta_fai).first()
    hla_reference = hla_fasta.combine(hla_fasta_fai)
    hla_kmers = Channel.fromPath(params.hla_kmers).first()
    hla_freq = Channel.fromPath(params.hla_freq).first()
    



    //**************************************************************************************************************************************
    // Pull Resources

    ctat_bundle = PULL_CTAT_RESOURCE_BUNDLE().ctat_resource_dir
    arriba_resources = PULL_ARRIBA_RESOURCES().resources
    vep_plugins = PULL_VEP_PVAC_PLUGINS().plugins

    // --vep_cache is optional: use the directory given, otherwise pull the release-115
    // GRCh38 cache once into ./vaximile_resources/vep_cache and reuse it on later runs.
    //
    // MUST be a value channel either way. VEP runs once per VCF, and a queue channel
    // holding a single item would be consumed by the first of those tasks, leaving the
    // rest with no cache. A process output that emits exactly once is already a value
    // channel, so both branches broadcast.
    vep_cache = params.vep_cache
        ? channel.value(file(params.vep_cache, checkIfExists: true))
        : PULL_VEP_CACHE().cache
    ascat_resources = PULL_ASCAT_RESOURCES(params.reference_includes_chr_prefix)


    //**************************************************************************************************************************************

    capture_kits = ch_capture_kits

    // Collapse samplesheet rows describing the same library before any work happens, so a
    // normal shared by several tumours is trimmed, aligned, marked, recalibrated, germline
    // called and HLA typed once rather than once per pair. Metas now carry a somatic_names
    // list; fan_out_pairs() restores the scalar somatic_name at each tumour/normal join.
    samplesheet_inputs = dedupe_libraries(ch_samplesheet)

    // Split DNA and RNA sequencing samples
    samples_branched = samplesheet_inputs
    .branch { meta, reads ->
        dna: meta.molecule == "DNA"
        rna: meta.molecule == "RNA"
    }

    dna_inputs = samples_branched.dna
    rna_inputs = samples_branched.rna

    //**************************************************************************************************************************************
    // Perform QC on DNA and RNA samples
    
    dna_fastp = DNA_QC_WORKFLOW(dna_inputs)
    rna_fastp = RNA_QC_WORKFLOW(rna_inputs)

    dna_fastqs = dna_fastp.fastp_fastqs
    rna_fastqs = rna_fastp.fastp_fastqs


    //**************************************************************************************************************************************
    // Split Capture Intervals

    // `as Integer` because a value from the command line arrives as a String, where the
    // same value in -params-file arrives as an Integer. num_intervals reaches
    // groupTuple(size:), which takes only an Integer, so `--scatter_count 20` would end the
    // run after alignment with
    //
    //   Value '20' cannot be used in in parameter 'size' for operator 'groupTuple'
    //
    // The schema types it as an integer, but nf-schema validates params rather than
    // rewriting them, so the coercion has to happen here.
    num_intervals = params.scatter_count as Integer
    intervals = SPLIT_INTERVALS(reference_genome, reference_dict, capture_kits, num_intervals, 0).interval_shards
    processed_regions = CAPTURE_KIT_BED_PROCESS(capture_kits).bed


    //**************************************************************************************************************************************

    // DO WGS/WES ALIGNMENT AND GATK BEST PRACTICES PREPROCESSING
    // ****************************************************************
    
    preproc_bam_workflow = DNA_ALIGN_AND_PREPROC(
        dna_fastqs,
        reference_genome,
        reference_dict,
        bwa_index,
        known_sites_dbsnp,
        known_sites_1000g_snps,
        known_indels,
        mills,
        common_germline,
        intervals,
        num_intervals
    )
    
    preproc_bams = preproc_bam_workflow.preproc_bams // GATK best practices pre-processing (meta, bam, bai)
    markdup_bams = preproc_bam_workflow.markdup_bams // Duplicates marked but no further pre-processing (meta, bam, bai)
    base_recal = preproc_bam_workflow.base_recal // Base recalibration tables from GATK (meta, recal_table)
    markdup_metrics = preproc_bam_workflow.markdup_metrics // samtools markdup duplicate stats (meta, metrics)
    pileup_summaries = preproc_bam_workflow.pileup_summaries // Pileup summaries from GATK (meta, pileup summary)
    flagstats = preproc_bam_workflow.flagstat // samtools flagstat output (meta, flagstat)
    coverage = preproc_bam_workflow.coverage // samtools coverage output (meta, coverage tsv)
    idxstats = preproc_bam_workflow.idxstats // samtools idxstats output (meta, stats tsv)
    
    //**************************************************************************************************************************************
    // HLA Typing
    
    // One HLA call set per library. Merging a pair's tumour and normal BAMs to type them
    // together used to happen here, but mhcflow rejects the result: it requires a single
    // read group and exits on the merge's two with
    //
    //   [helper.py] - [_check_single_rg:132] - [ERROR]: Found more than one read group
    //   information in BAM: [{'ID': 'N1_S1', ...}, {'ID': 'T1', ...}]
    hla_input = markdup_bams
    hla_workflow = HLA_TYPING_WORKFLOW(hla_input, params.reference_includes_chr_prefix,
        hla_reference, hla_kmers, hla_freq)
    
    // Typing output
    optitype = hla_workflow.optitype
    hlahd = hla_workflow.hlahd
    hlahd_tsv = hla_workflow.hlahd_tsv
    hla_pvac_input = hla_workflow.pvac_input

    //**************************************************************************************************************************************
    // bulkRNAseq Processing

    rna = RNASEQ_WORKFLOW(
        rna_fastqs,
        reference_genome,
        gtf,
        params.star_index,
        params.salmon_index,
        transcriptome_reference
    )

    star_bam = rna.star_bam // RNAseq BAM from STAR (meta, bam, bai)
    star_chimeric_out = rna.star_chimeric_out // (meta, chimeric out)
    salmon_tx = rna.salmon_tx // (meta, salmon quant.sf - transcript TPM)
    salmon_gene = rna.salmon_gene // (meta, salmon quant.genes.sf - gene TPM)
    salmon_dir = rna.salmon_dir // (meta, salmon run directory, for MultiQC)
    rna_strandedness = rna.rna_strand // (meta, rna strand prediction)
    
    //**************************************************************************************************************************************
    // Germline Calling on Normal samples
    
    // Branch to get only normal samples
    normal_tumor_split_bams = preproc_bams.branch { meta, bam, bai ->
        tumor: meta.sample_type == "TUMOR"
        normal: meta.sample_type == "NORMAL"
    }

    // Run Haplotype Caller
    // Should make a single germline workflow in the future and put this haplotype caller workflow in it
    germline = GERMLINE_WORKFLOW(
        normal_tumor_split_bams.normal,
        capture_kits,
        processed_regions,
        intervals,
        num_intervals,
        params.interval_padding,
        reference_genome,
        reference_dict,
        hapmap,
        mills,
        vep_cache,
        vep_plugins
    )

    germline_vcf = germline.germline_vcf
    germline_vep_report = germline.germline_vep
    germline_table = germline.germline_vcf_table

    
    //**************************************************************************************************************************************
    // Mutect2 - Somatic Variant Calling

    // (somatic metamap, tumor bam, tumor bai, normal bam, normal bai)
    somatic_samples = pair_tumor_normal(preproc_bams)

    // (somatic metamap, tumor pileups, normal pileups)
    somatic_pileups = pair_tumor_normal(pileup_summaries)
    
    mutect2 = MUTECT2(
        somatic_samples,
        somatic_pileups,
        intervals,
        num_intervals,
        params.interval_padding,
        reference_genome,
        reference_dict,
        gnomad,
        pon
    )
    mutect2_vcf = mutect2.mutect2_vcf
    
    //**************************************************************************************************************************************
    // Strelka - Somatic Variant Calling
    
    // Doesn't use the BQSR bams, just the bams with duplicates marked
    // Get into the format (somatic metadata, tumor bam, tumor bai, normal bam, normal bai)

    markdup_somatic_samples = pair_tumor_normal(markdup_bams)
    
    // Get vcf of strelka called snvs and indels
    strelka = STRELKA_WORKFLOW(markdup_somatic_samples, reference_genome, processed_regions)
    strelka_vcf = strelka.strelka_vcf
    
    
    //**************************************************************************************************************************************
    // Deepsomatic - Somatic Variant Calling

    deepsomatic = DEEPSOMATIC_WORKFLOW(markdup_somatic_samples, reference_genome, capture_kits)
    deepsomatic_vcf = deepsomatic.deepsomatic_vcf
    

    
    //**************************************************************************************************************************************
    // ASCAT CNV ANALYSIS
    
    somatic_samples_bed = markdup_somatic_samples
        .map{meta, tumor_bam, tumor_bai, normal_bam, normal_bai -> 
            tuple(meta.capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai) 
        }
        .combine(capture_kits, by:0)
        .map { kit, meta, tbam, tbai, nbam, nbai, bed ->
                tuple(meta, meta.tumor_meta.sex, nbam, nbai, tbam, tbai, bed)
        }

    
    ascat = ASCAT(somatic_samples_bed, 
        ascat_resources.alleles, 
        ascat_resources.loci, 
        reference_genome,
        ascat_resources.GC, 
        ascat_resources.RT,
    )

    // One publish target for the whole ASCAT result set rather than eight separate ones.
    // HLA LOH. Downstream of both HLA typing and ASCAT: it pairs each tumour with its own
    // normal, realigns it against that normal's HLA reference, and calls loss over the two.
    // Staged as a file rather than referenced through projectDir, so the tiling script
    // travels with the task the way any other input does.
    montage_script = Channel.value(file("${projectDir}/assets/montage_panels.py", checkIfExists: true))

    hla_loh_workflow = HLA_LOH_WORKFLOW(
        markdup_bams,
        hla_workflow.mhcflow_hla_ref,
        hla_workflow.mhcflow_realn_bam,
        ascat.purityploidy,
        hla_workflow.hla_bed,
        hla_kmers,
        hla_freq,
        montage_script
    )

    // ASCAT's plots and its two one-row result files, tiled and merged for the report.
    // join on the meta rather than mixing: the three emits are per pair and have to arrive
    // together, and ASCAT is allowed to fail here, which drops the pair from all three.
    ascat_mqc = ASCAT_MQC(
        ascat.png
            .join(ascat.purityploidy)
            .join(ascat.metrics)
            .map { meta, png, pp, metrics -> tuple(meta.somatic_name, meta, png, pp, metrics) },
        montage_script
    )

    ascat_results = ascat.segments
        .mix(ascat.cnvs, ascat.purityploidy, ascat.metrics, ascat.png,
             ascat.bafs, ascat.logrs, ascat.allelefreqs)




    //**************************************************************************************************************************************
    // Fusion calling with bulkRNAseq using Arriba and STARfusion
    
    fusions = FUSION_CALLING(
        star_bam, 
        star_chimeric_out, 
        reference_genome,
        gtf,
        ctat_bundle,
        arriba_resources
    )

    arriba_fusion = fusions.arriba_fusion
    star_fusion = fusions.star_fusion
    
    
    //**************************************************************************************************************************************
    // PVACtools Input Preparation
   

    // The three callers become one callset here, before anything downstream sees them.
    somatic_consensus = SOMATIC_CONSENSUS(
        mutect2_vcf,
        strelka_vcf,
        deepsomatic_vcf,
        reference_genome,
        reference_dict
    ).vcf

    pvac_input = PVAC_INPUT_PREP_WORKFLOW(
        somatic_consensus,
        markdup_bams,
        star_bam,
        salmon_tx,
        salmon_gene,
        reference_genome,
        reference_dict,
        vep_cache,
        vep_plugins,
        germline_vcf
    )


    somatic_vcf = pvac_input.somatic_vcf
    somatic_vcf_table = pvac_input.somatic_vcf_table
    phased_vcf = pvac_input.phased_vcf
    vep_report = pvac_input.vep_report
    
    //**************************************************************************************************************************************
    // PVACtools Neoantigen Prediction - Somatic Variants / RNA Fusions /

    // pVACtools wants one HLA call set per tumour/normal pair, which used to be the merged
    // sample's. With merging gone it comes from the normal: HLA type is germline, and the
    // normal is free of the tumour's LOH at the HLA locus and of its purity, so it is the
    // better of the two sources rather than only the available one. To predict against the
    // tumour's calls instead, switch the filter below to "TUMOR".
    //
    // fan_out_pairs re-keys by somatic_name, which is what puts a normal shared by two
    // tumours into both pairs.
    combined_hla_somatic_name = fan_out_pairs(
        hla_pvac_input.filter { meta, _calls -> meta.sample_type == "NORMAL" }
    )

    somatic_phased = somatic_vcf.join(phased_vcf).map{meta, somatic_vcf_final, somatic_index, phased_vcf_final, phased_index ->
        tuple(meta.somatic_name, meta, somatic_vcf_final, somatic_index, phased_vcf_final, phased_index)
    }

    pvacseq_input = somatic_phased.join(combined_hla_somatic_name)
    arriba_fusion_somatic_name = fan_out_pairs(arriba_fusion)
    star_fusion_somatic_name = fan_out_pairs(star_fusion)


    pvacfuse_input = arriba_fusion_somatic_name.join(star_fusion_somatic_name).join(combined_hla_somatic_name)

    pvactools = PVACTOOLS_WORKFLOW(pvacseq_input, pvacfuse_input, human_ref_peptides)

    pvacseq = pvactools.pvacseq
    pvacseq_mhc_i_combined = pvactools.pvacseq_mhc_i_combined
    pvacfuse = pvactools.pvacfuse


    somalier = SOMALIER(markdup_bams, reference_genome, somalier_sites)
    somalier_pairs = somalier.pairs
    somalier_samples = somalier.samples


    //**************************************************************************************************************************************
    // MultiQC Report Creation
    mqc_dna_fastp_reports = dna_fastp.fastp_reports.map{ meta, json -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name), meta.sample_name + "_" + meta.molecule ,json) }
    mqc_rna_fastp_reports = rna_fastp.fastp_reports.map{ meta, json -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name),meta.sample_name + "_" + meta.molecule ,json) }
    mqc_base_recal = preproc_bam_workflow.base_recal.map{ meta, table -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name),meta.sample_name + "_" + meta.molecule ,table) } 
    mqc_optitype = optitype.map{ meta, tsv_file, pdf -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name),meta.sample_name + "_" + meta.molecule, tsv_file) }
    mqc_salmon_tx = salmon_dir.map{ meta, quant -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name),meta.sample_name+ "_" + meta.molecule ,quant) }
    mqc_star_log = rna.star_final_log.map{ meta, log -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name),meta.sample_name+ "_" + meta.molecule, log) }
    mqc_vep_report = vep_report.map{meta, html -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name),meta.tumor_meta.sample_name+ "_" + meta.molecule, html) }
    mqc_flagstats = flagstats.map{meta, tsv -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name),meta.sample_name+ "_" + meta.molecule, tsv) }
    mqc_coverage = coverage.map{meta, tsv -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name),meta.sample_name+ "_" + meta.molecule, tsv) }
    mqc_idxstats = idxstats.map{meta, tsv -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name),meta.sample_name+ "_" + meta.molecule, tsv) }
    mqc_germline_vep = germline_vep_report.map{meta, html -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name), meta.sample_name + "_" + meta.molecule, html)}
    mqc_somalier_pairs = somalier_pairs.map{patient, pairs -> tuple(patient, null, null,pairs) }
    mqc_somalier_samples = somalier_samples.map{patient, samples -> tuple(patient, null, null,samples) }
    mqc_hlahd_tsv = hlahd_tsv.map{meta, tsv -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name), meta.sample_name + "_" + meta.molecule, tsv)}
    // New with SAMTOOLS_SORMADUP: MarkDuplicatesSpark was not run with --metrics-file, so
    // the report had no duplicate rate at all. MultiQC's samtools module parses markdup
    // text output.
    mqc_markdup = markdup_metrics.map{meta, metrics -> tuple(meta.patient, (meta.somatic_names ?: meta.somatic_name), meta.sample_name + "_" + meta.molecule, metrics)}
    // Pair-level and one row per image: transpose() because LOHHLA_PLOTS_MQC emits all of a
    // pair's plots as one list, and MultiQC's input is a flat file list. sample_name is null
    // as it is for somalier - these belong to a pair, not a library, and a non-null value
    // here would put a name into --replace-names that matches no sample.
    mqc_loh_plots = hla_loh_workflow.loh_plots_png
        .transpose()
        .map{ _somatic_name, meta, png -> tuple(meta.patient, meta.somatic_name, null, png) }
    mqc_loh_res = hla_loh_workflow.loh_res_mqc
        .map{ _somatic_name, meta, tsv -> tuple(meta.patient, meta.somatic_name, null, tsv) }
    mqc_ascat_plots = ascat_mqc.png
        .map{ _somatic_name, meta, png -> tuple(meta.patient, meta.somatic_name, null, png) }
    mqc_ascat_metrics = ascat_mqc.tsv
        .map{ _somatic_name, meta, tsv -> tuple(meta.patient, meta.somatic_name, null, tsv) }


    mqc_reports = mqc_dna_fastp_reports
        .mix(mqc_rna_fastp_reports)
        .mix(mqc_base_recal)
        .mix(mqc_optitype)
        .mix(mqc_salmon_tx)
        .mix(mqc_star_log)
        .mix(mqc_vep_report)
        .mix(mqc_flagstats)
        .mix(mqc_coverage)
        .mix(mqc_idxstats)
        .mix(mqc_germline_vep)
        .mix(mqc_somalier_pairs)
        .mix(mqc_somalier_samples)
        .mix(mqc_hlahd_tsv)
        .mix(mqc_markdup)
        .mix(mqc_loh_plots)
        .mix(mqc_loh_res)
        .mix(mqc_ascat_plots)
        .mix(mqc_ascat_metrics)
        .groupTuple()
    multiqc = MULTIQC(mqc_reports).html
    
    //**************************************************************************************************************************************

    emit:
    multiqc_reports = multiqc
    markdup_bams = markdup_bams               // duplicates marked, pre-recalibration
    preproc_bams = preproc_bams               // after BQSR
    star_bam = star_bam                       // RNA, coordinate-sorted and indexed
    germline_vcf_table = germline_table       // tabular germline consensus, as for somatic
    ascat_results = ascat_results             // segments, CNVs, purity/ploidy, plots
    hla_loh = hla_loh_workflow.loh_res        // per-pair HLA LOH results from lohhlamod
    hla_loh_plots = hla_loh_workflow.loh_plots // per-gene coverage/logR/BAF profiles
    somatic_vcf = somatic_vcf
    somatic_vcf_table = somatic_vcf_table
    germline_vcf = germline_vcf
    optitype_calls = optitype
    hlahd_calls = hlahd
    hla_pvac_input = hla_pvac_input
    pvacseq = pvacseq
    pvacseq_mhc_i_combined = pvacseq_mhc_i_combined
    pvacfuse = pvacfuse
    salmon_gene = salmon_gene
}
