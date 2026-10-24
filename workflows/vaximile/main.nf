// Helper processes used in the main workflow
include { SPLIT_INTERVALS } from "../../modules/local/split_intervals/main"
include { CAPTURE_KIT_BED_PROCESS } from "../../modules/local/capture_kit_bed_process/main"
include { PULL_VEP_PVAC_PLUGINS } from "../../modules/local/pull_vep_pvac_plugins/main"
include { PULL_VEP_CACHE } from "../../modules/local/pull_vep_cache/main"
include { PULL_CTAT_RESOURCE_BUNDLE } from "../../modules/local/pull_ctat_resource_bundle/main"
include { PULL_ARRIBA_RESOURCES } from "../../modules/local/pull_arriba_resources/main"
include { PULL_ASCAT_RESOURCES } from "../../modules/local/pull_ascat_resources/main"
include { ASCAT } from "../../modules/local/ascat/main"
include { dedupe_libraries; fan_out_pairs; pair_tumor_normal } from "../../subworkflows/pipeline_init.nf"

// Subworkflows
include { DNA_QC_WORKFLOW } from "../../subworkflows/dna_preprocessing.nf"
include { RNA_QC_WORKFLOW } from "../../subworkflows/rnaseq.nf"
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
include { QC } from "../../subworkflows/qc.nf"

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
    // Each comes from its parameter when given, otherwise it is downloaded once into ./vaximile_resources/.
    // All are value channels, so every task that needs one can read it.

    ctat_bundle = params.ctat_resource_lib
        ? channel.value(file(params.ctat_resource_lib, checkIfExists: true))
        : PULL_CTAT_RESOURCE_BUNDLE().ctat_resource_dir
    vep_plugins = params.vep_plugins
        ? channel.value(file(params.vep_plugins, checkIfExists: true))
        : PULL_VEP_PVAC_PLUGINS().plugins
    vep_cache = params.vep_cache
        ? channel.value(file(params.vep_cache, checkIfExists: true))
        : PULL_VEP_CACHE().cache

    // Arriba: any of the three files not given comes from the downloaded release.
    if ( params.arriba_blacklist && params.arriba_known_fusions && params.arriba_protein_domains ) {
        arriba_resources = channel.value(tuple(
            file(params.arriba_blacklist, checkIfExists: true),
            file(params.arriba_known_fusions, checkIfExists: true),
            file(params.arriba_protein_domains, checkIfExists: true)
        ))
    } else {
        arriba_resources = PULL_ARRIBA_RESOURCES().resources.map { blacklist, known_fusions, protein_domains ->
            tuple(
                params.arriba_blacklist ? file(params.arriba_blacklist, checkIfExists: true) : blacklist,
                params.arriba_known_fusions ? file(params.arriba_known_fusions, checkIfExists: true) : known_fusions,
                params.arriba_protein_domains ? file(params.arriba_protein_domains, checkIfExists: true) : protein_domains
            )
        }
    }

    // ASCAT: --ascat_resources holds the four folders the download creates.
    if ( params.ascat_resources ) {
        ascat_alleles = channel.value(file("${params.ascat_resources}/G1000_alleles_hg38", checkIfExists: true))
        ascat_loci = channel.value(file("${params.ascat_resources}/G1000_loci_hg38", checkIfExists: true))
        ascat_gc = channel.value(file("${params.ascat_resources}/GC_G1000_hg38", checkIfExists: true))
        ascat_rt = channel.value(file("${params.ascat_resources}/RT_G1000_hg38", checkIfExists: true))
    } else {
        ascat_download = PULL_ASCAT_RESOURCES(params.reference_includes_chr_prefix)
        ascat_alleles = ascat_download.alleles
        ascat_loci = ascat_download.loci
        ascat_gc = ascat_download.GC
        ascat_rt = ascat_download.RT
    }


    //**************************************************************************************************************************************

    capture_kits = ch_capture_kits

    // Deduplicate libraries first, so a shared normal is processed once.
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

    // CLI values arrive as strings; groupTuple(size:) needs an Integer.
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
    
    //**************************************************************************************************************************************
    // HLA Typing
    
    // HLA typing runs per library; mhcflow requires a single read group per BAM.
    hla_input = markdup_bams
    hla_workflow = HLA_TYPING_WORKFLOW(hla_input, params.reference_includes_chr_prefix)
    
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

    // Strelka and DeepVariant expect unrecalibrated base qualities, so they get the duplicate-marked BAMs.
    normal_markdup_bams = markdup_bams.filter { meta, _bam, _bai -> meta.sample_type == "NORMAL" }

    // Run Haplotype Caller
    // Should make a single germline workflow in the future and put this haplotype caller workflow in it
    germline = GERMLINE_WORKFLOW(
        normal_tumor_split_bams.normal,
        normal_markdup_bams,
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
        ascat_alleles,
        ascat_loci,
        reference_genome,
        ascat_gc,
        ascat_rt,
    )

    // HLA LOH: each tumour against its own normal's HLA reference, using ASCAT's purity and ploidy.
    // mhcflow needs novoalign, which is proprietary and no longer on bioconda, so this runs only
    // when --novoalign_dir points at the user's own copy.
    if ( params.novoalign_dir ) {
        hla_loh_workflow = HLA_LOH_WORKFLOW(
            markdup_bams,
            hla_reference,
            ascat.purityploidy,
            hla_workflow.hla_bed,
            hla_kmers,
            hla_freq,
            channel.value(file(params.novoalign_dir, checkIfExists: true))
        )
        loh_res = hla_loh_workflow.loh_res
        loh_plots = hla_loh_workflow.loh_plots
        loh_plots_png = hla_loh_workflow.loh_plots_png
        loh_res_mqc = hla_loh_workflow.loh_res_mqc
    } else {
        log.warn "--novoalign_dir not given: skipping HLA LOH. See the README for how to supply Novoalign."
        loh_res = channel.empty()
        loh_plots = channel.empty()
        loh_plots_png = channel.empty()
        loh_res_mqc = channel.empty()
    }

    // ASCAT plots and metrics for MultiQC, joined per pair.
    ascat_mqc = ASCAT_MQC(
        ascat.png
            .join(ascat.purityploidy)
            .join(ascat.metrics)
            .map { meta, png, pp, metrics -> tuple(meta.somatic_name, meta, png, pp, metrics) }
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
    somatic_consensus_workflow = SOMATIC_CONSENSUS(
        mutect2_vcf,
        strelka_vcf,
        deepsomatic_vcf,
        reference_genome,
        reference_dict
    )
    somatic_consensus = somatic_consensus_workflow.vcf

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

    // pVACtools uses the normal's HLA calls (HLA type is germline), fanned out to every pair,
    // unless --hla_overrides gives a pair its own alleles.
    hla_overrides = params.hla_overrides
        ? channel.fromPath(params.hla_overrides, checkIfExists: true)
            .splitCsv(header: true, quote: '"')
            .map { row -> tuple(row.somatic_name, row.alleles.replaceAll(/\s/, '')) }
        : channel.empty()

    combined_hla_somatic_name = fan_out_pairs(
        hla_pvac_input.filter { meta, _calls -> meta.sample_type == "NORMAL" }
    )
        .map { somatic_name, meta, calls -> tuple(somatic_name, meta, calls.text.trim()) }
        // remainder: true keeps the pairs without an override. Unknown pairs are rejected at launch.
        .join(hla_overrides, remainder: true)
        .map { somatic_name, meta, typed, override ->
            if ( override ) {
                log.info("${somatic_name}: using HLA alleles from --hla_overrides instead of the typed calls")
            }
            tuple(somatic_name, meta, override ?: typed)
        }

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


    //**************************************************************************************************************************************
    // Quality control and the MultiQC report

    // Reports from earlier steps that MultiQC reads as they are, as (meta, file or files).
    qc_reports = dna_fastp.fastp_reports
        .mix(
            rna_fastp.fastp_reports,
            base_recal,
            markdup_metrics,
            optitype.map { meta, tsv, _pdf -> tuple(meta, tsv) },
            hla_workflow.alleles_mqc,
            hlahd_tsv,
            salmon_dir,
            rna.star_final_log,
            vep_report,
            germline_vep_report,
            ascat_mqc.png.map { _name, meta, png -> tuple(meta, png) },
            ascat_mqc.tsv.map { _name, meta, tsv -> tuple(meta, tsv) },
            loh_plots_png.map { _name, meta, pngs -> tuple(meta, pngs) },
            loh_res_mqc.map { _name, meta, tsv -> tuple(meta, tsv) }
        )

    qc = QC(
        markdup_bams,
        capture_kits,
        reference_genome,
        reference_dict,
        somalier_sites,
        somatic_consensus_workflow.caller_vcfs,
        somatic_consensus,
        somatic_vcf,
        mutect2.contamination,
        pvactools.pvacseq_aggregated,
        germline.caller_vcfs,
        germline_vcf,
        qc_reports
    )
    multiqc = qc.multiqc
    
    //**************************************************************************************************************************************

    emit:
    multiqc_reports = multiqc
    markdup_bams = markdup_bams               // duplicates marked, pre-recalibration
    preproc_bams = preproc_bams               // after BQSR
    star_bam = star_bam                       // RNA, coordinate-sorted and indexed
    germline_vcf_table = germline_table       // tabular germline consensus, as for somatic
    ascat_results = ascat_results             // segments, CNVs, purity/ploidy, plots
    hla_loh = loh_res                         // per-pair HLA LOH results from lohhlamod
    hla_loh_plots = loh_plots                 // per-gene coverage/logR/BAF profiles
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
