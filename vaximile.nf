
// REQUIRED PARAMETERS
//**************************************************************************************************************************************
// Define default input and output files
params.outdir = "./neoantigen_vax_pipeline_out/"
params.samplesheet =  null
params.capture_kits = null
params.vep_cache =  null

//**************************************************************************************************************************************
// Optionally supply indices needed for alignment / quantification, otherwise they will be created.
params.kallisto_index = null
params.star_index = null
params.salmon_index = null
params.bwa_index = null

//**************************************************************************************************************************************
// Default reference genome, GTF, transcriptome for rnaseq, human reference proteome peptides for pvactools
params.reference_fa = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"
params.gtf = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
params.transcriptome_reference = "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.human_ref_peptides = "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"

//**************************************************************************************************************************************
// Default intervals file uses WGS regions from GATK to call over entire genome.
params.intervals_file = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list"
params.interval_padding = 100
params.scatter_count = 30

//**************************************************************************************************************************************
// GATK RESOURCE VCFS USED FOR MUTECT2, HAPLOTYPECALLER, DATA PRE-PROCESSING
params.common_germline = "gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz"
params.known_sites_dbsnp = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
params.known_sites_1000g_snps = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
params.known_indels = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
params.gnomad = "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
params.pon = "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
params.hapmap = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz"
params.mills = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.somalier_sites = "https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz"

//**************************************************************************************************************************************
// Helper processes used in the main workflow
include { SPLIT_INTERVALS; COMBINE_FASTQS; MERGE_BAMS; CAPTURE_KIT_BED_PROCESS } from "./modules/utilities.nf"
include { PULL_VEP_PVAC_PLUGINS; PULL_CTAT_RESOURCE_BUNDLE; PULL_ARRIBA_RESOURCES } from "./modules/utilities.nf"
include { MULTIQC } from "./modules/quality_control.nf"

// Subworkflows
include { DNA_QC_WORKFLOW; RNA_QC_WORKFLOW; SOMALIER } from "./workflows/qc_workflow.nf"
include { DNA_ALIGN_AND_PREPROC; BWA_INDEX; PREPARE_REFERENCE_FASTA } from "./workflows/align_and_preprocess_workflow.nf"
include { HLA_TYPING_WORKFLOW } from "./workflows/hla_typing_workflow.nf"
include { RNASEQ_WORKFLOW } from "./workflows/rnaseq_workflow.nf"
include { MUTECT2 } from "./workflows/mutect_workflow.nf"
include { STRELKA_WORKFLOW } from "./workflows/strelka_workflow.nf"
include { GERMLINE_WORKFLOW } from "./workflows/germline_workflow.nf"
include { FUSION_CALLING } from  "./workflows/fusion_calling_workflow.nf"
include { PVAC_INPUT_PREP_WORKFLOW } from "./workflows/pvac_vcf_preparation.nf"
include { PVACTOOLS_WORKFLOW } from "./workflows/pvactools_workflow.nf"
include { DEEPSOMATIC_WORKFLOW } from "./workflows/deepsomatic.nf"

//**************************************************************************************************************************************

// Make a (vcf, vcf index) channel automatically given a vcf path. Tbi must also exist
def make_vcf_channel(vcf_param) {
    return vcf_param
        ? Channel.value(tuple(file(vcf_param), file("${vcf_param}.tbi")))
        : Channel.empty()
}


// Main Workflow
workflow {
   
    main:

    //**************************************************************************************************************************************

    reference_fa = Channel.fromPath(params.reference_fa)
    prepared_reference = PREPARE_REFERENCE_FASTA(reference_fa)
    reference_genome = prepared_reference.fa_fai_pair.first()
    reference_dict = prepared_reference.dict.first()

    bwa_index_input = params.bwa_index
    if (bwa_index_input  == null) {
        bwa_index = BWA_INDEX(reference_genome, bwa_index_input)
    } else {
        bwa_index = Channel.fromPath("${params.bwa_index}/*.{amb,ann,bwt,pac,sa}", checkIfExists: true).collect()
    }

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

    //**************************************************************************************************************************************
    // Pull Resources

    ctat_bundle = PULL_CTAT_RESOURCE_BUNDLE().ctat_resource_dir
    arriba_resources = PULL_ARRIBA_RESOURCES().resources
    vep_plugins = PULL_VEP_PVAC_PLUGINS().plugins


    //**************************************************************************************************************************************
    // Read in Capture Kits

    capture_kits = Channel.fromPath(params.capture_kits)
        .splitCsv(header: true)
        .map { row -> tuple(row.kit, row.bed) }


    //**************************************************************************************************************************************
    // Read samplesheet into channels of form (sample metamap, (fastq1, fastq2))
    
    samplesheet_inputs = Channel.fromPath(params.samplesheet)
        | splitCsv ( header: true )
            | map { row ->

                // Be semi-flexible in parsing the sequencing type
                def molecule = (row.sequencing_type.toLowerCase() == 'rna') ? 'RNA' :
                    (row.sequencing_type.toLowerCase() in ['exome', 'genome', "exome_FFPE", "genome_FFPE"]) ? 'DNA' :
                    null

                def meta = [
                    somatic_name: row.somatic_name,
                    patient: row.patient,
                    sample_name: row.sample_name,
                    sample_type: row.sample_type,
                    sequencing_type: row.sequencing_type,
                    capture_kit: row.capture_kit,
                    molecule: molecule
                ]

                def reads = [
                    file(row.fastqr1, checkIfExists: true),
                    file(row.fastqr2, checkIfExists: true)

                ]

                return [meta, reads]
            }
    
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

    num_intervals = params.scatter_count
    intervals = SPLIT_INTERVALS(reference_genome, reference_dict, capture_kits, num_intervals, 0)
    processed_regions = CAPTURE_KIT_BED_PROCESS(capture_kits)


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
    pileup_summaries = preproc_bam_workflow.pileup_summaries // Pileup summaries from GATK (meta, pileup summary)
    flagstats = preproc_bam_workflow.flagstat // samtools flagstat output (meta, flagstat)
    coverage = preproc_bam_workflow.coverage // samtools coverage output (meta, coverage tsv)
    idxstats = preproc_bam_workflow.idxstats // samtools idxstats output (meta, stats tsv)
    
    //**************************************************************************************************************************************
    // HLA Typing
    
    // Group bams by meta.somatic_name to be merged
    bams_merged_input = markdup_bams.map{meta, bam, bai ->
        tuple(meta.somatic_name, meta.patient, meta, bam, bai)
    }
    .groupTuple(by: [0,1])
    .map{somatic_name, patient, metas, bams, bais ->
            def new_meta = [
                somatic_name: somatic_name,
                patient: patient,
                sample_name: "Merged_" + somatic_name,
                sample_type: "Merged",
                molecule: "DNA"
            ]

            tuple(new_meta, bams, bais)
    }
    
    combined_bams = MERGE_BAMS(bams_merged_input) // Merge BAMs sharing somatic_name

    // Call HLA alleles on individual samples AND merged samples
    hla_input = markdup_bams.mix(combined_bams)
    hla_workflow = HLA_TYPING_WORKFLOW(hla_input)
    
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
        params.kallisto_index,
        params.salmon_index,
        transcriptome_reference
    )

    star_bam = rna.star_bam // RNAseq BAM from STAR (meta, bam, bai)
    star_chimeric_out = rna.star_chimeric_out // (meta, chimeric out)
    kallisto_tx = rna.kallisto_tx // (meta, kallisto abundance)
    kallisto_gene = rna.kallisto_gene // (meta, kallisto gene abundance)
    salmon_tx = rna.salmon_tx // (meta, salmon quant abundance)
    rna_strandedness = rna.rna_strand // (meta, rna strand prediction)
    
    //**************************************************************************************************************************************
    // Germline Calling on Normal samples
    
    // Branch to get only normal samples
    normal_tumor_split_bams = preproc_bams.branch { meta, bam, bai ->
        tumor: meta.sample_type == "Tumor"
        normal: meta.sample_type == "Normal"
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
        params.vep_cache,
        vep_plugins
    )

    germline_vcf = germline.germline_vcf
    germline_vep_report = germline.germline_vep
    germline_table = germline.germline_vcf_table

    
    //**************************************************************************************************************************************
    // Mutect2 - Somatic Variant Calling

     
    // Groups samples into the form (somatic metamap, tumor bam, tumor bai, normal bam, normal bai)
    preproc_bams_type_branched = preproc_bams.map{meta, bam, bai ->
        tuple(meta.somatic_name, meta, bam, bai)
    }
    .branch {somatic_name, meta, bam, bai ->
        tumor: meta.sample_type == "Tumor"
        normal: meta.sample_type == "Normal"
    }

    paired_samples = preproc_bams_type_branched.tumor.join(preproc_bams_type_branched.normal)
    somatic_samples = paired_samples.map {somatic_name, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
        def somatic_meta = [
            somatic_name: somatic_name,
            patient: tumor_meta.patient,
            capture_kit: tumor_meta.capture_kit,
            tumor_meta: tumor_meta,
            normal_meta: normal_meta
        ]
        tuple(somatic_meta, tumor_bam, tumor_bai, normal_bam, normal_bai)
    }

    // Group pileups into the form (somatic metamap, tumor pileups, normal pileups)
    type_pileups = pileup_summaries.map{meta,pileup_table ->
        tuple(meta.somatic_name, meta, pileup_table)
    }
    .branch {somatic_name, meta, pileup_table ->
        tumor: meta.sample_type == "Tumor"
        normal: meta.sample_type == "Normal"
    }
    paired_pileups = type_pileups.tumor.join(type_pileups.normal)


    somatic_pileups = paired_pileups.map {somatic_name, tumor_meta, tumor_pileup, normal_meta, normal_pileup ->
        def somatic_meta = [
            somatic_name: somatic_name,
            patient: tumor_meta.patient,
            capture_kit: tumor_meta.capture_kit,
            tumor_meta: tumor_meta,
            normal_meta: normal_meta
        ]
        tuple(somatic_meta, tumor_pileup, normal_pileup)
    }
    
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

    markdup_bams_type_branched = markdup_bams.map{meta, bam, bai ->
        tuple(meta.somatic_name, meta, bam, bai)
    }
    .branch {somatic_name, meta, bam, bai ->
        tumor: meta.sample_type == "Tumor"
        normal: meta.sample_type == "Normal"
    }

    markdup_paired_samples = markdup_bams_type_branched.tumor.join(markdup_bams_type_branched.normal)

    markdup_somatic_samples = markdup_paired_samples.map {somatic_name, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
        def somatic_meta = [
            somatic_name: somatic_name,
            patient: tumor_meta.patient,
            capture_kit: tumor_meta.capture_kit,
            tumor_meta: tumor_meta,
            normal_meta: normal_meta
        ]
        tuple(somatic_meta, tumor_bam, tumor_bai, normal_bam, normal_bai)
    }

    // Get vcf of strelka called snvs and indels
    strelka = STRELKA_WORKFLOW(markdup_somatic_samples, reference_genome, processed_regions)
    strelka_vcf = strelka.strelka_vcf
    
    
    //**************************************************************************************************************************************
    // Deepsomatic - Somatic Variant Calling

    deepsomatic = DEEPSOMATIC_WORKFLOW(markdup_somatic_samples, reference_genome, capture_kits)
    deepsomatic_vcf = deepsomatic.deepsomatic_vcf

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

    pvac_input = PVAC_INPUT_PREP_WORKFLOW(
        mutect2_vcf,
        strelka_vcf,
        deepsomatic_vcf,
        markdup_bams,
        star_bam,
        kallisto_tx,
        kallisto_gene,
        reference_genome,
        reference_dict,
        params.vep_cache,
        vep_plugins,
        germline_vcf
    )

    somatic_vcf = pvac_input.somatic_vcf
    somatic_vcf_table = pvac_input.somatic_vcf_table
    phased_vcf = pvac_input.phased_vcf
    vep_report = pvac_input.vep_report
    
    //**************************************************************************************************************************************
    // PVACtools Neoantigen Prediction - Somatic Variants / RNA Fusions /

    hla_branched = hla_pvac_input.branch{meta, calls ->
        merged: meta.sample_type == "Merged"
        tumor: meta.sample_type == "Tumor"
        normal: meta.sample_type == "Normal"
    }

    combined_hla_somatic_name = hla_branched.merged.map{meta, calls ->
        tuple(meta.somatic_name, meta, calls)
    }

    somatic_phased = somatic_vcf.join(phased_vcf).map{meta, somatic_vcf_final, somatic_index, phased_vcf_final, phased_index ->
        tuple(meta.somatic_name, meta, somatic_vcf_final, somatic_index, phased_vcf_final, phased_index)
    }

    pvacseq_input = somatic_phased.join(combined_hla_somatic_name)
    arriba_fusion_somatic_name = arriba_fusion.map{meta, tsv -> tuple(meta.somatic_name, meta, tsv) }
    star_fusion_somatic_name = star_fusion.map{meta, tsv -> tuple(meta.somatic_name, meta, tsv) }


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
    mqc_dna_fastp_reports = dna_fastp.fastp_reports.map{ meta, json -> tuple(meta.patient, meta.somatic_name, meta.sample_name + "_" + meta.molecule ,json) }.unique{meta, som, sample, x -> sample}
    mqc_rna_fastp_reports = rna_fastp.fastp_reports.map{ meta, json -> tuple(meta.patient, meta.somatic_name,meta.sample_name + "_" + meta.molecule ,json) }.unique{meta, som,sample, x -> sample}
    mqc_base_recal = preproc_bam_workflow.base_recal.map{ meta, table -> tuple(meta.patient, meta.somatic_name,meta.sample_name + "_" + meta.molecule ,table) } .unique{meta, som,sample, x -> sample} 
    mqc_optitype = optitype.map{ meta, tsv_file, pdf -> tuple(meta.patient, meta.somatic_name,meta.sample_name + "_" + meta.molecule, tsv_file) }.unique{meta, som,sample, x -> sample}
    mqc_salmon_tx = salmon_tx.map{ meta, quant -> tuple(meta.patient, meta.somatic_name,meta.sample_name+ "_" + meta.molecule ,quant) }.unique{meta, som,sample, x -> sample}
    mqc_star_log = rna.star_final_log.map{ meta, log -> tuple(meta.patient, meta.somatic_name,meta.sample_name+ "_" + meta.molecule, log) }.unique{meta, som,sample, x -> sample}
    mqc_vep_report = vep_report.map{meta, html -> tuple(meta.patient, meta.somatic_name,meta.tumor_meta.sample_name+ "_" + meta.molecule, html) }.unique{meta, som,sample, x -> sample}
    mqc_flagstats = flagstats.map{meta, tsv -> tuple(meta.patient, meta.somatic_name,meta.sample_name+ "_" + meta.molecule, tsv) }.unique{meta, som,sample, x -> sample}
    mqc_coverage = coverage.map{meta, tsv -> tuple(meta.patient, meta.somatic_name,meta.sample_name+ "_" + meta.molecule, tsv) }.unique{meta, som,sample, x -> sample}
    mqc_idxstats = idxstats.map{meta, tsv -> tuple(meta.patient, meta.somatic_name,meta.sample_name+ "_" + meta.molecule, tsv) }.unique{meta, som,sample, x -> sample}
    mqc_germline_vep = germline_vep_report.map{meta, html -> tuple(meta.patient, meta.somatic_name, meta.sample_name + "_" + meta.molecule, html)}.unique{meta, som,sample, x -> sample}
    mqc_somalier_pairs = somalier_pairs.map{patient, pairs -> tuple(patient, null, null,pairs) }
    mqc_somalier_samples = somalier_samples.map{patient, samples -> tuple(patient, null, null,samples) }
    mqc_hlahd_tsv = hlahd_tsv.map{meta, tsv -> tuple(meta.patient, meta.somatic_name, meta.sample_name + "_" + meta.molecule, tsv)}.unique{meta, som,sample, x -> sample}


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
        .groupTuple()
    multiqc = MULTIQC(mqc_reports)
    
    //**************************************************************************************************************************************

    publish:
        multiqc_reports = multiqc
        somatic_vcf = somatic_vcf
        somatic_vcf_table = somatic_vcf_table
        germline_vcf = germline_vcf
        optitype_calls = optitype 
        hlahd_calls =  hlahd
        hla_pvac_input = hla_pvac_input
        pvacseq = pvacseq
        pvacseq_mhc_i_combined = pvacseq_mhc_i_combined
        pvacfuse = pvacfuse
        kallisto_gene = kallisto_gene
}

output {
    multiqc_reports {
        path { patient, report -> "${params.outdir}/${patient}/multiqc/" }
    }
    somatic_vcf {
        path { meta, vcf, vcf_index -> "${params.outdir}/${meta.patient}/${meta.somatic_name}/variants" }
    }
    somatic_vcf_table {
        path { meta, table -> "${params.outdir}/${meta.patient}/${meta.somatic_name}/variants" }
    }
    optitype_calls {
        path { meta, tsv, pdf -> "${params.outdir}/${meta.patient}/${meta.somatic_name}/hla/optitype/" }
    }
    hlahd_calls {
        path { meta, calls -> "${params.outdir}/${meta.patient}/${meta.somatic_name}/hla/hlahd/" }
    }
    hla_pvac_input {
        path { meta, calls -> "${params.outdir}/${meta.patient}/${meta.somatic_name}/hla/" }
    }
    pvacseq {
        path { meta, pvacseq_dir -> "${params.outdir}/${meta.patient}/${meta.somatic_name}/pvactools/" }
    }
    pvacseq_mhc_i_combined {
        path { patient, report -> "${params.outdir}/${patient}/pvactools_report" }
    }
    pvacfuse {
        path { meta, pvacfuse_dir -> "${params.outdir}/${meta.patient}/${meta.somatic_name}/pvactools" }
    }
    
    germline_vcf {
        path { meta, vcf, tbi -> "${params.outdir}/${meta.patient}/${meta.somatic_name}/germline/" }
    }
    kallisto_gene {
        path { meta, gene_abundance -> "${params.outdir}/${meta.patient}/${meta.somatic_name}/kallisto" }
    }

}
