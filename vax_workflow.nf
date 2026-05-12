// Default parameter input

params.outdir = "./neoantigen_vax_pipeline_out/"
params.samplesheet =  null

params.kallisto_index = null
params.star_index = null
params.salmon_index = null
params.bwa_index = null
params.vep_cache =  null

params.interval_padding = 100
params.scatter_count = 30

params.reference_fa = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
params.gtf = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
params.transcriptome_reference = "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.human_ref_peptides = "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"

params.intervals_file = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list"
params.common_germline = "gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz"
params.common_germline_index = "gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi"
params.known_sites_dbsnp = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
params.known_sites_dbsnp_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
params.known_sites_1000g_snps = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
params.known_sites_1000g_snps_index = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
params.known_indels = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
params.known_indels_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
params.gnomad = "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
params.gnomad_index = "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi"
params.pon = "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
params.pon_index = "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi"
params.hapmap = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz"
params.hapmap_index = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi"
params.mills = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.mills_index = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"


// Helper processes included in the main workflow
include { SPLIT_INTERVALS; COMBINE_FASTQS; MULTIQC } from "./modules/utilities.nf"
include { PULL_VEP_PVAC_PLUGINS; PULL_CTAT_RESOURCE_BUNDLE; PULL_ARRIBA_RESOURCES } from "./modules/utilities.nf"

// Subworkflows 
include { DNA_QC_WORKFLOW; RNA_QC_WORKFLOW } from "./workflows/qc_workflow.nf"
include { DNA_ALIGN_AND_PREPROC } from "./workflows/align_and_preprocess_workflow.nf"
include { HLA_TYPING_WORKFLOW } from "./workflows/hla_typing_workflow.nf"
include { RNASEQ_WORKFLOW } from "./workflows/rnaseq_workflow.nf"
include { MUTECT2 } from "./workflows/mutect_workflow.nf"
include { STRELKA_WORKFLOW } from "./workflows/strelka_workflow.nf"
include { HAPLOTYPE_CALLER } from "./workflows/germline_workflow.nf"
include { FUSION_CALLING } from  "./workflows/fusion_calling_workflow.nf"
include { PVAC_INPUT_PREP_WORKFLOW } from "./workflows/pvac_vcf_preparation.nf"
include { PVACTOOLS_WORKFLOW } from "./workflows/pvactools_workflow.nf"


// Make a (vcf, vcf index) channel automatically given a vcf path. Tbi must also exist
def make_vcf_channel(vcf_param) {
    return vcf_param
        ? Channel.value(tuple(file(vcf_param), file("${vcf_param}.tbi")))
        : Channel.empty()
}


workflow {
   
    main:
    
    
    // Parse Parameters
    // ****************************************************************
    reference_fa = Channel.fromPath(params.reference_fa)
    reference_fai = Channel.fromPath("${params.reference_fa}.fai")
    reference_dict = Channel.fromPath(params.reference_fa.replace(".fasta",".dict"))
    reference_genome = reference_fa.combine(reference_fai).combine(reference_dict).first()
    common_germline = make_vcf_channel(params.common_germline)
    known_sites_dbsnp = make_vcf_channel(params.known_sites_dbsnp)
    known_sites_1000g_snps = make_vcf_channel(params.known_sites_1000g_snps)
    known_indels = make_vcf_channel(params.known_indels)
    mills = make_vcf_channel(params.mills)
    gnomad = make_vcf_channel(params.gnomad)
    pon = make_vcf_channel(params.pon)
    hapmap = make_vcf_channel(params.hapmap)
    intervals_file = Channel.fromPath(params.intervals_file).first()
    transcriptome_reference = Channel.fromPath(params.transcriptome_reference).first()
    gtf = Channel.fromPath(file(params.gtf)).first()
    human_ref_peptides = Channel.fromPath(file(params.human_ref_peptides)).first()
    bwa_index = params.bwa_index
    vep_plugins = PULL_VEP_PVAC_PLUGINS().plugins
    ctat_bundle = PULL_CTAT_RESOURCE_BUNDLE().ctat_resource_dir
    arriba_resources = PULL_ARRIBA_RESOURCES().resources


    // READ IN SAMPLE DATA FROM SAMPLESHEET
    // ****************************************************************
    samplesheet_inputs = Channel.fromPath(params.samplesheet)
        | splitCsv ( header: true )
            | map { row ->

                def molecule = (row.sequencing_type.toLowerCase() == 'rna') ? 'RNA' :
                    (row.sequencing_type.toLowerCase() in ['exome', 'genome', "exome_FFPE", "genome_FFPE"]) ? 'DNA' :
                    null

                def meta = [
                    somatic_name: row.somatic_name,
                    sample_name: row.sample_name,
                    sample_type: row.sample_type,
                    sequencing_type: row.sequencing_type,
                    molecule: molecule
                ]

                def reads = [
                    file(row.fastqr1, checkIfExists: true),
                    file(row.fastqr2, checkIfExists: true)

                ]

                return [meta, reads]
            }

    // Split DNA and RNA sequencing samples
    samples_branched = samplesheet_inputs.branch { meta, reads ->
        dna: meta.molecule == "DNA"
        rna: meta.molecule == "RNA"
        }
    dna_inputs = samples_branched.dna
    rna_inputs = samples_branched.rna

    
    // RUN FASTP QC ON ALL SAMPLES
    // ****************************************************************
    
    dna_fastp = DNA_QC_WORKFLOW(dna_inputs)
    rna_fastp = RNA_QC_WORKFLOW(rna_inputs)

        
    // SPLIT CAPTURE INTERVALS
    // ****************************************************************

    num_intervals = params.scatter_count
    intervals = SPLIT_INTERVALS(reference_genome, intervals_file, num_intervals, 0)
        .flatten()
        .map { file -> tuple(file.baseName, file) }
    
    // DO WGS/WES ALIGNMENT AND GATK BEST PRACTICES PREPROCESSING
    // ****************************************************************
    
    
    // Combine dna fastqs from the same somatic sample ie Tumor + Normal
    // Currently only works with two samples
    sample_grouped_fastqs = dna_fastp.fastp_fastqs
        .map {meta, fastq1, fastq2 ->
            def new_meta = [
                somatic_name: meta.somatic_name,
                sample_name: meta.somatic_name,
                sample_type: "Tumor_Normal",
                sequencing_type: meta.sequencing_type,
                molecule: meta.molecule
            ]
            tuple(new_meta, fastq1, fastq2)
        }
        .groupTuple()
        .map { meta, fastqs_r1, fastqs_r2 ->
            tuple(meta, fastqs_r1.sort(), fastqs_r2.sort())
        }
    combined_fastqs = COMBINE_FASTQS(sample_grouped_fastqs)


    all_dna_fastqs = dna_fastp.fastp_fastqs.mix(combined_fastqs)

    preproc_bam_workflow = DNA_ALIGN_AND_PREPROC(
        all_dna_fastqs,
        reference_genome,
        bwa_index,
        known_sites_dbsnp,
        known_sites_1000g_snps,
        known_indels,
        mills,
        common_germline,
        intervals,
        num_intervals
    )
    
    
    // GATK best practices preprocess bams ie BQSR (meta, bam, bai)
    preproc_bams_all = preproc_bam_workflow.preproc_bams
    preproc_bams_samples = preproc_bams_all.branch{meta, bam, bai ->
        merged: meta.sample_type == "Tumor_Normal"
        single: meta.sample_type == "Tumor" || meta.sample_type == "Normal"
    }
    preproc_bams = preproc_bams_samples.single

    // Bam with duplicates marked but no further processing (meta, bam, bai)
    markdup_bams_all = preproc_bam_workflow.markdup_bams
    markdup_bams_samples = markdup_bams_all.branch{meta, bam, bai ->
        merged: meta.sample_type == "Tumor_Normal"
        single: meta.sample_type == "Tumor" || meta.sample_type == "Normal"
    }
    markdup_bams = markdup_bams_samples.single

    
    // Base recalibration tables from GATK (meta, recal_table)
    base_recal_all = preproc_bam_workflow.base_recal
    base_recal_samples = base_recal_all.branch{meta, recal_table ->
        merged: meta.sample_type == "Tumor_Normal"
        single: meta.sample_type == "Tumor" || meta.sample_type == "Normal"
    }
    base_recal = base_recal_samples.single


    // Pileup summaries from GATK (meta, pileup summary)
    pileup_summaries_all = preproc_bam_workflow.pileup_summaries     
    pileup_summaries_samples = pileup_summaries_all.branch{meta, summary ->
        merged: meta.sample_type == "Tumor_Normal"
        single: meta.sample_type == "Tumor" || meta.sample_type == "Normal"
    }
    pileup_summaries = pileup_summaries_samples.single 

    // HLA TYPING
    // ****************************************************************
    
    hla_workflow = HLA_TYPING_WORKFLOW(markdup_bams_all)
    
    optitype = hla_workflow.optitype
    hlahd = hla_workflow.hlahd
    hla_pvac_input = hla_workflow.pvac_input

    
    // RNASEQ PROCESSING
    // ****************************************************************

    
    rna = RNASEQ_WORKFLOW(
        rna_fastp.fastp_fastqs,
        reference_genome,
        gtf,
        params.star_index,
        params.kallisto_index,
        params.salmon_index,
        transcriptome_reference
    )

    star_bam = rna.star_bam // RNAseq BAM from STAR (meta, bam, bai)
    kallisto_tx = rna.kallisto_tx // (meta, kallisto abundance)
    kallisto_gene = rna.kallisto_gene // (meta, kallisto gene abundance)
    salmon_tx = rna.salmon_tx // (meta, salmon quant abundance)
    rna_strandedness = rna.rna_strand // (meta, rna strand prediction)

    
    
    // GERMLINE CALLING
    // ****************************************************************
    
    // Run germline calling only on Normal samples
    bams_type_split = preproc_bams.branch { meta, bam, bai ->
        tumor: meta.sample_type == "Tumor"
        normal: meta.sample_type == "Normal"
    }

    // Run basic haplotype calling on individual samples
    haplotype_caller = HAPLOTYPE_CALLER(
        bams_type_split.normal,
        intervals,
        num_intervals,
        params.interval_padding,
        reference_genome,
        hapmap,
        mills
    )

    germline_vcf = haplotype_caller.germline_vcf
        
    
    
    // MUTECT2 VARIANT CALLING
    // ****************************************************************

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
            tumor_meta: tumor_meta,
            normal_meta: normal_meta
        ]
        tuple(somatic_meta, tumor_pileup, normal_pileup)
    }
    
    // Run mutect2 on paired Tumor-Normal samples
    mutect = MUTECT2(
        somatic_samples,
        somatic_pileups,
        intervals,
        num_intervals,
        params.interval_padding,
        reference_genome,
        gnomad,
        pon
    )

    // ****************************************************************
    // STRELKA CALLING
    
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
            tumor_meta: tumor_meta,
            normal_meta: normal_meta
        ]
        tuple(somatic_meta, tumor_bam, tumor_bai, normal_bam, normal_bai)
    }
    
    // Get vcf of strelka called snvs and indels
    strelka = STRELKA_WORKFLOW(markdup_somatic_samples, reference_genome, intervals_file)

    
    // ****************************************************************
    // FUSION CALLING
   
    /*
    fusions = FUSION_CALLING(
        rna.star_bam, 
        rna.star_chimeric_out, 
        reference_genome,
        gtf,
        ctat_bundle,
        arriba_resources
    )
    */
    

    // ****************************************************************
    // PVACTOOLS INPUT FILE PREPARATION

    pvac_input = PVAC_INPUT_PREP_WORKFLOW(
        mutect.mutect2_vcf,
        strelka.strelka_vcf,
        markdup_bams,
        rna.star_bam,
        rna.kallisto_tx,
        rna.kallisto_gene,
        reference_genome,
        params.vep_cache,
        vep_plugins,
        germline_vcf
    )

    somatic_vcf = pvac_input.somatic_vcf
    phased_vcf = pvac_input.phased_vcf
    vep_report = pvac_input.vep_report

    
    /*

    // ****************************************************************
    // PVACTOOLS NEOANTIGEN PREDICTIONS

    hla_branched = hla_pvac_input.branch{meta, calls ->
        combined: meta.sample_type == "Tumor_Normal"
        tumor: meta.sample_type == "Tumor"
        normal: meta.sample_type == "Normal"
    }

    combined_hla = hla_branched.combined.map{meta, calls ->
        tuple(meta.somatic_name, meta, calls)
    }

    somatic_phased = somatic_vcf.join(phased_vcf).map{meta, somatic_vcf, somatic_index, phased_vcf, phased_index ->
        tuple(meta.somatic_name, meta, somatic_vcf, somatic_index, phased_vcf, phased_index)
    }

    pvacseq_input = somatic_phased.join(combined_hla)

    //pvactools = PVACTOOLS_WORKFLOW(pvacseq_input, human_ref_peptides)

    
    mqc_dna_fastp_reports = dna_fastp.fastp_reports.map{ meta, json -> tuple(meta.somatic_name, json) } // sample_meta
    mqc_rna_fastp_reports = rna_fastp.fastp_reports.map{ meta, json -> tuple(meta.somatic_name, json) }  //sample_meta
    mqc_base_recal = preproc_bam_workflow.base_recal.map{ meta, table -> tuple(meta.somatic_name, table) }  // Base recalibration tables from GATK (meta, recal_table)
    mqc_optitype = optitype.map{ meta, tsv_file, pdf -> tuple(meta.somatic_name, tsv_file) }
    mqc_salmon_tx = salmon_tx.map{ meta, quant -> tuple(meta.somatic_name, quant) }
    mqc_star_log = rna.star_final_log.map{ meta, log -> tuple(meta.somatic_name, log) }
    mqc_vep_report = vep_report.map{meta, html -> tuple(meta.somatic_name, html) }

    mqc_reports = mqc_dna_fastp_reports
        .mix(mqc_rna_fastp_reports)
        .mix(mqc_base_recal)
        .mix(mqc_optitype)
        .mix(mqc_salmon_tx)
        .mix(mqc_star_log)
        .mix(mqc_vep_report)
        .groupTuple()



    multiqc = MULTIQC(mqc_reports)

    */
    
    publish:
        dna_fastp_fastqs = dna_fastp.fastp_fastqs
}

output {
    dna_fastp_fastqs{
        path { meta, fastq1, fastq2 -> "${params.outdir}/${meta.somatic_name}/qc/fastp/${meta.sample_name}" }
    }

}
