// Default parameter input

params.outdir = "./neoantigen_vax_pipeline_out/"
params.dna_sample_sheet =  null
params.rna_sample_sheet =  null

params.kallisto_index = null
params.star_index = null
params.bwa_index = null
params.vep_cache =  null
params.vep_plugins = null
params.ctat_resource_dir = null
params.arriba_blacklist = null 
params.arriba_known_fusions = null
params.arriba_protein_domains = null

params.reference_fa = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
params.reference_index_dir = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38{.fasta.fai,.dict}"
params.reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
params.gencode_gtf = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
params.kallisto_reference =  "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.human_ref_peptides = "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"

params.scatter_count = 30
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


include { SPLIT_INTERVALS; SPLIT_INTERVALS_PADDED; COMBINE_FASTQS } from "./modules/utilities.nf"
//include { FASTP } from "./modules/qc.nf"
//include { CREATE_STAR_INDEX; STAR_ALIGN; STAR_INDEX_BAM; CREATE_BWA_INDEX; BWA_MAP; PREPROCESS_BAM } from "./modules/alignment.nf"
//include { STAR_FUSION; ARRIBA_FUSION } from "./modules/fusion_calling.nf"
//include { VT_POSTPROCESS_GERMLINE; VEP_ANNOTATE_GERMLINE; INDEX_FINAL_VCF_GERMLINE; POSTPROCESS_HAPLOTYPE_SCATTER; HAPLOTYPE_CALLER_SCATTER } from "./modules/germline.nf"
//include { HLAHD_HLA_CALLS; OPTITYPE_HLA_CALLS; POSTPROCESS_OPTITYPE; OPTITYPE; HLAHD; HLA_COMBINE_FASTQS } from "./modules/hla_typing.nf"
//include { PVACSEQ; PVACFUSE } from "./modules/pvactools.nf"
//include { BAMREADCOUNT; VEP_ANNOTATE; ANNOTATE_VCF_COVERAGE; ANNOTATE_VCF_EXPRESSION; PHASE_VCF_SELECT_VARIANTS; 
//    PHASE_VCF_COMBINE_VARIANTS; PHASE_VCF_SORT_VCF; PHASE_VCF_RENAME; PHASE_VCF_RBPHASING; PHASE_VCF_VEP; PHASE_VCF_INDEX } from "./modules/pvactools_vcf_prep.nf"
//include { KALLISTO_INDEX; KALLISTO_QUANT; KALLISTO_TXIMPORT } from "./modules/rnaseq.nf"
//include { DEEPSOMATIC; STRELKA; POSTPROCESS_STRELKA; MUTECT2_SCATTER; POSTPROCESS_MUTECT2_SCATTER; SPLIT_INTERVALS } from "./modules/somatic_calling.nf"
//include { INDEX_FINAL_VCF; VCF_TO_TABLE; MERGE_SOMATIC_VCFS; FILTER_VCF; ADD_VCF_GT_FIELD; VEP_FILTER; VT_SOMATIC_POSTPROCESS } from "./modules/somatic_postprocess.nf"
    
include { DNA_QC_WORKFLOW; RNA_QC_WORKFLOW } from "./workflows/qc_workflow.nf"
include { DNA_ALIGNMENT_WORKFLOW } from "./workflows/alignment_workflow.nf"
include { HLA_TYPING_WORKFLOW } from "./workflows/hla_typing_workflow.nf"

def make_vcf_channel(vcf_param) {
    return vcf_param
        ? Channel.value(tuple(file(vcf_param), file("${vcf_param}.tbi")))
        : Channel.empty()
}

workflow {
   
    main:
    // INPUT PARSING
    // READ IN SAMPLE DATA FROM SAMPLESHEET

    dna_inputs = Channel.fromPath(params.dna_sample_sheet)
        | splitCsv( header: true )
            | map { row ->

                def meta = [
                    somatic_name: row.somatic_name,
                    sample_name: row.sample_name,
                    sample_type: row.sample_type,
                    sequencing_type: row.sequencing_type,
                    molecule: "DNA"
                ]

                def reads = [
                    file(row.fastqr1, checkIfExists:true),
                    file(row.fastqr2, checkIfExists:true)
                ]

                return [meta, reads]
            }


    rna_inputs = Channel.fromPath(params.rna_sample_sheet)
        | splitCsv ( header: true )
            | map { row ->
                def meta = [
                    somatic_name: row.somatic_name,
                    sample_name: row.sample_name,
                    sample_type: row.sample_type,
                    strand: row.strand,
                    molecule: "RNA"
                ]
                
                def reads = [
                    file(row.fastqr1, checkIfExists: true),
                    file(row.fastqr2, checkIfExists: true)
                ]

                return [meta, reads]

            }


    
    // OLD BLOCK
    /*
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
    */
    
    
    // ****************************************************************
    // PULL REFERENCE FASTA AND REFERENCE FASTA SUPPLEMENTAL FILES
    
    //reference_fa = Channel.fromPath(params.reference_fa).first()
    //reference_index_files = Channel.fromPath(params.reference_index_dir).collect()
    //reference_dict = Channel.fromPath(params.reference_dict).first()
    
    reference_fa = Channel.fromPath(params.reference_fa)
    reference_fai = Channel.fromPath("${params.reference_fa}.fai")
    reference_dict = Channel.fromPath(params.reference_fa.replace(".fasta",".dict"))

    // Combine them into a tuple channel if needed together
    reference_genome = reference_fa.combine(reference_fai).combine(reference_dict).first()


    common_germline = make_vcf_channel(params.common_germline)
    known_sites_dbsnp = make_vcf_channel(params.known_sites_dbsnp)
    known_sites_1000g_snps = make_vcf_channel(params.known_sites_1000g_snps)
    known_indels = make_vcf_channel(params.known_indels)
    mills = make_vcf_channel(params.mills)
    gnomad = make_vcf_channel(params.gnomad)
    pon = make_vcf_channel(params.pon)
    hapmap = make_vcf_channel(params.hapmap)

    
    /*
    //common_germline = Channel.fromPath(params.common_germline).first()
    //common_germline_index = Channel.fromPath(params.common_germline_index).first()
    
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
    */
    
    intervals_file = Channel.fromPath(params.intervals_file).first()
    
    kallisto_reference = Channel.fromPath(params.kallisto_reference).first()
    gencode_gtf = Channel.fromPath(file(params.gencode_gtf)).first()
    
    ctat_resource_dir = Channel.fromPath(file(params.ctat_resource_dir)).first()
   
    arriba_blacklist = Channel.fromPath(file(params.arriba_blacklist)).first()
    arriba_known_fusions = Channel.fromPath(file(params.arriba_known_fusions)).first()
    arriba_protein_domains = Channel.fromPath(file(params.arriba_protein_domains)).first()

    human_ref_peptides = Channel.fromPath(file(params.human_ref_peptides)).first()
    
    //bwa_index = Channel.fromPath(params.bwa_index).collect()
    bwa_index = params.bwa_index


    // ****************************************************************
    // RUN FASTP QC ON ALL SAMPLES
    
    dna_fastp = DNA_QC_WORKFLOW(dna_inputs)
    rna_fastp = RNA_QC_WORKFLOW(rna_inputs)
    
    
    // Split Into

    num_intervals = params.scatter_count
    intervals_padded = SPLIT_INTERVALS_PADDED(reference_genome, intervals_file, num_intervals)
        .flatten()
        .map { file -> tuple(file.baseName, file) }
    
    intervals = SPLIT_INTERVALS(reference_genome, intervals_file, num_intervals)
        .flatten()
        .map { file -> tuple(file.baseName, file) }


    preproc_bam_workflow = DNA_ALIGNMENT_WORKFLOW(
        dna_fastp.fastp_fastqs,
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

    preproc_bams = preproc_bam_workflow.preproc_bams
    markdup_bams = preproc_bam_workflow.markdup_bams
    base_recal = preproc_bam_workflow.base_recal

    
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

    hla_fastq_input = dna_fastp.fastp_fastqs.mix(combined_fastqs)

   
    hla_workflow = HLA_TYPING_WORKFLOW(
        hla_fastq_input
    )

    /*
    // ****************************************************************
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
    
    // ****************************************************************
    // ALIGNMENT AND PREPROCESSING OF DNA

    if (params.bwa_index) {
        bwa_index = Channel.fromPath(params.bwa_index).collect()
    } else {
        bwa_index = CREATE_BWA_INDEX(reference_fa, reference_index_files)
    }

    bwa_sam = BWA_MAP(fastp_by_molec.dna, reference_fa, reference_index_files, bwa_index)
    bwa_mapped = BWA_POSTPROCESS(bwa_sam)

    
    // BAM PREPROCESSING OF DNA: GATK BEST PRACTICES

    preprc = PREPROCESS_BAM(bwa_mapped,
                            reference_fa, reference_index_files,
                            known_sites_dbsnp, known_sites_dbsnp_index,
                            known_sites_1000g_snps, known_sites_1000g_snps_index,
                            known_indels, known_indels_index,
                            mills, mills_index,
                            common_germline, common_germline_index)


    preproc_bams_type_branched = preprc.preproc_bams
                                   | branch {meta, bam, bai ->
                                        normal: meta.sample_type == "Normal"
                                        tumor: meta.sample_type == "Tumor"
                                        }
    
    preproc_bams_by_sample = preprc.preproc_bams.map { meta, bam, bai -> tuple(meta.somatic_sample, [meta, bam, bai]) } 
                                                     | groupTuple
                                                     | map { somatic_id, samples ->
                                                                def tumor  = samples.find { it[0].sample_type == 'Tumor' }
                                                                def normal = samples.find { it[0].sample_type == 'Normal' }
                                                                                     
                                                                return [ somatic_id, 
                                                                         tumor[0], tumor[1], tumor[2],
                                                                         normal[0], normal[1], normal[2]]
                                                     }
    

    // ****************************************************************
    // ALIGNMENT OF RNA SEQUENCING USING STAR
    
    if (params.star_index) {
        star_index = Channel.fromPath(params.star_index)
    }
    else {
        star_index = CREATE_STAR_INDEX(reference_fa, reference_index_files, gencode_gtf).star_index
    }
    
    star_rna_align = STAR_ALIGN(fastp_by_molec.rna, params.star_index)
    star_rna = STAR_INDEX_BAM(star_rna_align)

    // ****************************************************************
    // Transcript-fusion detection using STAR_FUSION and ARRIBA

    star_fusion = STAR_FUSION(star_rna_align.chimeric_out, ctat_resource_dir)
    
    arriba_fusions = ARRIBA_FUSION(star_rna_align.star_bam,
                                    reference_fa,
                                    reference_index_files,
                                    gencode_gtf,
                                    arriba_blacklist,
                                    arriba_known_fusions,
                                    arriba_protein_domains)

    
    // ****************************************************************

    // RNA: TRANSCRIPT ABUNDANCE ESTIMATION
    
    if (params.kallisto_index) {
        kallisto_index = params.kallisto_index
    } else {
        kallisto_index = KALLISTO_INDEX(kallisto_reference)

    }

    kallisto = KALLISTO_QUANT(fastp_by_molec.rna,
                              kallisto_index)

    kallisto = kallisto.abundance

    kallisto_gene = KALLISTO_TXIMPORT(kallisto, gencode_gtf)

    
    // ****************************************************************
    // INTERVAL CREATION FOR VARIANT CALLING

    intervals = SPLIT_INTERVALS(reference_fa, reference_index_files,
                                intervals_file, params.scatter_count)


    
   
    
    // ****************************************************************

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

    
    // ****************************************************************
    // SOMATIC VARIANT CALLING
    

    
    // MUTECT2
   
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
   
    

    // ****************************************************************
    // POST-PROCESS SOMATIC VCFS

    // Variant Decomposition & Normalization
    somatic_vcfs_vt = VT_SOMATIC_POSTPROCESS(somatic_filtered_vcfs, reference_fa, reference_index_files)
        | map { meta, vcf, tbi  -> tuple(meta.somatic_name, [meta, vcf, tbi]) }
        | groupTuple
        

    merged_vcf = MERGE_SOMATIC_VCFS(somatic_vcfs_vt,
                                      reference_fa,
                                      reference_index_files)

    
    // ****************************************************************
    // PVACtools VCF PREPARATION
    
    // Annotate VCF with VEP
    vep_annot = VEP_ANNOTATE(merged_vcf,
                       reference_fa,
                       params.vep_cache,
                       params.vep_plugins)

    vep = VEP_FILTER(vep_annot, params.vep_cache, params.vep_plugins)


    // Put all bams together (DNA + RNA) for input into bamreadcount
    all_bams = preprc.preproc_bams.concat(star_rna.star_bam)
        | map { meta, bam, bai -> tuple(meta.sample_name, [meta, bam, bai]) }


    vcf_prepared = vep.flatMap { vcf_meta, vcf ->
        [
            [ vcf_meta.tumor_metamap,  vcf_meta, vcf ],
            [ vcf_meta.normal_metamap, vcf_meta, vcf ]
        ]
    } 
    | map { sample_meta, vcf_meta, vcf  ->
                tuple(sample_meta.sample_name, [vcf_meta, vcf])
    }
    | combine(all_bams, by:0)
    | map { sample_name, vcf_list, bam_list ->
            def (vcf_meta, vcf) = vcf_list
            def (bam_meta, bam, bai) = bam_list
            return [vcf_meta, vcf, bam_meta, bam, bai]
    }



    // Add Read Coverage to VCF with Bamreadcount
    bamreadcount = BAMREADCOUNT(vcf_prepared,
                              reference_fa)

    bamreadcount_by_sample = bamreadcount.map { vcf_meta, vcf, bam_meta, brc_indel, brc_snv ->
                                                tuple (vcf_meta.somatic_name, [vcf_meta, vcf, bam_meta, brc_indel, brc_snv])
                                } | groupTuple 
                                | map { somatic_name, samples  ->

                                    def tumor = samples.find { it[2].sample_type == "Tumor" && it[2].molecule == 'DNA' }
                                    def normal = samples.find { it[2].sample_type == "Normal" && it[2].molecule == "DNA" }
                                    def tumor_rna = samples.find { it[2].sample_type == "Tumor" && it[2].molecule == "RNA" }

                                    def vcf_metadata = tumor[0]
                                    def vcf_file     = tumor[1]
                                    
                                    def tumor_files  = [ tumor[3], tumor[4] ]
                                    def normal_files = [ normal[3], normal[4] ]
                                    def tumor_rna_files    = [ tumor_rna[3], tumor_rna[4] ]
                                    
                                    return tuple(somatic_name, vcf_metadata, vcf_file, tumor_files, normal_files, tumor_rna_files)
                                }



    vcf_rna_annotated_coverage = ANNOTATE_VCF_COVERAGE(bamreadcount_by_sample)
    
    
    // Add transcript abundance estimation from kallisto
    

    kallisto = kallisto | map { meta, abundance -> tuple(meta.sample_name, [meta, abundance]) }
    kallisto_gene = kallisto_gene | map { meta, gene_tpm -> tuple (meta.sample_name, [meta, gene_tpm]) }
    vcf_annot = vcf_rna_annotated_coverage | map { meta, vcf -> tuple(meta.tumor_metamap.sample_name, [meta, vcf]) }
    
    vcf_annotated = vcf_annot
        | join(kallisto)
        | join(kallisto_gene)
        | map { id, vcf_info, kallisto_tx, kallisto_gene ->
                def (vcf_meta, vcf) = vcf_info
                def (kallisto_tx_meta, tx_abundance) = kallisto_tx
                def (kallisto_gene_meta, gene_abundance) = kallisto_gene
                return [vcf_meta, vcf, kallisto_tx_meta, tx_abundance, gene_abundance ]
        }

    vcf_annotated_expression = ANNOTATE_VCF_EXPRESSION(vcf_annotated)

    
    // Index final somatic vcf for input to PVACseq
    vcf_final = INDEX_FINAL_VCF(vcf_annotated_expression)
    vcf_tables = VCF_TO_TABLE(vcf_final)
    
    
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

            
    
    // ****************************************************************
    // RUN PVACTOOLS SUITE (PVACSEQ + PVACFUSE)

    PVACSEQ = PVACSEQ(pvacseq_input, human_ref_peptides)


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


    PVACFUSE = PVACFUSE(pvacfuse_input, human_ref_peptides)

    */

    publish:
        dna_fastp_fastqs = dna_fastp.fastp_fastqs
}

output {
    dna_fastp_fastqs{
        path { meta, fastq1, fastq2 -> "${params.outdir}/${meta.somatic_name}/qc/fastp/${meta.sample_name}" }
    }

}
