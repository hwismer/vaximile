include { FASTP } from "../modules/quality_control.nf"
include { SOMALIER_EXTRACT; SOMALIER_RELATE } from "../modules/quality_control.nf"

workflow DNA_QC_WORKFLOW {

    take:
        reads_channel
        
    main:
        fastp = FASTP(reads_channel)

    emit:
        fastp_fastqs = fastp.fastqs
        fastp_reports = fastp.reports
}

workflow RNA_QC_WORKFLOW {

    take:
        reads_channel

    main:
        fastp = FASTP(reads_channel)

    emit:
        fastp_fastqs = fastp.fastqs
        fastp_reports = fastp.reports

}

workflow SOMALIER {

    take:
        sample_bams
        reference_genome
        somalier_sites_vcf

    main:
        extract = SOMALIER_EXTRACT(sample_bams, reference_genome, somalier_sites_vcf)
        extract_unique = extract.unique { meta, files -> 
            meta.sample_name
        }
        extract_by_patient = extract_unique.map{meta, extracted ->
            tuple(meta.patient, meta, extracted)
        }
        .groupTuple()

        relate = SOMALIER_RELATE(extract_by_patient)
    
    emit:
        extract = extract
        pairs = relate.pairs
        samples = relate.samples
        groups = relate.groups
        html = relate.html
        

}
