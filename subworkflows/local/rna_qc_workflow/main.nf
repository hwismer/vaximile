include { FASTP } from "../../../modules/local/fastp/main"

workflow RNA_QC_WORKFLOW {

    take:
        reads_channel

    main:
        fastp = FASTP(reads_channel)

    emit:
        fastp_fastqs = fastp.fastqs
        fastp_reports = fastp.reports

}
