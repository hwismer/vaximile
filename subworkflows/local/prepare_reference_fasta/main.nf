include { PREPARE_FASTA } from "../../../modules/local/prepare_fasta/main"
include { INDEX_FASTA } from "../../../modules/local/index_fasta/main"
include { MAKE_FASTA_DICT } from "../../../modules/local/make_fasta_dict/main"

workflow PREPARE_REFERENCE_FASTA {
    
    take:
        fasta
    main:
        
        fasta_proc = PREPARE_FASTA(fasta)
        fasta_plus_fai = INDEX_FASTA(fasta_proc).fai
        dict = MAKE_FASTA_DICT(fasta_plus_fai).dict

    emit:
        fa_fai_pair = fasta_plus_fai
        dict  = dict
        
}
