include { CREATE_BWA_INDEX } from "../../../modules/local/create_bwa_index/main"

workflow BWA_INDEX {
    
    take:
        reference_genome
        bwa_index

    main:
        if ( bwa_index ) {
           bwa_index_ch = Channel.fromPath("${bwa_index}/*").collect()
        } else {
            bwa_index_ch = CREATE_BWA_INDEX(reference_genome)
        }


    emit:
        bwa_index = bwa_index_ch

}
