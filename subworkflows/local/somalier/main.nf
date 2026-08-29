include { SOMALIER_EXTRACT } from "../../../modules/local/somalier_extract/main"
include { SOMALIER_RELATE } from "../../../modules/local/somalier_relate/main"

workflow SOMALIER {

    take:
        sample_bams
        reference_genome
        somalier_sites_vcf

    main:
        extract = SOMALIER_EXTRACT(sample_bams, reference_genome, somalier_sites_vcf).somalier
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
