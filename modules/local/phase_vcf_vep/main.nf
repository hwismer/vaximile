process PHASE_VCF_VEP {

    /*

        Use VEP to annotated the phased vcf.

    */
    cpus 8
    memory "32GB"
    container "ensemblorg/ensembl-vep:release_115.0"

    tag "VEP on phased vcf $phased_vcf"

    input:
        tuple val(meta), path(phased_vcf)
        path reference_fa
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("phased_vcf_vep.vcf"), emit: vcf

    script:
        """
        vep \
            --input_file $phased_vcf  \
            --output_file "phased_vcf_vep.vcf" \
            --format vcf --vcf --symbol --terms SO --tsl --biotype \
            --hgvs --fasta $reference_fa  \
            --offline --cache $vep_cache \
            --plugin Frameshift --plugin Wildtype --plugin Downstream \
            --pick \
            --fork ${task.cpus} \
            --dir_plugins $vep_plugins \
            --transcript_version

        """

}
