process BAMREADCOUNT {

    /*

    Run bamreadcount to add coverage information to a vcf using a BAM file. The sample names in the VCF
    must match the read groups and sample names in the BAM.

    */

    cpus 4
    memory "32GB"
    cache "lenient"

    container "mgibio/bam_readcount_helper-cwl:1.2.1"

    input:
        tuple val(somatic_name), val(somatic_meta), path(vcf), val(sample_meta), path(bam), path(bai)
        tuple path(reference_fa), path(reference_index), path(reference_dict)

    output:

        tuple val(somatic_meta), val(sample_meta), path("${sample_meta.sample_name}_${sample_meta.molecule}_bamrc_helper/*indel.tsv"), path("${sample_meta.sample_name}_${sample_meta.molecule}_bamrc_helper/*snv.tsv"), emit: brc_files
    script:
        
        """
        mkdir ${sample_meta.sample_name}_${sample_meta.molecule}_bamrc_helper
        bam_readcount_helper.py \
            $vcf \
            ${sample_meta.sample_name} \
            $reference_fa \
            $bam \
            ${sample_meta.molecule} \
            ${sample_meta.sample_name}_${sample_meta.molecule}_bamrc_helper
        """

}

process ANNOTATE_VCF_COVERAGE {

    /*

        Use the output of bamreadcount to annotate coverage given a particular sample.
    */

    cpus 4
    memory "16GB"

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_meta), val(sample_meta), path(indels), path(snvs), path(vcf)
    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${sample_meta.sample_name}_${sample_meta.molecule}_coverage.vcf")
        
    script:
        """
        vcf-readcount-annotator \
            $vcf \
            $indels \
            ${sample_meta.molecule} \
            -s ${sample_meta.sample_name} \
            -t indel \
            -o vcf1.vcf

        vcf-readcount-annotator \
            vcf1.vcf \
            $snvs \
            ${sample_meta.molecule} \
            -s ${sample_meta.sample_name} \
            -t snv \
            -o ${somatic_meta.somatic_name}_${sample_meta.sample_name}_${sample_meta.molecule}_coverage.vcf
        """

}

process VEP_ANNOTATE {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    cpus 4
    memory "16GB"
    
    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(somatic_meta), path(vcf)
        tuple path(reference_fa), path(reference_index), path(reference_dict)
        path vep_cache
        path vep_plugins

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_vep.vcf"), emit: vcf
        tuple val(somatic_meta), path("*.html"), emit: report

    script:
        """
        vep \
            --input_file $vcf  \
            --output_file ${somatic_meta.somatic_name}_vep.vcf \
            --everything \
            --format vcf --vcf --symbol --terms SO --mane_select --canonical --tsl --biotype --hgvs \
            --fasta $reference_fa  \
            --offline --cache \
            --plugin Frameshift --plugin Wildtype \
            --pick \
            --dir_plugins $vep_plugins \
            --dir_cache $vep_cache \
            --transcript_version
        """
}

process VEP_POPULATION_FILTER {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    cpus 2
    memory "16GB"
    cache "lenient"


    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(somatic_meta), path(vep_vcf)
        path vep_cache
        path vep_plugins

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_vep_filter.vcf"), emit: vcf

    script:
        """
        filter_vep -i $vep_vcf \
            -o "${somatic_meta.somatic_name}_vep_filter.vcf" \
            --format vcf \
            --filter "gnomADe_AF < 0.001 or not gnomADe_AF"
        """
}

process ANNOTATE_VCF_TRANSCRIPT_EXPRESSION {

    /*

    Use the abundance estimates from kallist to annotate transcript expression in a vcf file.

    */
    cpus 2
    memory "16GB"

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_name), val(somatic_meta), path(vcf), val(sample_meta), path(tx_abundance)

    output:
        tuple val(somatic_meta), path("${somatic_name}_tx_expression.vcf")

    script:

        """
        vcf-expression-annotator \
            $vcf \
            -s ${sample_meta.sample_name} \
            $tx_abundance \
            kallisto transcript \
            -o ${somatic_name}_tx_expression.vcf
        """

}

process ANNOTATE_VCF_GENE_EXPRESSION {

    /*

    Use the abundance estimates from kallist to annotate transcript expression in a vcf file.

    */
    cpus 2
    memory "16GB"

    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_name), val(somatic_meta), path(vcf), val(sample_meta), path(gene_abundance)


    output:
        tuple val(somatic_meta), path("${somatic_name}_gene_expression.vcf")

    script:

        """
        vcf-expression-annotator \
            $vcf \
            $gene_abundance \
            custom gene \
            -i ENSEMBLID \
            -e TPM \
            -s ${sample_meta.sample_name} \
            --ignore-ensembl-id-version \
            -o "${somatic_name}_gene_expression.vcf"
        """

}

