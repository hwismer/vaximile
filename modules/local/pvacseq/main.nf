process PVACSEQ {

    // Run pVACseq on a tumour/normal pair's somatic and phased germline VCFs.

    label 'process_very_high'
    conda "bioconda::pvactools=7.0.1"
    container "griffithlab/pvactools:7.0.1"

    tag "pVACseq on ${somatic_name}"
    input:
        tuple val(somatic_name),
            val(somatic_meta), val(tumor_sample_name), val(normal_sample_name),
            path(somatic_vcf), path(somatic_vcf_index),
            path(phased_vcf), path(phased_vcf_index), 
            val(hla_meta), path(hla_pvac_input)
        path(human_ref_peptides)
    output:
        tuple val(somatic_meta), path("*_pvacseq"), emit: pvacseq_dir
        tuple val(somatic_meta), path("*_pvacseq/MHC_Class_I/*MHC_I.all_epitopes.aggregated.tsv"), emit: pvaseq_mhc_i_aggr

    script:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        def args = task.ext.args ?: ''
        """
        pvacseq run \
            $somatic_vcf \
            ${tumor_sample_name} \
            \$(head $hla_pvac_input -n 1) \
            all \
            "${prefix}_pvacseq" \
            $args \
            --phased-proximal-variants-vcf $phased_vcf \
            --normal-sample-name ${normal_sample_name} \
            --iedb-install-directory /opt/iedb \
            --peptide-fasta $human_ref_peptides \
            -t $task.cpus
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        mkdir -p ${prefix}_pvacseq/MHC_Class_I
        touch ${prefix}_pvacseq/MHC_Class_I/${prefix}_MHC_I.all_epitopes.aggregated.tsv
        """
}
