process PVACSEQ {

    /*

    Run PVACseq on a single sample.

    Inputs:
        *(See PVACseq input preparation documents for more info.)
        Somatic VCF
        Phased Germline VCF
        Tumor Sample Metadata
        Normal Sample Metadata

    Output:
        Pvacseq folder containing:
            Combined neoantigen predictions
            MHC I neoantigen predictions
            MHC II neoantigen predictions

    */

    label 'process_very_high'
    container "griffithlab/pvactools:7.0.1"

    tag "pVACseq on ${somatic_name}"
    input:
        tuple val(somatic_name), 
            val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
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
            ${somatic_meta.tumor_meta.sample_name} \
            \$(head $hla_pvac_input -n 1) \
            all \
            "${prefix}_pvacseq" \
            $args \
            --phased-proximal-variants-vcf $phased_vcf \
            --normal-sample-name ${somatic_meta.normal_meta.sample_name} \
            --iedb-install-directory /opt/iedb \
            --peptide-fasta $human_ref_peptides \
            -t $task.cpus
        """
}
