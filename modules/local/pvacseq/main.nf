process PVACSEQ {

    // Run pVACseq on a tumour/normal pair's somatic and phased germline VCFs.

    label 'process_very_high'
    // Container only: pVACtools runs the IEDB predictors from /opt/iedb, which only the image has.
    container "griffithlab/pvactools:7.1.5"

    tag "pVACseq on ${somatic_name}"
    input:
        tuple val(somatic_name),
            val(somatic_meta), val(tumor_sample_name), val(normal_sample_name),
            path(somatic_vcf), path(somatic_vcf_index),
            path(phased_vcf), path(phased_vcf_index), 
            val(hla_meta), val(hla_alleles)
        path(human_ref_peptides)
    output:
        tuple val(somatic_meta), path("*_pvacseq"), emit: pvacseq_dir
        tuple val(somatic_meta), path("*.MHC_I.report.tsv"), emit: mhc_i_report
        tuple val(somatic_meta), path("*.MHC_II.report.tsv"), optional: true, emit: mhc_ii_report

    script:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        def args = task.ext.args ?: ''
        """
        pvacseq run \
            $somatic_vcf \
            ${tumor_sample_name} \
            "${hla_alleles}" \
            all \
            "${prefix}_pvacseq" \
            $args \
            --phased-proximal-variants-vcf $phased_vcf \
            --normal-sample-name ${normal_sample_name} \
            --iedb-install-directory /opt/iedb \
            --peptide-fasta $human_ref_peptides \
            -t $task.cpus

        # The class I aggregated report for later steps: the ML classifier's copy, which adds its
        # call, when it ran (it needs class II results too), otherwise pVACseq's own.
        report="${prefix}_pvacseq/MHC_Class_I/${tumor_sample_name}.MHC_I.all_epitopes.aggregated"
        if [ -e "\$report.ML_predict.tsv" ]; then
            cp "\$report.ML_predict.tsv" ${tumor_sample_name}.MHC_I.report.tsv
        else
            cp "\$report.tsv" ${tumor_sample_name}.MHC_I.report.tsv
        fi
        # The class II aggregated report, when class II alleles were given.
        report="${prefix}_pvacseq/MHC_Class_II/${tumor_sample_name}.MHC_II.all_epitopes.aggregated.tsv"
        if [ -e "\$report" ]; then
            cp "\$report" ${tumor_sample_name}.MHC_II.report.tsv
        fi
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        mkdir -p ${prefix}_pvacseq/MHC_Class_I
        touch ${prefix}_pvacseq/MHC_Class_I/${tumor_sample_name}.MHC_I.all_epitopes.aggregated.tsv
        touch ${tumor_sample_name}.MHC_I.report.tsv ${tumor_sample_name}.MHC_II.report.tsv
        """
}
