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

    cpus 16
    memory "64GB"
    container "griffithlab/pvactools:6.0.3"

    tag "pVACseq on ${somatic_name}"
    input:
        tuple val(somatic_name), 
            val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            path(phased_vcf), path(phased_vcf_index), 
            val(hla_meta), path(hla_pvac_input)
        path(human_ref_peptides)
    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_pvacseq"), emit: pvacseq_dir

    script:
        """
        pvacseq run \
            $somatic_vcf \
            ${somatic_meta.tumor_meta.sample_name} \
            \$(head $hla_pvac_input -n 1) \
            all \
            "${somatic_meta.somatic_name}_pvacseq" \
            -e1 8,9,10,11 \
            -e2 12,13,14,15,16,17,18 \
            --phased-proximal-variants-vcf $phased_vcf \
            --normal-sample-name ${somatic_meta.normal_meta.sample_name} \
            --iedb-install-directory /opt/iedb \
            --pass-only \
            --run-reference-proteome-similarity \
            --peptide-fasta $human_ref_peptides \
            -m2 percentile \
            -m median \
            -a sample_name \
            --problematic-amino-acids P:-2 \
            -t $task.cpus
        """
}

process PVACFUSE {
    
    cpus 8
    memory "64GB"
    container "griffithlab/pvactools:6.0.3"

    
    tag "pVACfuse on ${somatic_name}"
    input:
        tuple val(somatic_name), val(arriba_meta), path(arriba_fusions), val(star_meta), path(starfusion_calls), val(hla_meta), path(hla_pvac_input)
        path human_ref_peptides

    output:
        tuple val(star_meta), path("${somatic_name}_pvacfuse"), emit: pvacfuse_dir
    
    script:
        """
        pvacfuse run \
            $arriba_fusions \
            "${somatic_name}" \
            \$(head $hla_pvac_input -n 1) \
            all \
            "${somatic_name}_pvacfuse" \
            --starfusion-file $starfusion_calls \
            -e1 8,9,10,11 \
            -e2 12,13,14,15,16,17,18 \
            --iedb-install-directory /opt/iedb \
            --run-reference-proteome-similarity \
            --peptide-fasta $human_ref_peptides \
            -m median \
            -m2 percentile \
            -t $task.cpus
        """

}
