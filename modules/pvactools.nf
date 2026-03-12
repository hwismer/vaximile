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

    cpus 12
    memory "128GB"

    container "griffithlab/pvactools:6.0.3"

    publishDir "${params.outdir}/${somatic_meta.somatic_name}/pvactools/", mode: "copy"

    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),path(phased_vcf), path(phased_vcf_index), path(hla_pvac_input)
        path(human_ref_peptides)
    output:
        path("${somatic_meta.somatic_name}_pvacseq"), emit: pvacseq_dir

    script:
        """
        pvacseq run \
            $somatic_vcf \
            ${somatic_meta.tumor_metamap.sample_name} \
            \$(head $hla_pvac_input -n 1) \
            all \
            "${somatic_meta.somatic_name}_pvacseq" \
            -e1 8,9,10,11 \
            -e2 12,13,14,15,16,17,18 \
            --phased-proximal-variants-vcf $phased_vcf \
            --normal-sample-name ${somatic_meta.normal_metamap.sample_name} \
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
    
    cpus 12
    memory "128GB"

    container "griffithlab/pvactools:6.0.3"

    publishDir "${params.outdir}/${sample_meta.somatic_sample}/pvactools/", mode: "copy"

    input:
        tuple val(sample_meta), path(arriba_fusions), path(hla_pvac_input), path(starfusion_calls)
        path human_ref_peptides

    output:
        path("${sample_meta.somatic_sample}_pvacfuse"), emit: pvacfuse_dir
    
    script:
        """
        pvacfuse run \
            $arriba_fusions \
            "${sample_meta.somatic_sample}" \
            \$(head $hla_pvac_input -n 1) \
            all \
            "${sample_meta.somatic_sample}_pvacfuse" \
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
