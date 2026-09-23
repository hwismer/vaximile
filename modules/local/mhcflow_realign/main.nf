process MHCFLOW_REALIGN {

    label 'process_max'

    conda "${moduleDir}/environment.yml"

    tag "Realigning ${meta.sample_name} to ${somatic_name} normal HLA reference"

    input:
        // Realign the tumour against its normal's HLA reference, so both BAMs share one reference for LOH.
        tuple val(somatic_name), val(meta), path(bam), path(bai), path(hla_fasta), path(hla_nix)
        path hla_bed
        path hla_kmers
        path hla_freqs

    output:
        // realigner/, not finalizer/: --realn-only returns before typing, so this run has
        // no finalizer stage. The reference it aligned against is the normal's, which is
        // what makes this BAM comparable to the normal's own finalizer BAM.
        tuple val(somatic_name), val(meta), path("${meta.sample_name}_realn/realigner/*.hla.realn.bam"), path("${meta.sample_name}_realn/realigner/*.hla.realn.bam.bai"), emit: realn_bam

    script:
    def args = task.ext.args ?: ''
    """
    mhcflow --bam $bam \
        --ref $hla_fasta \
        --bed $hla_bed \
        --tag $hla_kmers \
        --freq $hla_freqs \
        --nproc $task.cpus \
        --realn-only \
        $args \
        --outdir ${meta.sample_name}_realn
    """

    stub:
    """
    mkdir -p ${meta.sample_name}_realn/realigner
    touch ${meta.sample_name}_realn/realigner/${meta.sample_name}.hla.realn.bam
    touch ${meta.sample_name}_realn/realigner/${meta.sample_name}.hla.realn.bam.bai
    """
}
