process MHCFLOW_REALIGN {

    label 'process_max'

    conda "${moduleDir}/environment.yml"

    tag "Realigning ${meta.sample_name} to ${somatic_name} normal HLA reference"

    input:
        // The tumour library, and the HLA reference mhcflow inferred for the normal it is
        // paired with. Step 2 of the LOH workflow: the tumour is realigned against the
        // NORMAL's alleles rather than the full HLA reference, because lohhlamod compares
        // coverage between the two BAMs and that is only meaningful over one common,
        // subject-specific set of sequences.
        tuple val(somatic_name), val(meta), path(bam), path(bai), path(hla_fasta), path(hla_nix)
        path hla_bed
        path hla_kmers
        path hla_freqs

    output:
        // realigner/, not finalizer/: --realn-only returns before typing, so this run has
        // no finalizer stage. The reference it aligned against is the normal's, which is
        // what makes this BAM comparable to the normal's own finalizer BAM.
        tuple val(somatic_name), val(meta), path("${meta.sample_name}_realn/realigner/*.hla.realn.bam"), path("${meta.sample_name}_realn/realigner/*.hla.realn.bam.bai"), emit: realn_bam
        path "versions.yml", topic: versions

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
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unknown: unknown
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.sample_name}_realn/realigner
    touch ${meta.sample_name}_realn/realigner/${meta.sample_name}.hla.realn.bam
    touch ${meta.sample_name}_realn/realigner/${meta.sample_name}.hla.realn.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unknown: unknown
    END_VERSIONS
    """
}
