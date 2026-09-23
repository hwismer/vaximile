process LOHHLAMOD {

    // Runs LOHHLA via LOHHLAmod

    label 'process_high'

    container "hmitzinger/lohhlamod@sha256:8d37ec357a5a5e4f04b6adbd7c69fc73156273527de0132a8f7d5eeb7398a293"

    tag "Detecting HLA LOH in ${somatic_name}"

    input:
        tuple val(somatic_name), val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(hla_fasta), path(purityploidy)

    output:
        tuple val(somatic_name), val(meta), path("${somatic_name}_lohhla"), emit: loh_dir
        tuple val(somatic_name), val(meta), path("${somatic_name}_lohhla/${somatic_name}.loh.res.tsv"), emit: loh_res
        tuple val(somatic_name), val(meta), path("${somatic_name}_lohhla/*.rds"), emit: loh_rds

    script:
    def args = task.ext.args ?: ''
    def tstates = purityploidy ? "--tstates tstates.tsv" : ''
    def build_tstates = purityploidy ? """
    printf 'SampleID\\tTumorPloidy\\tTumorPurityNGS\\n' > tstates.tsv
    awk 'NR == 2 { print "${somatic_name}" "\\t" \$2 "\\t" \$1 }' ${purityploidy} >> tstates.tsv
    """ : ''
    """
    $build_tstates
    lohhlamod --subject ${somatic_name} \\
        --tbam $tumor_bam \\
        --nbam $normal_bam \\
        --hlaref $hla_fasta \\
        $tstates \\
        --threads $task.cpus \\
        $args \\
        --outdir ${somatic_name}_lohhla

    """

    stub:
    """
    mkdir -p ${somatic_name}_lohhla
    touch ${somatic_name}_lohhla/${somatic_name}.loh.res.tsv
    touch ${somatic_name}_lohhla/hla_a.rds ${somatic_name}_lohhla/hla_b.rds ${somatic_name}_lohhla/hla_c.rds
    """
}
