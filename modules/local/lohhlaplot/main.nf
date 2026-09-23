process LOHHLAPLOT {

    label 'process_single'

    container "hmitzinger/lohhlamod@sha256:8d37ec357a5a5e4f04b6adbd7c69fc73156273527de0132a8f7d5eeb7398a293"

    tag "Plotting HLA LOH profiles for ${somatic_name}"

    input:
        // The result table and the per-gene .rds files, staged flat into the task directory
        // rather than as the LOHHLAMOD output directory. lohhlaplot writes its plots INTO
        // --loh-dir, and that directory would be a symlink to the LOHHLAMOD task: plotting
        // would write into another task's work directory, which breaks under `scratch` and
        // leaves one task's outputs depending on another task not having been cleaned.
        // Staged flat, --loh-dir is this task's own directory and the plots land here.
        tuple val(somatic_name), val(meta), path(loh_res), path(loh_rds)

    output:
        tuple val(somatic_name), val(meta), path("${somatic_name}_plots"), emit: plots

    script:
    def args = task.ext.args ?: ''
    """
    lohhlaplot --sample ${somatic_name} \\
        --loh-res $loh_res \\
        --loh-dir . \\
        $args

    """

    stub:
    """
    mkdir -p ${somatic_name}_plots
    touch ${somatic_name}_plots/hla_a.logR.pdf
    """
}
