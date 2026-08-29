process KALLISTO_QUANT {


    /*

        Get transcript abundance estimations using kallisto quant.

    */

    label 'process_high'

    conda "bioconda::kallisto=0.51.1"
    
    tag "Kallisto quant on ${meta.sample_name}"

    input:
        tuple val(meta), path(read1), path(read2)
        path kallisto_index

    output:
        tuple val(meta), path("*_kallisto/abundance.tsv"), emit: abundance
        tuple val(meta), path("*_kallisto"), emit: kallisto_dir
        // Removed: `path("*")` globbed the whole work directory, so staged inputs
        // (FASTQs, index dirs, the decompressed GTF) and versions.yml were emitted
        // as results. Nothing consumed it.
        path "versions.yml", topic: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        kallisto quant -i $kallisto_index -o ${prefix}_kallisto -t ${task.cpus} $read1 $read2 > "${prefix}_kallist_stdout.out"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallisto: \$(kallisto version 2>&1 | sed 's/kallisto, version //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        mkdir -p ${prefix}_kallisto
        touch ${prefix}_kallisto/abundance.tsv
        touch ${prefix}_kallist_stdout.out
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallisto: 0.51.1
        END_VERSIONS
        """
}
