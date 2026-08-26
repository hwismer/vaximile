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
        tuple val(meta), path("*"), emit: tutto

    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        kallisto quant -i $kallisto_index -o ${prefix}_kallisto -t ${task.cpus} $read1 $read2 > "${prefix}_kallist_stdout.out"
        """
}
