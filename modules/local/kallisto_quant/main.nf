process KALLISTO_QUANT {


    /*

        Get transcript abundance estimations using kallisto quant.

    */

    cpus 8
    memory "32GB"

    conda "bioconda::kallisto=0.51.1"
    
    tag "Kallisto quant on ${meta.sample_name}"

    input:
        tuple val(meta), path(read1), path(read2)
        path kallisto_index

    output:
        tuple val(meta), path("${meta.sample_name}_kallisto/abundance.tsv"), emit: abundance
        tuple val(meta), path("${meta.sample_name}_kallisto"), emit: kallisto_dir
        tuple val(meta), path("*"), emit: tutto

    script:
        """
        kallisto quant -i $kallisto_index -o ${meta.sample_name}_kallisto -t ${task.cpus} $read1 $read2 > "${meta.sample_name}_kallist_stdout.out" 
        """
}
