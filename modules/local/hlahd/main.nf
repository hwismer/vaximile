process HLAHD {

    // HLA typing with HLA-HD. FASTQs are decompressed first to avoid bowtie start-up failures.
    
    label 'process_high'
    container "griffithlab/hlahd:1.0"

    tag "HLA-HD on ${meta.sample_name}"
    
    input:
        tuple val(meta), path(fastq1), path(fastq2)
    
    output:
        tuple val(meta), path("*/result/*_final.result.txt"), emit: final_hla_calls
        // No trailing slash: a glob ending in "/" never matches, so `*/result/` fails
        // with "Missing output file(s)" even when the directory is there.
        tuple val(meta), path("*/result", type: 'dir'), emit: result_dir

    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        /opt/hlahd/bin/hlahd.1.6.1.sh \
            -f /opt/hlahd/freq_data \
            -t $task.cpus \
            $fastq1 $fastq2 \
            /opt/hlahd/HLA_gene.split.txt \
            /opt/hlahd/dictionary \
            "${prefix}" \
            ./
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        mkdir -p ${prefix}/result
        touch ${prefix}/result/${prefix}_final.result.txt
        """

}
