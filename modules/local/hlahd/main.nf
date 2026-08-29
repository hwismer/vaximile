process HLAHD {

    /*

    Run HLA-HD to perform HLA typing on patient fastqs

    This has given me some trouble before with temp directories and issues where bowtie never initializes.
    This should be fixed by unzipping fastqs at the beginning of the script.
    One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.


    */
    
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
        path "versions.yml", topic: versions

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
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hlahd: 1.0
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        mkdir -p ${prefix}/result
        touch ${prefix}/result/${prefix}_final.result.txt
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hlahd: 1.0
        END_VERSIONS
        """

}
