process INDEX_FASTA {

    // Indexes a fasta file with samtools faidx

    label 'process_medium'
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    tag "Indexing $fasta"

    input:
    	path fasta

    output:
    	tuple path(fasta), path("${fasta}.fai"), emit: fai
    	path "versions.yml", topic: versions

    script:
        """
        set -euo pipefail
    	samtools faidx $fasta
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
        END_VERSIONS
        """

    stub:
        """
        touch ${fasta}.fai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: 1.23.1
        END_VERSIONS
        """
}
