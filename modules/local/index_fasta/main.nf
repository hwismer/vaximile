process INDEX_FASTA {

    // Indexes a fasta file with samtools faidx

    label 'process_medium'
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    tag "Indexing $fasta"

    input:
    	path fasta

    output:
    	tuple path(fasta), path("${fasta}.fai"), emit: fai

    script:
        """
        set -euo pipefail
    	samtools faidx $fasta
        """

    stub:
        """
        touch ${fasta}.fai
        """
}
