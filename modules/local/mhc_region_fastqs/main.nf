process MHC_REGION_FASTQS {
    
    label 'process_medium'

    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"
    
    tag "Extracting MHC regions from ${meta.sample_name}"

    input:
    	tuple val(meta), path(bam), path(bai)

    output:
    	tuple val(meta), path("*_R1.fastq"), path("*_R2.fastq"), emit: reads

    script:
    def mhc_region = 'chr6:28510120-33480577'
    def prefix = task.ext.prefix ?: "${meta.sample_name}"

    """
    set -euo pipefail

    # Build BED of non-primary contigs
    samtools idxstats --threads $task.cpus $bam \
      | awk '
          \$1 != "*" &&
          \$1 != "chr6" &&
          \$1 !~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/ &&
          \$1 != "MT"
        ' \
      | awk '{print \$1 "\\t0\\t" \$2}' \
      > nonprimary.bed

    # chr6 MHC interval
    samtools view --threads $task.cpus \
        -b \
        -F 0x904 \
        ${bam} \
        ${mhc_region} \
        > chr6_mhc.bam

    # unmapped reads
    samtools view --threads $task.cpus \
        -b \
        -f 4 \
        -F 0x904 \
        ${bam} \
        > unmapped.bam

    # reads on non-primary contigs
    samtools view --threads $task.cpus \
        -b \
        -F 0x904 \
        -L nonprimary.bed \
        ${bam} \
        > nonprimary.bam

    # merge selected reads
    samtools merge \
        -f \
        merged.bam \
        chr6_mhc.bam \
        unmapped.bam \
        nonprimary.bam

    # BAM -> paired FASTQ
    samtools collate -Ou merged.bam \
      | samtools fastq \
            -1 ${prefix}_R1.fastq \
            -2 ${prefix}_R2.fastq \
			-0 /dev/null \
			-s /dev/null \
            -
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    touch ${prefix}_R1.fastq
    touch ${prefix}_R2.fastq
    """
}
