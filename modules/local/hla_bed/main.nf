process HLA_BED {
    // Create a BED file containing the regions of GRCh38 that harbor HLA genes.
    // For downstream use with mhcflow.
    
    label 'process_single'

    conda "conda-forge::coreutils=9.3"

    input:
    val chr_prefix

    output:
    path "hla_region.bed", emit: bed
    
    script:
    def contig = chr_prefix ? 'chr6' : '6'
    def start = 28510119
    def end = 33480577
    """
    printf '%s\\t%d\\t%d\\t%s\\n' '${contig}' ${start} ${end} 'MHC' > hla_region.bed
    """

    stub:
    """
    touch hla_region.bed
    """
}
