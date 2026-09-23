process LOHHLA_MQC {

    // Loads LOHHLAmod results into MultiQC

    label 'process_single'

    conda "conda-forge::poppler=26.09.0 conda-forge::pillow=12.3.0 python=3.12"

    tag "Preparing ${somatic_name} HLA LOH results for MultiQC"

    input:
        tuple val(somatic_name), val(meta), path(plot_dir), path(loh_res)
        path montage_script

    output:
        tuple val(somatic_name), val(meta), path("*_hla_loh_mqc.png"), emit: png
        tuple val(somatic_name), val(meta), path("*_lohres.tsv"), emit: tsv

    script:
    // Rasterise the LOH PDFs and tile them into one sheet per pair for MultiQC (genes across, plots down).
    // The tiling is assets/montage_panels.py, staged as an input.
    """
    mkdir -p panels
    for pdf in ${plot_dir}/*.pdf; do
        [ -e "\$pdf" ] || continue
        pdftoppm -png -r 100 -singlefile "\$pdf" "panels/\$(basename "\$pdf" .pdf)"
    done

    python3 $montage_script panels/*.png \\
        --out ${somatic_name}_hla_loh_mqc.png \\
        --layout grid \\
        --row-order t_dp,n_dp,tn_dp,logR,baf

    # The result table is keyed by HLA gene alone - hla_a, hla_b, hla_c - which is unique
    # within a pair and not across them. MultiQC reads the first column as the row key and
    # merges every file it matches into one table, so without the pair in the key a second
    # pair overwrites the first rather than appending to it.
    awk -F'\t' -v OFS='\t' -v pair="${somatic_name}" '
        NR == 1 { \$1 = "Pair_HLAGene"; print; next }
        { \$1 = pair "_" \$1; print }
    ' ${loh_res} > ${somatic_name}_lohres.tsv

    """

    stub:
    """
    touch ${somatic_name}_hla_loh_mqc.png
    touch ${somatic_name}_lohres.tsv
    """
}
