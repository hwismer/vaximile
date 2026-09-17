process LOHHLA_MQC {

    label 'process_single'

    conda "conda-forge::poppler=26.09.0"

    tag "Preparing ${somatic_name} HLA LOH results for MultiQC"

    input:
        tuple val(somatic_name), val(meta), path(plot_dir), path(loh_res)

    output:
        // MultiQC picks up any *_mqc.png as a custom-content image and titles the section
        // from the file name, so the pair and the plot are encoded there: it is the only
        // place that survives into the report.
        tuple val(somatic_name), val(meta), path("*_mqc.png"), emit: png
        tuple val(somatic_name), val(meta), path("*_lohres.tsv"), emit: tsv
        path "versions.yml", topic: versions

    script:
    // lohhlaplot writes PDFs, which MultiQC cannot embed - it takes png, jpg and jpeg only -
    // so they are rasterised here rather than lost from the report. -r 150 keeps the axis
    // labels legible at report width without making the HTML enormous; -singlefile because
    // each plot is one page and pdftoppm otherwise appends a page number to the name.
    //
    // The dots in lohhlaplot's names (hla_a.logR.pdf) are translated to underscores, and
    // that is load-bearing rather than cosmetic. MultiQC derives an image section's name
    // from the file name and strips what looks like an extension off it, using a prefix
    // match: `.logR` matches its built-in `.log`, so hla_a.logR_mqc.png reports as plain
    // `hla_a` and collides with the gene's other plots. With underscores all five per gene
    // come through as distinct sections.
    """
    for pdf in ${plot_dir}/*.pdf; do
        [ -e "\$pdf" ] || continue
        name=\$(basename "\$pdf" .pdf | tr '.' '_')
        pdftoppm -png -r 150 -singlefile "\$pdf" "${somatic_name}_\${name}_mqc"
    done

    # The result table is keyed by HLA gene alone - hla_a, hla_b, hla_c - which is unique
    # within a pair and not across them. MultiQC reads the first column as the row key and
    # merges every file it matches into one table, so without the pair in the key a second
    # pair overwrites the first rather than appending to it.
    awk -F'\t' -v OFS='\t' -v pair="${somatic_name}" '
        NR == 1 { \$1 = "Pair_HLAGene"; print; next }
        { \$1 = pair "_" \$1; print }
    ' ${loh_res} > ${somatic_name}_lohres.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pdftoppm: \$(pdftoppm -v 2>&1 | head -1 | sed 's/pdftoppm version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${somatic_name}_hla_a_logR_mqc.png
    touch ${somatic_name}_lohres.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pdftoppm: 26.09.0
    END_VERSIONS
    """
}
