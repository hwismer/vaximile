process MULTIQC {

    // Per-patient MultiQC report, including custom HLA and copy-number sections and how the run was launched.

    label 'process_low'
    conda "bioconda::multiqc=1.34-0"

    tag "Running MultiQC on ${patient}"

    input:
        tuple val(patient), val(somatic_names), val(sample_names), path(files)
        val(run_info)

    output:
        tuple val(patient), path("*_report.html"), emit: html

    script:
    
    def prefix = task.ext.prefix ?: "${patient}"
    def clean_sample_names = sample_names.findAll { it != null }.unique().sort { -it.size() }

    def rename_tsv = (clean_sample_names)
        .collect { s -> "^${s}.*\t${s}" }
        .join('\n')
    def run_info_json = groovy.json.JsonOutput.toJson(run_info)

    """
    printf '%s\n' "${rename_tsv}" > ${patient}_rename.tsv

    # Run Information section: how the pipeline was launched, and every parameter it used.
    cat > run_info.json <<'JSON'
${run_info_json}
JSON
    python3 - <<'PY'
import html
import json

info = json.load(open("run_info.json"))


def table(entries, code_keys):
    rows = "".join(
        f"<tr><th>{html.escape(k)}</th><td>{'<code>' + html.escape(v) + '</code>' if k in code_keys else html.escape(v)}</td></tr>"
        for k, v in entries.items()
    )
    return f'<table class="table table-condensed">{rows}</table>'


with open("run_info_mqc.html", "w") as fh:
    fh.write('<!--\\nid: "run_info"\\nsection_name: "Run Information"\\ndescription: "How this run of the pipeline was launched, and every parameter it used, whether set on the command line or left at its default."\\n-->\\n')
    fh.write("<h4>Run</h4>" + table(info["details"], {"Command"}))
    fh.write("<h4>Parameters</h4>" + table({f"--{k}": v for k, v in info["params"].items()}, set()) + "\\n")
PY

	# Optitype custom content
    for f in *_result.tsv; do
        [ -e "\$f" ] || continue
        s=\$(basename "\$f" _result.tsv)
        awk -F'\t' -v OFS='\t' -v s="\$s" 'NR == 1 { \$1 = "Sample"; print; next } { \$1 = s; print }' "\$f" > "\${s}_optitype.tsv"
    done
   cat > multiqc_config.yaml <<EOF
exclude_modules:
  - optitype

mosdepth_config:
  general_stats_coverage: [20, 50, 100]

# Draw the variant scatter plots straight away instead of behind a "Show plot" button.
# MultiQC 1.34 still applies the deprecated num_datasets_plot_limit as well.
plots_defer_loading_numseries: 2000
num_datasets_plot_limit: 2000

custom_content:
  order:
    - run_info
    - somatic
    - germline
    - hla
    - cnv
    - neoantigens

custom_data:
  hla_calls:
    file_format: tsv
    section_name: "HLA-HD"
    description: "HLA class I and II allele calls from HLA-HD, per library"
    parent_id: hla
    parent_name: "HLA"
    parent_description: "HLA typing, and loss of heterozygosity at the HLA locus in each tumour against its own normal"
    plot_type: table
    pconfig:
      id: "hla_calls"
      title: "HLA-HD Calls"

  optitype_calls:
    file_format: tsv
    section_name: "OptiType"
    description: "HLA class I allele calls from OptiType, per library"
    parent_id: hla
    parent_name: "HLA"
    plot_type: table
    pconfig:
      id: "optitype_calls"
      title: "OptiType Calls"

  hla_loh:
    file_format: tsv
    section_name: "LOH metrics"
    description: "Allele-level copy number and loss-of-heterozygosity statistics from lohhlamod, one row per HLA gene per tumour/normal pair. The remaining columns - copy number bounds, the four median logR columns, bin counts and the per-allele loss percentages - are hidden by default and can be shown from Configure Columns."
    parent_id: hla
    parent_name: "HLA"
    plot_type: table
    pconfig:
      id: "hla_loh"
      title: "HLA LOH"
    headers:
      HLA_A1:
        title: "Allele 1"
      HLA_A2:
        title: "Allele 2"
      HLA_A1_CN:
        title: "A1 CN"
        format: "{:,.2f}"
      HLA_A2_CN:
        title: "A2 CN"
        format: "{:,.2f}"
      MM_LogR_Paired_Pvalue:
        title: "MM logR p"
        format: "{:,.2e}"
      Median_BAF:
        title: "Median BAF"
        format: "{:,.3f}"
      Pct_CN_Diff_Supporting_Bins:
        title: "% bins CN diff"
        format: "{:,.1f}"
      HLA_A1_CN_Lower: { hidden: True }
      HLA_A1_CN_Upper: { hidden: True }
      HLA_A2_CN_Lower: { hidden: True }
      HLA_A2_CN_Upper: { hidden: True }
      HLA_A1_Median_LogR: { hidden: True }
      HLA_A2_Median_LogR: { hidden: True }
      HLA_A1_MM_Median_LogR: { hidden: True }
      HLA_A2_MM_Median_LogR: { hidden: True }
      Num_MM: { hidden: True }
      Num_Bins: { hidden: True }
      Num_MM_Bins: { hidden: True }
      Pct_A1_Loss_Supporting_Bins: { hidden: True }
      Pct_A2_Loss_Supporting_Bins: { hidden: True }

  somatic_vaf:
    section_name: "VAF"
    description: "Histogram of tumor DNA VAF of the final somatic calls, after the two-caller consensus and the population filter, one panel per tumor/normal pair. A clonal peak sits near half the tumor purity; a pile-up at very low VAF is mostly calls at the limit of detection."
    parent_id: somatic
    parent_name: "Somatic"
    plot_type: image

  somatic_normal_vaf:
    section_name: "VAF - tumor vs normal"
    parent_id: somatic
    parent_name: "Somatic"
    plot_type: image

  germline_vaf:
    section_name: "VAF"
    description: "Histogram of the VAF of the consensus germline calls, one panel per normal. Heterozygous variants should peak at 0.5 and homozygous ones at 1. A heterozygous peak that is broad or shifted away from 0.5 suggests contamination or a mixed sample."
    parent_id: germline
    parent_name: "Germline"
    plot_type: image

  somatic_rna_vaf:
    section_name: "VAF - RNA vs DNA"
    description: "Tumor DNA against RNA VAF of the final somatic calls that have RNA coverage, one panel per tumor/normal pair. Mutations near the diagonal are expressed in proportion to their DNA fraction; those on the x-axis are not expressed."
    parent_id: somatic
    parent_name: "Somatic"
    plot_type: image

  ascat_metrics:
    file_format: tsv
    section_name: "ASCAT metrics"
    description: "Tumour purity, ploidy and ASCAT's own QC metrics, one row per tumour/normal pair"
    parent_id: cnv
    parent_name: "Copy Number"
    parent_description: "Allele-specific copy number from ASCAT, whose purity and ploidy the HLA LOH analysis uses"
    plot_type: table
    pconfig:
      id: "ascat_metrics"
      title: "ASCAT Metrics"
EOF

    # One section per pair for the tiled sheets, written here rather than declared above
    # because the pairs are only known from what was staged. A single custom_data entry
    # matching every sheet would read better in the config and lose data in the report:
    # MultiQC keeps one image per section, so the second pair's sheet would replace the
    # first. Naming the sections here also keeps the pair name intact, where letting
    # MultiQC derive it from the file name turns TP_M1_N1 into "TP M1 N1".
    for sheet in *_hla_loh_mqc.png; do
        [ -e "\$sheet" ] || continue
        pair=\$(basename "\$sheet" _hla_loh_mqc.png)
        printf '  hla_loh_plots_%s:\n    section_name: "LOH plots - %s"\n    parent_id: hla\n    parent_name: "HLA"\n    plot_type: image\n' "\$pair" "\$pair" >> multiqc_config.yaml
    done
    for sheet in *_ascat_mqc.png; do
        [ -e "\$sheet" ] || continue
        pair=\$(basename "\$sheet" _ascat_mqc.png)
        printf '  ascat_plots_%s:\n    section_name: "ASCAT plots - %s"\n    parent_id: cnv\n    parent_name: "Copy Number"\n    plot_type: image\n' "\$pair" "\$pair" >> multiqc_config.yaml
    done

    cat >> multiqc_config.yaml <<EOF
sp:
  hla_calls:
    fn: "*_hlahd.tsv"
  optitype_calls:
    fn: "*_optitype.tsv"
  hla_loh:
    fn: "*_lohres.tsv"
  ascat_metrics:
    fn: "*_ascatmetrics.tsv"
  somatic_vaf:
    fn: "*_somatic_vaf_mqc.png"
  somatic_normal_vaf:
    fn: "*_somatic_normal_vaf_mqc.png"
  germline_vaf:
    fn: "*_germline_vaf_mqc.png"
  somatic_rna_vaf:
    fn: "*_somatic_rna_vaf_mqc.png"
EOF

    for sheet in *_hla_loh_mqc.png; do
        [ -e "\$sheet" ] || continue
        pair=\$(basename "\$sheet" _hla_loh_mqc.png)
        printf '  hla_loh_plots_%s:\n    fn: "%s"\n' "\$pair" "\$sheet" >> multiqc_config.yaml
    done
    for sheet in *_ascat_mqc.png; do
        [ -e "\$sheet" ] || continue
        pair=\$(basename "\$sheet" _ascat_mqc.png)
        printf '  ascat_plots_%s:\n    fn: "%s"\n' "\$pair" "\$sheet" >> multiqc_config.yaml
    done

    printf '\nsample_names_replace_regex: true\n' >> multiqc_config.yaml
    
    cat multiqc_config.yaml
    
    multiqc \
        -n ${prefix}_report.html \
        --replace-names ${patient}_rename.tsv \
        -c multiqc_config.yaml \
        -i "${patient} - Vaximile Pipeline Metrics" \
        .
    """

    stub:
    def prefix = task.ext.prefix ?: "${patient}"
    """
    touch ${prefix}_report.html
    """
}
