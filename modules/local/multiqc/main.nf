process MULTIQC {

    /*
    Runs MultiQC on gathered files on a per-patient basis.
    Currently also imports HLA-HD calls into table format.

    Replaces sample names with sample names from metadata and merged samples that start with Merge
    */

    label 'process_low'
    conda "bioconda::multiqc=1.34-0"

    tag "Running MultiQC on ${patient}"

    input:
        tuple val(patient), val(somatic_names), val(sample_names), path(files)

    output:
        tuple val(patient), path("*_report.html"), emit: html
        path "versions.yml", topic: versions

    script:
    
    def prefix = task.ext.prefix ?: "${patient}"
    def clean_sample_names = sample_names.findAll { it != null }.unique().sort { -it.size() }

    def rename_tsv = (clean_sample_names)
        .collect { s -> "^${s}.*\t${s}" }
        .join('\n')

    """
    printf '%s\n' "${rename_tsv}" > ${patient}_rename.tsv
   
   cat > multiqc_config.yaml <<EOF
custom_data:
  hla_calls:
    file_format: tsv
    section_name: "HLA-HD"
    description: "HLA Allele Calls"
    plot_type: table
    pconfig:
      id: "hla_calls"
      title: "HLA-HD Calls"

  hla_loh:
    file_format: tsv
    section_name: "HLA LOH"
    description: "Allele-level copy number and loss-of-heterozygosity statistics from lohhlamod, one row per HLA gene per tumour/normal pair. The remaining columns - copy number bounds, the four median logR columns, bin counts and the per-allele loss percentages - are hidden by default and can be shown from Configure Columns."
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

sp:
  hla_calls:
    fn: "*_hlahd.tsv"
  hla_loh:
    fn: "*_lohres.tsv"

sample_names_replace_regex: true
EOF
    
    cat multiqc_config.yaml
    
    multiqc \
        -n ${prefix}_report.html \
        --replace-names ${patient}_rename.tsv \
        -c multiqc_config.yaml \
        -i "${patient} - UCSF Custom Immunoprofiler CustomVax Pipeline Metrics" \
		-b "Info | Patient: ${patient} \n | VEP Outputs: Germline (normal sample name) and Somatic (tumor/normal pair, e.g. Patient1_T1_N1)" \
        .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${patient}"
    """
    touch ${prefix}_report.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.34
    END_VERSIONS
    """
}
