process MULTIQC {

    /*
    Runs MultiQC on gathered files on a per-patient basis.
    Currently also imports HLA-HD calls into table format.

    Replaces sample names with sample names from metadat and merged samples that start with Merge
    */

    label 'process_low'
    conda "bioconda::multiqc=1.34-0"

    tag "Running MultiQC on ${patient}"

    input:
        tuple val(patient), val(somatic_names), val(sample_names), path(files)

    output:
        tuple val(patient), path("*_report.html")

    script:
    
    def prefix = task.ext.prefix ?: "${patient}"
    def clean_sample_names = sample_names.findAll { it != null }.unique().sort { -it.size() }
    def clean_somatic_names = somatic_names.findAll { it != null }.unique().sort { -it.size() }

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

sp:
  hla_calls:
    fn: "*_hlahd.tsv"

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${patient}"
    """
    touch ${prefix}_report.html
    """
}
