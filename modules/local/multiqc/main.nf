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
    description: "Allele-level copy number and loss-of-heterozygosity statistics from lohhlamod, one row per HLA gene per tumour/normal pair"
    plot_type: table
    pconfig:
      id: "hla_loh"
      title: "HLA LOH"

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
