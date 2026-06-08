process FASTP {

    cpus 8
    memory "16GB"
    conda "bioconda::fastp=1.0.1"

    tag "FastP on ${meta.sample_name} w/ ${meta.molecule}"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_R1_fastp.fastq.gz"), path("${meta.sample_name}_${meta.molecule}_R2_fastp.fastq.gz"), emit: fastqs
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}.json"), emit: reports
    

    script:
        """
        fastp --thread $task.cpus \
              -i ${reads[0]} \
              -I ${reads[1]} \
              -o "${meta.sample_name}_${meta.molecule}_R1_fastp.fastq.gz" \
              -O "${meta.sample_name}_${meta.molecule}_R2_fastp.fastq.gz" \
              -R "${meta.sample_name}_${meta.molecule}_fastp_report" \
              -h "${meta.sample_name}_${meta.molecule}.html" \
              -j "${meta.sample_name}_${meta.molecule}.json" \

        """
}

process SAMTOOLS_FLAGSTAT {

    cpus 8
    memory "24GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"
    
    tag "Samtools flagstat on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("${meta.sample_name}_${meta.molecule}.flagstat")

    script:
    """
    samtools flagstat --threads $task.cpus $bam > ${meta.sample_name}_${meta.molecule}.flagstat
    """

}

process SAMTOOLS_COVERAGE {
    
    cpus 8
    memory "24GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"
    
    tag "Samtools coverage on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("${meta.sample_name}_${meta.molecule}_coverage.tsv")

    script:
    """
    samtools coverage $bam > ${meta.sample_name}_${meta.molecule}_coverage.tsv
    """
}


process SAMTOOLS_IDXSTATS {
    
    cpus 4
    memory "16GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    tag "Samtools idxstats on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("${meta.sample_name}_${meta.molecule}_idxstats.tsv")

    script:
    """
    samtools idxstats $bam > ${meta.sample_name}_${meta.molecule}_idxstats.tsv
    """
}


process SOMALIER_EXTRACT {


    cpus 4
    memory "32GB"
    conda "bioconda::somalier=0.3.2-0 bioconda::htslib=1.23.1"

    tag "somalier extract on ${meta.sample_name} at ${sites}"

    input:
        tuple val(meta), path(bam), path(bai)
        tuple path(reference_fa), path(fai)
        path(sites)

    output:
        tuple val(meta), path("${meta.sample_name}.somalier")
        

    script:
    """
    somalier extract \
        -s $sites \
        -f $reference_fa \
        $bam 
    """
}

process SOMALIER_RELATE {


    cpus 4
    memory "32GB"
    conda "bioconda::somalier=0.3.2-0 bioconda::htslib=1.23.1"

    tag "somalier relate on patient ${patient}"

    input:
        tuple val(patient), val(metas),  path(files)

    output:
        tuple val(patient), path("*.pairs.tsv"),  emit: pairs
        tuple val(patient), path("*.samples.tsv"), emit: samples
        tuple val(patient), path("*.groups.tsv"), emit: groups
        tuple val(patient), path("*.html"), emit:html

    script:
    """
    export SOMALIER_REPORT_ALL_PAIRS=1
    somalier relate \
        -o $patient \
        *.somalier
    """
}

process MULTIQC {

    /*
    Runs MultiQC on gathered files on a per-patient basis.
    Currently also imports HLA-HD calls into table format.

    Replaces sample names with sample names from metadat and merged samples that start with Merge
    */

    cpus 2
    memory "16GB"
    conda "bioconda::multiqc=1.34-0"

    tag "Running MultiQC on ${patient}"

    input:
        tuple val(patient), val(somatic_names), val(sample_names), path(files)

    output:
        tuple val(patient), path("${patient}_report.html")

    script:
    
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
        -n ${patient}_report.html \
        --replace-names ${patient}_rename.tsv \
        -c multiqc_config.yaml \
        -i "${patient} - UCSF Custom Immunoprofiler CustomVax Pipeline Metrics" \
		-b "Info | Patient: ${patient} \n | VEP Outputs: Germline (normal sample name) and Somatic (tumor/normal pair, e.g. Patient1_T1_N1)" \
        .
    """
}
