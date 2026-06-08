# UCSF Custom Immunoprofiler - Tumor Neoantigen Vaccine Pipeline

## Prerequisites

  ### VEP Data Cache
    Download the VEP data cache that corresponds to your chosen genome.
    https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/
  
  ### Optional
  #### BWA2 Index
    Prebuilt index files compatible with bwa2 for the chosen reference genome. If not supplied, an index will automatically be created.
  #### STAR Index
    Prebuilt index files compatible with STAR 2.7.10 for the chosen reference genome. If not supplied, an index will automatically be created.
  #### Kallisto Index
    Prebuilt index files compatible with kallisto for the chosen reference genome. If not supplied, an index will automatically be created.
  

## Install
    1. Clone the repository
    2. Install Nextflow,
    3. Prepare required input files and optional parameters
    4. Prepare capture kit samplesheet
    4. Prepare patient samplesheet

## Running

nextflow -C ./aux_files/vax_config.nf run \
    vax_workflow.nf \
    -params-file ./aux_files/params.json \
    -resume \
    -ansi-log true

## Input

## Output

