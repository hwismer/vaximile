# UCSF Custom Immunoprofiler - Tumor Neoantigen Vaccine Pipeline

## Prerequisites

  ### VEP Data Cache
    Download the VEP data cache that corresponds to your chosen genome.
    https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/
  ### VEP Plugins
    In order to run PVACtools, install the Frameshift and Wildtype plugins into the VEP plugins folder.
    See the PVACtools instructions here: https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/vep.html
  ### Arriba Data
    Must supply an a blacklist, known fusions, and protein domans file that can be found in the tarball here: 
    https://github.com/suhrig/arriba/releases
  ### CTAT Resource Directory
    Install the CTAT plug-and-play resource bundle for your chosen genome.
    https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
    Unzip and supply the ctat_genome_lib_dir directory to the ctat_resource_directory parameter.
  
  ### Optional
  #### BWA2 Index
    Prebuilt index files compatible with bwa2 for the chosen reference genome. If not supplied, an index will automatically be created.
  #### STAR Index
    Prebuilt index files compatible with STAR 2.7.10 for the chosen reference genome. If not supplied, an index will automatically be created.
  #### Kallisto Index
    Prebuilt index files compatible with kallisto for the chosen reference genome. If not supplied, an index will automatically be created.
  

## Install
'Simply clone the repo'

## Running

## Input

## Output

