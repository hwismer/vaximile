#!/bin/bash

nextflow -C ./aux_files/vax_config.nf run \
    vax_workflow.nf \
    -params-file ./aux_files/params.json \
    -resume \
    -ansi-log true
