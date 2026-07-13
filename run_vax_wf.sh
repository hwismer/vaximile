#!/bin/bash

nextflow -C ./inputs/vax_config.nf run \
    vax_workflow.nf \
    -params-file ./inputs/params.json \
    -resume \
    -ansi-log true
