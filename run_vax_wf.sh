#!/bin/bash

nextflow -C vax_config.nf run \
    -offline vax_workflow.nf \
    -params-file params.json \
    -resume \
    -ansi-log true
