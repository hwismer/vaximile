#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=14-00:00:00

WORKDIR=/c4/home/hwismer/nextflow/work/

ulimit -n 1024

nextflow -C vax_config.nf run \
    -offline vax_workflow.nf \
    -params-file params.json \
    -resume \
    -ansi-log true
