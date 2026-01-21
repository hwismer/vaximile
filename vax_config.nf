
process {
    executor = "slurm"
    cpus = 32
    memory = "64GB"
    time = 14.d
    cache = "lenient"
    queue = "freecycle,krummellab,common"
    scratch = true
    clusterOptions = '--gres=scratch:750G'
}

apptainer {
    enabled = true
    autoMounts= true
    runOptions = '--writable-tmpfs --cleanenv -e --no-home --env PYTHONNOUSERSITE=1 -B $SINGULARITY_TMPDIR:/tmp -B $SINGULARITY_TMPDIR:/scratch'
    pullTimeout="30min"
}

executor {
    name = "slurm"
    queueSize = 50
    queue = "freecyle,krummellab,common"
}

plugins {
    id 'nf-google'
}
conda.enabled = true
