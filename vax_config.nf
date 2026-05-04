
process {
    executor = "slurm"
    cpus = 32
    memory = "64GB"
    time = 14.d
    cache = "lenient"
    scratch = true
    clusterOptions = '--gres=scratch:500G'
    errorStrategy = 'ignore'
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
}

plugins {
    id 'nf-google'
}
google {
    enabled = true
}
conda.enabled = true

outputDir = 'results'
