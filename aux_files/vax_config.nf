
process {
    executor = "slurm"
    time = 14.d
    scratch = true
    queue = "krummellab,common"
    clusterOptions = '--gres=scratch:250G'
    errorStrategy = 'ignore'
    maxRetries = 2
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
workflow {
    output {
        mode = "copy"
    }
}
conda.enabled = true

outputDir = '.'

