
process {
    executor = "slurm"
    time = 14.d
    //cache = "lenient"
    scratch = true
    queue = "krummellab"
    clusterOptions = '--gres=scratch:250G --nodelist=c4-n18,c4-n35 --exclude=c4-n1,c4-n2,c4-n4'
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
conda.enabled = true

outputDir = '.'

