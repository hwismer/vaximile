// process.scratch = true
// process.errorStrategy = 'finish'

workDir = '/c4/home/hwismer/nextflow/work/'
cacheDir = '/c4/home/hwismer/nextflow/cache/'

process {
    executor = "slurm"
    cpus = 32
    memory = "64GB"
    time = 14.d
    cache = "lenient"
    queue = "freecycle,krummellab,common"
    //clusterOptions = '--gres=scratch:500G'
}

apptainer {
    enabled = true
    autoMounts= true
    runOptions = '--writable-tmpfs'
    runOptions = '-B $SINGULARITY_TMPDIR:/tmp -B $SINGULARITY_TMPDIR:/scratch'
    //envWhitelist = ['SINGULARITY_TMPDIR']
    runOptions = '-B ${workDir}:/tmp/'
    runOptions = "-B ${workDir}:/scratch"
    //home = '${project_dir}/work'
    pullTimeout="30min"
}

executor {
    name = "slurm"
    queueSize = 30
    queue = "freecyle,krummellab,common"
}

plugins {
    id 'nf-google'
}
conda.enabled = true
