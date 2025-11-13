// process.scratch = true
// process.errorStrategy = 'finish'

process {
    executor = "slurm"
    cpus = 30
    memory = "128GB"
    time = 14.d
    maxForks = 30
    //clusterOptions '--gres=scratch:500G'
}

apptainer {
    enabled = true
    autoMounts= true
}

executor {
    name = "slurm"
    queueSize = 30
}
conda.enabled = true
