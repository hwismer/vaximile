// process.scratch = true
// process.errorStrategy = 'finish'

process {
    executor = "slurm"
    cpus = 6
    memory = "128GB"
    time = 14.d
    maxForks = 30
    //clusterOptions '--gres=scratch:500G'
}

apptainer {
    enabled = true
    autoMounts= true
    //runOptions = '--writable-tmpfs -B ${projectDir}/work/:/tmp'
    runOptions = "-B ${projectDir}/work:/scratch"
    //home = '${project_dir}/work'
}

executor {
    name = "slurm"
    queueSize = 30
}
conda.enabled = true
