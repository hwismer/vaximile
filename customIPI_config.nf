// process.scratch = true
// process.errorStrategy = 'finish'

process {
    executor = "slurm"
    cpus = 32
    memory = "64GB"
    time = 14.d
    //clusterOptions '--gres=scratch:500G'
}

apptainer {
    enabled = true
    autoMounts= true
    runOptions = '--writable-tmpfs -B ${projectDir}/work/:/tmp'
    runOptions = "-B ${projectDir}/work:/scratch"
    //home = '${project_dir}/work'
    pullTimeout="30min"
}

executor {
    name = "slurm"
    queueSize = 30
}

plugins {
    id 'nf-google'
}
conda.enabled = true
