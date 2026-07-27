// nextflow.config

process {
  // Reasonable defaults for most workflows
  cpus   = 1
  memory = '2 GB'
  time   = '2h'

  errorStrategy = 'ignore'
  maxRetries    = 1
}

profiles {

  // Run on the submission machine
  local {
    process.executor = 'local'
  }

  // Common HPC pattern: SLURM + Apptainer/Singularity
  slurm {
    process.executor = 'slurm'
    singularity.enabled = true
    apptainer.enabled = true   // use this instead on newer setups
    // process.queue = 'normal'   // uncomment and set your cluster default
  }

  // PBS/Torque clusters
  pbs {
    process.executor = 'pbs'
    singularity.enabled = true
  }

  // LSF clusters
  lsf {
    process.executor = 'lsf'
    singularity.enabled = true
  }

  // AWS Batch
  awsbatch {
    process.executor = 'awsbatch'
    // Usually use Docker-compatible containers in cloud
    docker.enabled = true
  }

  // Google Cloud Batch
  gcp {
    process.executor = 'google-batch'
    docker.enabled = true
    google {
      project  = 'PROJECT_ID'
      location = 'us-central1'
    }
  }
}

// Optional per-task tuning by label
process {
  withLabel: bigmem {
    memory = '128 GB'
    time   = '48h'
  }

  withLabel: gpu {
    // Example only; exact directive depends on your environment
    // accelerator = 1
  }
}












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

