# vaximile: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased]

### Changed

- Reorganised the repository into the nf-core directory layout: `main.nf` entry point,
  `workflows/vaximile/`, one directory per module under `modules/local/`, one directory per
  subworkflow under `subworkflows/local/`, and parameter defaults moved to
  `nextflow.config` and `nextflow_schema.json`. Process bodies are unchanged. See
  [docs/nf-core-migration.md](docs/nf-core-migration.md) for what remains.

- Added verified conda specs to 39 container-only modules, so `-profile conda` now covers
  70 of 94 local modules instead of 31. Containers are retained; every spec was checked
  against bioconda with `conda search` and pins the container's version (exceptions:
  `bwa-mem2=2.2.1`, inferred because the container tag does not state a version, and
  `ensembl-vep=115`, which is how bioconda publishes release 115.0). The 20 still
  container-only have real blockers — GATK3's licensed jar, unpackaged DeepSomatic and
  HLA-HD, major version gaps, and absolute container paths in scripts — tabulated in
  docs/usage.md.
- Converted all 94 local modules to nf-core process conventions. Each now carries a
  `process_*` label instead of its own `cpus`/`memory` (tiers in `conf/base.config`), 67
  support `task.ext.prefix`, and 19 read `task.ext.args` from `conf/modules.config`.
  Output filenames are unchanged. 27 modules cannot use `ext.prefix` without risking a
  glob that captures their own inputs; those are enumerated in
  [docs/nf-core-migration.md](docs/nf-core-migration.md).
- Resource requests now come from eight tiers. Each module was placed in the smallest tier
  meeting or exceeding its previous CPU and memory request, so 79 of 94 request somewhat
  more than before and none requests less.

### Fixed

- `HLAHD` declared its directory output as the glob `*/result/`. A glob ending in `/`
  never matches, so the task failed with "Missing output file(s) `*/result/`" even though
  HLA-HD had run and written its results. Introduced when the explicit
  `./${meta.sample_name}/result/` path became a glob during the ext.prefix conversion. Now
  `path("*/result", type: 'dir')`. No other module has a trailing-slash glob.
- `SALMON_QUANT` passed `--libType` after `-1/-2`, which salmon rejects outright
  ("The (--libType/-l) option must precede the input files"). Introduced when
  `--libType A` and `--validateMappings` — which sat on opposite sides of the read
  files — were collapsed into a single `$args` placed after them. `$args` now precedes
  `-1/-2`. The other 18 modules using `ext.args` were audited for the same reordering;
  10 also had non-contiguous flags but all are order-insensitive option-only CLIs
  (VEP, GATK, bcftools, DeepVariant/DeepSomatic) or place `$args` after every
  positional (pVACseq, pVACfuse).
- Disabled `timeline`, `report` and `trace` in `nextflow.config`. Enabling them makes
  Nextflow inject `command -v ps || exit 1` into every task wrapper, which killed tasks
  running in containers without `procps` (STAR, DeepVariant, DeepSomatic, bcftools) before
  their tool ran. The symptom is an empty `.command.out` and only the "Command 'ps'
  required by nextflow" line in `.command.err`. Request the reports per-run with
  `-with-report`/`-with-timeline`/`-with-trace` instead. `dag` is unaffected and stays on.
- `nextflow.config` interpolated `${manifest.name}` inside the `validation.help.command`
  string, where `manifest` is not in scope. This failed config parsing outright, so every
  `nextflow` invocation from the repository root aborted before compiling anything.

- `mhcflow` referenced its conda environment as `./envs/mhcflow.yml`, which resolved
  against the launch directory. The file now lives beside the module and is referenced with
  `${moduleDir}`.
