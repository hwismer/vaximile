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
