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

### Fixed

- `mhcflow` referenced its conda environment as `./envs/mhcflow.yml`, which resolved
  against the launch directory. The file now lives beside the module and is referenced with
  `${moduleDir}`.
