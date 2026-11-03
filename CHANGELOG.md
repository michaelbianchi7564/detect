# Changelog

All notable changes to the **DETECT** pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec(v2.0.0.html).

---

## [Unreleased]

### Planned for future releases:
- **Dynamic FASTQ input Discovery:** Support for flexible sample naming suffixes (`_1`/`_2`, `_R1`/`_R2`) and configurable trimming bypass (`skip_trimming: true`).
- **Config Sync Protection:** Automated `config.yaml` replication check within `run_pipeline.sh` to prevent acciental overwrites.
- **DAG & Report Visualization:** Automatic generation of Snakemake HTML execution reports and rule graph exports (`dag.pdf`).
- **Container Support:** Docker and Singularity container support (`v1.2.0`).
- **Metrics Integration:** Picard `CollectRnaSeqMetrics` integration (`v1.2.0`).

---
## [1.1.1] 2026-10-??

### Fixed
- **featureCounts Paired-End Handling:** Removed redundant hard-coded paired-end flags to ensure clean and correct argument generation directly from the configuration file.
- **DESeq2 Multithreading:** Explicitly passed configured thread counts to the `differential_exp` execution environment, enabling parallel processing in DESeq2.

---

## [1.1.0] 2026-10-16

### Added
- **MultiQC Integration:** Added a new `multiqc` rule generating a single interactive HTML report (`MULTIQC_REPORT/multiqc_report.html`) summarizing QC and alignment metrics, included in `rule all`.
- **featureCounts Log Redirection:** Correction of shell logic to capture execution logs (`2> featureCounts.log`)

### Changed
- **License Governance:** Officially migrated repository license from GPLv3 to **AGPLv3** (`AGPL-3.0-or-later`) to close cloud execution loopholes and enforce proper commercial licensing for enterprise R&D.
- **Repository Metadata:** Updated `README.md` and `CITATION.cff` to reflect the personal IP authorship, personal contact channels, and AGPLv3/Commercial dual-licensing structure.
- **Wrapper Script:** Shared snakemake arguments are now defined once in `run_pipeline.sh`, and the Snakemake output is also logged to a temporary file for error detection.

### Fixed
- **Pipeline Execution Lock:** Resolved Snakemake execution locking behavior to ensure recovery and prevent file locks after interrupted runs.
- **Interrupted Run Recovery:** Added `--rerun-incomplete` to regenerate incomplete files after an abrupt interruption.
- **fastq Extra Parameters:** Fixed `dict_to_cli` emitting a literal `--extra` option and duplicating the `--extra` flags, which made fastp fail.

---

## [1.0.0] - 2026/07/25

### Added
- Initial release of DETECT pipeline.
- Core Snakemake rules for FASTQ quality control (`fastp`), genome alignment (`HISAT2`), read counting (`featureCounts`), and differential expression analysis (`DESeq2`).
- Automated heatmap generation and gene clustering script in R.
