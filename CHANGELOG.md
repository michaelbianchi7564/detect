# Changelog

All notable changes to the **DETECT** pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec(v2.0.0.html).

---

## [Unreleased]

### Planned for v1.1.0
- **MultiQC Integration:** MultiQC rule to generate a single interactive HTML report (`multiqc_report.html`) summarizing QC and alignment metrics.
- **Dynamic FASTQ input Discovery:** Support for flexible sample naming suffixes (`_1`/`_2`, `_R1`/`_R2`) and configurable trimming bypass (`skip_trimming: true`).
- **Config Sync Protection:** Automated `config.yaml` replication check within `run_pipeline.sh` to prevent acciental overwrites.
- **DAG & Report Visualization:** Automatic generation of Snakemake HTML execution reports and rule graph exports (`dag.pdf`)

## [1.1.0] 2026-10-16

### Added
- **featureCounts Log Redirection:** Correction of shell logic to capture execution logs (`2> featureCounts.log`)

### Changed
- **License Governance:** Officially migrated repository license from GPLv3 to **AGPLv3** (`AGPL-3.0-or-later`) to close cloud execution loopholes and enforce proper commercial licensing for enterprise R&D.
- **Repository Metadata:** Updated `README.md` and `CITATION.cff` to reflect the personal IP authorship, personal contact channels, and AGPLv3/Commercial dual-licensing structure.

### Fixed
- **Log Direction:** Corrected `featureCounts` shell execution logic to properly capture execution logs (`2> featureCounts.log`)

---

## [1.0.0] - 2026/07/25

### Added
- Initial release of DETECT pipeline.
- Core Snakemake rules for FASTQ quality control (`fastp`), genome alignment (`HISAT2`), read counting (`featureCounts`), and differential expression analysis (`DESeq2`).
- Automated heatmap generation and gene clustering script in R.
