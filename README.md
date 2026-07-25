<div align="center">

# 🧬🐍 DETECT 🐍🧬

### **D**ifferential **E**xpression **T**ranscriptomic **E**valuation and **C**lustering **T**ool

[![Workflow: Snakemake](https://img.shields.io/badge/Workflow-Snakemake-emerald?style=flat-square)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/Version-v1.0.0-green?style=flat-square)](https://github.com/michaelbianchi7564/detect)
---

A scalable, reproducible **Snakemake pipeline** for genome-referenced differential expression analysis, unsupervised clustering, and automated heatmap generation.

</div>

---
## 1. 🧬 OVERVIEW

This repository contains a fully automated, end-to-end, and highly reproducible **RNA-seq pipeline** built using **Snakemake**, **Conda**, and **R**.

Designed following the FAIR data principles (Findable, Accessible, Interoperable and Reusable), the pipeline clearly separates workflow logic and user data. Researchers can shift from raw sequencing reads to statistical differential expression, unsupervised gene clustering, and publication-ready visualizations without modifying the underlying codebase.

### Core Phylosophy & Architecture
* **Zero-Code Execution:** Users interact exclusively with localized configuration files and data directories. In this manner the pipeline logic remains unmodified and standardized.
* **Isolated Environments:** Software dependencies are managed automatically through isolated Conda environments and R package versions, allowing the pipeline to run reproducibly.
* **On-the-fly Data Validation:** the R `cli` package, the pipeline performs rigorous structural checks on metadata sample naming consistency, and count matrix integrity, aborting early with informative alerts if blocking anomalies are detected.
* **Run Autonomy**: The pipeline architecture permits complete isolation between different projects. It needs a dedicated subdirectory for each run, independent datasets can be processed using distinct configuration profiles without risk of overwriting previous results.
* **Gene-Level Quantitation & Integrity:** This pipeline's version is strictly at the **gene-level**. By passing unique feature identifiers from `featureCounts`.

### Pipeline Workflow & Integrated Software
This workflow integrates widely used bioinformatics tools into a single reproducible pipeline, from raw sequencing reads to differential expression, gene clustering and publication-ready results:

`FASTQ` → `FastQC` / `fastp` → `HISAT2` → `Samtools` → `featureCounts` → `DESeq2` (R) → Clustering & Heatmaps

1. **Quality Control & Trimming:**
* `FastQC`for raw read quality control.
* `fastp`for ultra-fast, automated adapter trimming and quality filtering.
2. **Alignment & Quantification:**
* `HISAT2`for splice-aware genomic alignment.
* `samtools`& `Picard`for structural BAM sorting, indexing, and alignment metrix extraction.
* `featureCounts`(Subread) for rapid genomic geature-level quantification.
3. **Differential Expression & Clustering (R):**
* `DESeq2`for robust differential expression analysis, utilizing `ashr`for log2 fold-change shrinkage.
* `ggplot2` & `pheatmap`for automated generation of PCA plots, Volcano plots, and expression heatmaps.
* `dynamicTreecut`for unsupervised hierarchical clustering and automated gene cluster assignment.

### Input & Output at first glance
* **Inputs:** Raw paired-end or single-end `FASTQ` files, a sample `metadata.csv`, a `config.yaml`, and reference genomes files (`FASTA` / `GTF` / `GFF`).
* **Outputs:** Results in `XLSX`and `CSV`tables for differential expressions, detailed `PDF` plots for visual inspection, and R objects (`.rds`) to allow immediate downstream customization in R.

## 2. 📋 REQUIREMENTS

Before installing and running the pipeline, ensure your system meets the following prerequisites:

### Operating System & Compatibility
The underlying bioinformatic binaries (such as `HISAT2`, `SamTools`, and `Subread`) are natively compiled for Linux architectures.
* **Linux (Highly Recommended):** Fully tested and natively supported.
* **Windows:** Not natively supported. Windows users can run the pypeline perfectly by enabling **WSL2** (Windows Subsystem for Linux) and installing Miniconda inside the Linux subsystem.
* **MacOS:** Not natively supported on Apple silicon (M1/M2/M3/M4) due to missing ARM64 binaries for legacy bioinformatic tools on Bioconda. Intel-based Macs may work.

### Hardware Resource Allocations
* **CPU Cores:** The execution script provides dynamic core allocation. By default, it allocates **4 cores** if no parameter is provided, but it can scale up to match your hardware capacity.
* **Memory (RAM):** A minimum of 16GB RAM is recommended. Depending on the size of your reference genome and data, 32GB or more is highly preferred to handle the alignment index (`HISAT2`) efficiently.

### Software Requirements
* **Miniconda:** A working installation of Miniconda (or Anaconda) is required to manage and build the workflow's isolated software environments.

## 3. 🧬 REFERENCE GENOME PREPARATION

Before launching a new run, you must manually place your reference sequence and annotation files into the  folder `reference_genomes/`. To ensure a correct integration with the alignment and quantification rules, follow this mandatory naming convention:
1. **Reference Genome:** Save the genomic assembly file in FASTA format as `genome.fna`.
2. **Annotation File:** Save the corresponding structural gene annotation in GFF format as `annotation.gff`.

Your project tree for references must look exactly like this:

reference_genomes/
├── genome.fna              # Reference genome
├── annotation.gff
└── annotation.gtf

## 4. 🚀 QUICK START

To give both quick environment verification and full biological validation, the repository provides two distinct testing tiers located in the `tests/`directory.

### Option A: 2-Minutes Quick-Start (Recommended)
To verify that all pipeline dependencies, paths, and Snakemake rules are correctly configured, you can run a mini-test using a pre-subsampled *E. coli* dataset.

Step 1. **Prepare the test dataset:**
```bash
chmod +x tests/download_quickstart_data.sh
./tests/download_quickstart_data.sh
```

Step 2: **Launch the quick-test run:**
`bash ./run_pipeline.sh E_coli_test 4`

### Option B: Full Biological Benchmark and Validation

This pipeline has been rigorously validated using part of a public *Escherichia coli* RNA-seq dataset from the NCBI Short Read Archive (SRA), across the following accessions: `SRR35290504`, `SRR35290506`, `SRR35290507`, `SRR35290509`, `SRR35290510`, and `SRR35290512`.

> **Prerequisites:** Replicating this full benchmark requires `sra-tools` (providing `fasterq-dump`) and `gffread` to be pre-installed and available in your active environment:
> ```bash
> conda activate snakemake
> conda install -c bioconda -c conda-forge sra-tools gffread -y
> ```

To download the raw data and run the complete biological analysis (approx 2-3 hours total execution time, depending on your hardware specs):

### Step 1: Download the Test Dataset
Run the automated download script. It pulls the complete reference files, build the `.gtf` annotation, and fetch the raw dataset from the NCBI Cloud Storage mirror directly into `runs/E_coli_benchmark/DATA/`:
```bash
chmod +x download_benchmark_data.sh
./download_full_benchmark_data.sh
```

### Step 2: Launch the Benchmark Run
Once the FASTQ files are downloaded and compressed into the data directory, launch the pipeline wrapper script. The script automatically locates the configuration and metadata profiles dedicated to the `E_coli_benchmark` dataset:

```bash
# Run the test pipeline specifying 8 CPU cores
./run_pipeline.sh E_coli_benchmark 8
```

### Step 3: Inspect the results
The execution will take approximately 2 or 3 hours to complete, depending on your hardware specs and the reference genome build process. Once finished, you can inspect the generated quality reports, differential expression tables (`.xlsx` and `.csv`), and visualization plots (`.pdf`) in the output directory:

```bash
ls -lh runs/E_coli_benchmark/results/
```

## 5. 📁  REPOSITORY ARCHITECTURE

The structure of this repository is projected to maintain a strict separation between pipeline logic, configuration parameters, and run-specific execution data. This layout allows the pipeline to be portable and easy adaptable to new datasets.

.
├── config/
│   └── config.yaml                 # Global pipeline parameters and software settings
├── reference_genomes/              # Contains reference genomes
│   ├── your_organism_name/         # Each folder contains reference genome and annotation files of a species
├── runs/                           # Isolated run environments
│   ├── E_coli_test/                # Standard Quick-Start run directory
│   │   ├── DATA/                   # Input FASTQ files (subsampled)
│   │   ├── results/                # Pipeline's outputs
│   │   └── metadata.tsv            # Experimental design & sample mapping sheet
│   └── E_coli_benchmark/           # Full biological benchmark run directory   
│   │   ├── DATA/                   # Input FASTQ files (subsampled)
│   │   ├── results/                # Pipeline's outputs
│   │   └── metadata.tsv            # Experimental design & sample mapping sheet
├── scripts/                        # Helper scripts for downstream analyses
│   ├── Differential_Expression_Analysis.R  # DESeq2 integration & plotting script
│   └── utils.py                    # Python utility functions
├── tests/                          # Validation and setup scripts
│   ├── download_quickstart_data.sh      # Stages the ~5MB Quick-Start dataset
│   └── download_full_benchmark_data.sh  # Stages the full biological SRA dataset
├── workflow/                       # Core Snakemake pipeline logic
│   ├── envs/                       # Isolated Conda environment definition files
│   │   ├── fastp.yaml              # Quality filtering (fastp) env
│   │   ├── hisat2.yaml             # Alignment (HISAT2) env
│   │   ├── R.yaml                  # R/DESeq2 env
│   │   └── subread.yaml            # Quantitation (featureCounts) env
│   └── RNA_seq.smk                 # Master Snakefile defining the workflow rules
├── run_pipeline.sh                 # Master execution wrapper bash script
└── README.md                       # Project documentation and quick-start guide

### Architecture Highlights
* `workflow/RNA_seq.smk`: The core of the pipeline. It defines the logical flow of execution from raw FASTQ pre-processing to gene expression counts.
* `workflow/envs`: To guarantee absolute reproducibility, each rule of this pipeline is run inside an isolated Conda environment defineby by these `.yaml`files.
* `config/config.yaml`: Contains global configuration (e.g. reference genome, software flags, and tool parameters). It is designed to avoid harcoding of settings inside the rules.
* `runs/`: Every test, benchmark, or project runs inside its own subdirectory under `runs/`. This prevents data clashing and guarantee that root repository remains clean. Each run directory carries its own `DATA/`, `REF/`, and `metadata.tsv` files. 

## 6. ⚙️ CONFIGURATION & METADATA

To run the pipeline on your biological dataset, you only need to configure two files within your run directory (e.g. `runs/your_experiment/`): the sample metadata file (`metadata.tsv`) and the run-specific configuration file (`config.yaml`).

### A. The Metadata File (`metadata.tsv`)
This pipeline uses a tab-separated values (`tsv`) file to map raw FASTQ files to their respective biological conditions. This file must be placed directly in your run folder.

#### Format Requirements:
* **Format:** Tab-separated (commas or semicolons are not accepted).
* **Required Columns:**
    * `SAMPLE`: The prefix of your FASTQ files, representing their ID.
    * `CONDITION`: The experimental condition.

#### Example `metadata.tsv`:
```
SAMPLE  CONDITION
SRR35290504	TREAT_T
SRR35290506	TREAT_T
SRR35290507	TREAT_E
SRR35290509	TREAT_E
SRR35290510	TREAT_D
SRR35290512	TREAT_D
```

### B. The Configuration Parameters (`config.yaml`)
The pipeline behavior is controlled through a two-level system. Default parameters are located in `config/config.yaml`, while run-specific parameters can be customized in `runs/your_run/config.yaml`.
Instead of editing the global pipeline defaults in `config/config.yaml`, copy the template to your run directory and modify only the necessary.

```
bash

cp config/config.yaml runs/your_run/config.yaml
```

Open your local `runs/your_run/config.yaml` and update the essential paths and reference settings:

```
yaml

input_dir: "DATA"              # Folder containing your FASTQ files inside the run directory
reference_condition: "TREAT_D" # The baseline/control group for DESeq2 contrasts

reference:
  genome: "REF/your_genome.fna"
  annotation_tab: "REF/your_annotation.gtf"
```
## 7. ⚡ PIPELINE EXECUTION MODE

The execution of the pipeline is standardized through the wrapper script `run_pipeline.sh`. This script manages environment initialization and automated parameter passing to Snakemake.

### A. Dry-Run
To validate your updated config files, paths, and sample naming conventions without triggering heavy computation, perform a dry-run:
```bash
# General syntax: ./run_pipeline.sh <run_name> <cores> --dry-run
./run_pipeline.sh E_coli_benchmark 4 --dry-run
```

### B. Full Pipeline Run
To launch the actual computation:
```bash
./run_pipeline.sh <run_name> <cores>
```

## 8. 📊 DETAILED OUTPUT ARCHITECTURE

At the end of the execution, all processed files are structured inside your dedicated run directory.

runs/your_run_name/results/
├── QC/
│   └── *_fastqc.{html,zip}          # Quality control reports for raw reads
├── TRIMMING/
│   ├── *.fastq.gz                   # Cleaned and quality-filtered reads
│   ├── fastp.{json,html}            # Trimming statistics and adapter metrics
│   ├── QC/ 
│   │   └── *_trimmed_fastqc.html    # Quality control reports after trimming validation
├── ALIGNMENT/
│   ├── *.bam                        # Coordinate-sorted genomic alignments
│   └── *.bam.bai                    # Coordinate indexes for IGV browser inspection
├── QUANTIFICATION/
│   ├── counts.txt                   # Raw gene-level count matrix from featureCounts
│   └── counts.txt.summary           # Assignment metrics (Assigned, Unassigned, etc.)
└── DIFF_EXPR/
    ├── rds_files/
    │   └── *.rds                    # Serialized DESeq2 R objects for downstream analysis
    ├── tables/
    │   ├── *_VS_*.csv               # Statistical tables per contrast
    │   ├── significative_genes.xlsx # Summary of significant DEGs
    │   └── log2FC_matrix.csv        # Cleaned Log2FC matrix
    └── plots/
        ├── PCA.pdf                  # Unsupervised Principal Component Analysis
        ├── heatmap_samples_vst.pdf  # VST profile heatmap of DEGs
        └── heatmap_conditions_log2FC.pdf  # Expression profile heatmap of DEGs

## 9. 🔧 TROUBLESHOOTING NOTES

### A. Automatic Dot Rectification
During data acquisition, the R's native `read.delim` automatically converts non-standard characters (e.g., `/`) with dots (`.`) in column names. The script `scripts/Differential_Expression_Analysis.R` handels this behaviour by matching the converted column names using regular expressions. 

```R
correct_pattern <- paste0("^runs\\.", config$run_name, "\\.", config$output_dir, "\\.ALIGNMENT\\.(.*)\\.bam$")
colnames(counts) <- sub(correct_pattern, "\\1", colnames(counts))
```

This guarantees that column identification perfectly matches your `metadata.tsv`sample list.

## 10. 💬 SUPPORT & CONTACT
If you encounter any issues, broken links, or bugs while running **DETECT**, please feel free to:
* Open an **[issue](http://github.com/michaelbianchi7564/detect/issues)** directly on this repository
* Contact me at: `michael.bianchi@unitus.it`

*Note: Benchmark download scripts are currently undergoing maintenance for updated dataset URLs and will be updated in upcoming releases*

## 11. 📜 LICENSE & CITATION
### Legal Deposit & Licensing
This pipeline code and its automated infrastructure are deposited and archived under the **GNU General Public License v3.0 (GPLv3)**. You are free to run, study, share, and modify this pipeline for research and open-source use, provided that any derivative work preserves the original copyright notice and is distributed under the same GPLv3 licence. See the `LICENSE` file in the root repository for full details.

For proprietary, closed-source, or commercial enterprise licensing options, please contact the creator & maintainer directly at: `michael.bianchi@unitus.it`

### Citation
If you use this pipeline architecture, automated wrapper script, or downstream DESeq2 components in your publication or research project, please cite this repository:

* **Bianchi, M.** (2026). *DETECT: Automated and Isolated End-to-End RNA-seq Pipeline using Snakemake and DESeq2*. GitHub Repository: [https://github.com/michaelbianchi7564/detect](https://github.com/michaelbianchi7564/detect)
