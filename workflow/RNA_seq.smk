
import os
import sys

from snakemake.io import glob_wildcards                                 # Import utility to extract wildcard values from file patterns
from snakemake.shell import shell                                       # Import utility to use shell during hisat-building

CWD = os.getcwd()
if CWD not in sys.path:
    sys.path.append(CWD)

from scripts.utils import (                                                     # Import helper functions to get parameters from configfile
    validate_config,
    dict_to_cli,
    add_extra,
    build_hisat2_cli,
    featurecounts_cli
)

configfile: "config.yaml"

if "run_name" not in config:
    raise ValueError(
        "\n\n[ERROR] 'run_name' parameter is not defined.\n"
        "Do not run the snakemake file directly, use the wrapper script instead:\n"
        "./run_pipeline.sh {RUN} n_cores\n")

validate_config(config['featurecounts'])

RUN_DIR = f"runs/{config['run_name']}"
INPUT_DIR = f"{RUN_DIR}/{config['input_dir']}"
OUTPUT_DIR = f"{RUN_DIR}/{config['output_dir']}"
DIFF_EXPR_DIR = f"{OUTPUT_DIR}/DIFF_EXPR"
REFERENCE_DIR = f"reference_genomes/{config['species']}"

SAMPLES_R1 = set(glob_wildcards(f"{INPUT_DIR}/{{sample}}_R1.fq.gz").sample)    # Detect sample IDs from R1 reads FASTQ files in DATA/
SAMPLES_R2 = set(glob_wildcards(f"{INPUT_DIR}/{{sample}}_R2.fq.gz").sample)    # Detect sample IDs from R2 reads FASTQ files in DATA/
SAMPLES = sorted(SAMPLES_R1 & SAMPLES_R2)                               # Keep only samples that have both R1 and R2 files and sort them

print(f"Detected {len(SAMPLES)} paired-end samples: {', '.join(SAMPLES)}")      # Print the number of samples detected

if not SAMPLES:
    raise ValueError(f"No paired-end samples found in {INPUT_DIR}/")            # Stop execution if no samples are found

rule all:
    input:
        expand(f"{OUTPUT_DIR}/QC/{{sample}}_R1_fastqc.html", sample = SAMPLES),
        expand(f"{OUTPUT_DIR}/QC/{{sample}}_R2_fastqc.html", sample = SAMPLES),
        expand(f"{RUN_DIR}/HISAT_INDEX/{{species}}.{{i}}.ht2", species = config['species'], i = range(1,9)),
        expand(f"{OUTPUT_DIR}/ALIGNMENT/QC_PICARD/{{sample}}.metrics", sample = SAMPLES),
        f"{OUTPUT_DIR}/QUANTIFICATION/counts.txt",
        f"{DIFF_EXPR_DIR}/tables/significative_genes.xlsx",
        f"{DIFF_EXPR_DIR}/tables/log2FC_matrix.csv",
        f"{DIFF_EXPR_DIR}/tables/vsdata_matrix.csv",
        f"{DIFF_EXPR_DIR}/tables/cluster_assignment.csv"

rule quality_check:
    input:
        R1 = f"{INPUT_DIR}/{{sample}}_R1.fq.gz",
        R2 = f"{INPUT_DIR}/{{sample}}_R2.fq.gz",
    output:
        R1 = f"{OUTPUT_DIR}/QC/{{sample}}_R1_fastqc.html",
        R2 = f"{OUTPUT_DIR}/QC/{{sample}}_R2_fastqc.html",
    threads:
        config["threads"]["fastqc"]
    params:
        fastqc_params = add_extra("", config["fastqc"])
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/QC
        fastqc -t {threads} -o {OUTPUT_DIR}/QC/ {params.fastqc_params} {input.R1} {input.R2}
        """

rule trimming:
    input:
        R1 = f"{INPUT_DIR}/{{sample}}_R1.fq.gz",
        R2 = f"{INPUT_DIR}/{{sample}}_R2.fq.gz",
    output:
        R1P = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R1_p.fastq",
        R2P = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R2_p.fastq",
        R1U = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R1_u.fastq",
        R2U = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R2_u.fastq",
        REPORT = f"{OUTPUT_DIR}/TRIMMING/{{sample}}.html",
        JSON = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_report.json",
    threads:
        config["threads"]["fastp"]
    params:
        fastp_params = add_extra(dict_to_cli(config["fastp"]), config["fastp"])
    conda:
        "envs/fastp.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/TRIMMING/
        fastp -i {input.R1} -I {input.R2} -o {output.R1P} -O {output.R2P} --unpaired1 {output.R1U} --unpaired2 {output.R2U} -h {output.REPORT} -j {output.JSON} -w {threads} {params.fastp_params}
        """

rule qc_after_trimming:
    input:
        R1P = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R1_p.fastq",
        R2P = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R2_p.fastq",
        R1U = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R1_u.fastq",
        R2U = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R2_u.fastq",
    output:
        R1P = f"{OUTPUT_DIR}/TRIMMING/QC/{{sample}}_R1_p_fastqc.html",
        R2P = f"{OUTPUT_DIR}/TRIMMING/QC/{{sample}}_R2_p_fastqc.html",
        R1U = f"{OUTPUT_DIR}/TRIMMING/QC/{{sample}}_R1_u_fastqc.html",
        R2U = f"{OUTPUT_DIR}/TRIMMING/QC/{{sample}}_R2_u_fastqc.html",
    threads:
        config["threads"]["fastqc"]
    params:
        fastqc_params = add_extra("", config["fastqc"])
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/TRIMMING/QC/
        fastqc -t {threads} -o {OUTPUT_DIR}/TRIMMING/QC/ {params.fastqc_params} {input}
        """

rule index_building:
    input:
        genome=f"{REFERENCE_DIR}/{config['reference']['genome']}",
        annotation=f"{REFERENCE_DIR}/{config['reference']['annotation']}",
    output:
        expand(f"{RUN_DIR}/HISAT_INDEX/{{species}}.{{i}}.ht2", species=config['species'], i=range(1,9))
    params:
        index_prefix=f"{RUN_DIR}/HISAT_INDEX/{config['species']}",
        building_params = add_extra(dict_to_cli(config["index_building"]), config["index_building"]),
    conda:
        "envs/hisat2.yaml"
    shell:
        """
        mkdir -p {RUN_DIR}/HISAT_INDEX/
        hisat2-build {input.genome} {params.index_prefix} {params.building_params}
        echo "index built"
        """

rule align:
    input:
        R1P = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R1_p.fastq",
        R2P = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R2_p.fastq",
        R1U = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R1_u.fastq",
        R2U = f"{OUTPUT_DIR}/TRIMMING/{{sample}}_R2_u.fastq",
        INDEX = expand(f"{RUN_DIR}/HISAT_INDEX/{{species}}.{{i}}.ht2", species=config['species'], i=range(1,9))
    output:
        SAM_FILE = temp(f"{OUTPUT_DIR}/ALIGNMENT/{{sample}}.sam"),
        SUMMARY = f"{OUTPUT_DIR}/ALIGNMENT/REPORTS/{{sample}}.txt"
    threads:
        config["threads"]["hisat2"]
    params:
        index_prefix = f"{RUN_DIR}/HISAT_INDEX/{config['species']}",
        alignment_params = build_hisat2_cli(config["hisat2"])
    conda:
        "envs/hisat2.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/ALIGNMENT/REPORTS/
        hisat2 -q --very-sensitive -p {threads} -x {params.index_prefix} -1 {input.R1P} -2 {input.R2P} -U {input.R1U},{input.R2U} -S {output.SAM_FILE} --summary-file {output.SUMMARY} {params.alignment_params}
        """

rule binary:
    input:
        SAM_FILE = f"{OUTPUT_DIR}/ALIGNMENT/{{sample}}.sam"
    output:
        SORTED = f"{OUTPUT_DIR}/ALIGNMENT/{{sample}}.bam",
        INDEX = f"{OUTPUT_DIR}/ALIGNMENT/{{sample}}.bam.bai"
    threads:
        config["threads"]["samtools_sort"]
    params:
        sort_params = add_extra("", config["samtools_sort"]),
        index_params = add_extra("", config["samtools_index"])
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output.SORTED} {input.SAM_FILE} {params.sort_params}
        samtools index -@ {threads} {output.SORTED} {params.index_params}
        """

rule alignment_quality_check:
    input:
        SORTED = f"{OUTPUT_DIR}/ALIGNMENT/{{sample}}.bam",
        REFERENCE = f"{REFERENCE_DIR}/{config['reference']['genome']}"
    output:
        REPORT = f"{OUTPUT_DIR}/ALIGNMENT/QC_PICARD/{{sample}}.metrics"
    conda:
        "envs/picard.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/ALIGNMENT/QC_PICARD/
        picard CollectAlignmentSummaryMetrics I={input.SORTED} O={output.REPORT} R={input.REFERENCE}
        """

rule quantification:
    input:
        SORTED = expand(f"{OUTPUT_DIR}/ALIGNMENT/{{sample}}.bam", sample = SAMPLES),
        ANNOTATION = f"{REFERENCE_DIR}/{config['reference']['annotation_tab']}"
    output:
        COUNT_MATRIX = f"{OUTPUT_DIR}/QUANTIFICATION/counts.txt"
    threads:
        config["threads"]["featureCounts"]
    params:
        featurecounts_params = add_extra(featurecounts_cli(config["featurecounts"]), config["featurecounts"])
    conda:
        "envs/subread.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/QUANTIFICATION/
        featureCounts -p -T {threads} -a {input.ANNOTATION} -o {output.COUNT_MATRIX} {input.SORTED} {params.featurecounts_params} 2> {OUTPUT_DIR}/QUANTIFICATION/featureCounts.log
        """

rule differential_exp:
    input:
        COUNTS = f"{OUTPUT_DIR}/QUANTIFICATION/counts.txt",
        METADATA = f"{RUN_DIR}/metadata.tsv",
        CONFIG = f"{RUN_DIR}/config.yaml",
    output:
        PCA = f"{DIFF_EXPR_DIR}/plots/PCA.pdf",
        SIG_GENES = f"{DIFF_EXPR_DIR}/tables/significative_genes.xlsx",
        LOG2FC_MATRIX = f"{DIFF_EXPR_DIR}/tables/log2FC_matrix.csv",
        VST_MATRIX = f"{DIFF_EXPR_DIR}/tables/vsdata_matrix.csv",
        CLUSTERS = f"{DIFF_EXPR_DIR}/tables/cluster_assignment.csv",
    threads:
        config["threads"]["deseq2"]
    conda:
        "envs/R.yaml"
    shell:
        """
        mkdir -p {DIFF_EXPR_DIR}/tables {DIFF_EXPR_DIR}/genes {DIFF_EXPR_DIR}/plots/volcano_plots {DIFF_EXPR_DIR}/debug {DIFF_EXPR_DIR}/rds_files
        Rscript scripts/Differential_Expression_Analysis.R {input.COUNTS} {input.METADATA} {input.CONFIG}
        """

