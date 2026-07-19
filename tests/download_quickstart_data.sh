#!/bin/bash

# Esci immediatamente se un comando fallisce
set -e

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_ROOT"

echo "The project root is: $PROJECT_ROOT"

# Definizione delle cartelle di destinazione corrette
DATA_DIR="runs/E_coli_test/DATA"
REF_DIR="runs/E_coli_test/REF"

mkdir -p "$DATA_DIR"
mkdir -p "$REF_DIR"

# 1. Copia dei file di riferimento locali
echo "[1/2] Staging Reference Genome & Annotation..."
cp resources/Reference_genome_E_coli_test/genome.fna "$PROJECT_ROOT/$REF_DIR/e_coli_k12_mg1655.fna"
cp resources/Reference_genome_E_coli_test/annotation.gff "$PROJECT_ROOT/$REF_DIR/e_coli_k12_mg1655.gff"
cp resources/Reference_genome_E_coli_test/annotation.gtf "$PROJECT_ROOT/$REF_DIR/e_coli_k12_mg1655.gtf"

# 2. Generazione locale delle reads sintetiche tramite Python
echo "[2/2] Generating Synthetic RNA-seq Reads (No internet required!)..."

python3 - <<EOF
import os
import gzip
import random

# Percorsi dei file
ref_path = "$PROJECT_ROOT/$REF_DIR/e_coli_k12_mg1655.fna"
data_dir = "$PROJECT_ROOT/$DATA_DIR"

# Leggiamo il genoma di riferimento
print("Reading reference genome...")
genome_seq = []
with open(ref_path, "r") as f:
    for line in f:
        if not line.startswith(">"):
            genome_seq.append(line.strip().upper())
ref_seq = "".join(genome_seq)
genome_len = len(ref_seq)
print(f"Loaded reference genome of {genome_len} bp.")

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(comp.get(base, 'N') for base in reversed(seq))

def generate_sample(sample_id, num_reads=5000, read_len=100, fragment_len=250, seed=42):
    random.seed(seed)
    r1_path = os.path.join(data_dir, f"{sample_id}_R1.fastq.gz")
    r2_path = os.path.join(data_dir, f"{sample_id}_R2.fastq.gz")
    
    print(f"Generating {num_reads} paired-end reads for {sample_id}...")
    
    with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
        for i in range(num_reads):
            # Campionamento casuale sul genoma
            start = random.randint(0, genome_len - fragment_len)
            fragment = ref_seq[start : start + fragment_len]
            
            # Read 1 e Read 2 (reverse complement)
            r1_seq = fragment[:read_len]
            r2_seq = reverse_complement(fragment[-read_len:])
            
            # Qualità fittizia alta (tutti 'I' = Phred 40)
            qual = "I" * read_len
            
            # Scrittura formato FASTQ
            r1.write(f"@{sample_id}_{i}/1\n{r1_seq}\n+\n{qual}\n")
            r2.write(f"@{sample_id}_{i}/2\n{r2_seq}\n+\n{qual}\n")

# Generiamo i 4 campioni (2 controlli, 2 trattati) con seed diversi per simulare la varianza biologica
generate_sample("SRR35290504", num_reads=10000, seed=10) # WT 1
generate_sample("SRR35290506", num_reads=10000, seed=20) # WT 2
generate_sample("SRR35290507", num_reads=10000, seed=30) # Mutant 1
generate_sample("SRR35290509", num_reads=10000, seed=40) # Mutant 2

EOF

echo "========================================================="
echo " SUCCESS: Local E. coli test environment is ready!"
echo " Offline synthetic reads generated successfully!"
echo " Ready to run Snakemake!"
echo "========================================================="
