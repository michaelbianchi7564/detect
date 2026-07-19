#!/bin/bash

# Esci immediatamente se un comando fallisce
set -e

# Poiché lo script si trova dentro tests/, andiamo alla radice del progetto
PROJECT_ROOT="$(cd "$(dirname "$0")" && cd .. && pwd)"
cd "$PROJECT_ROOT"

# Destinazione isolata per il benchmark biologico completo
DATA_DIR="runs/E_coli_benchmark/DATA"
REF_DIR="runs/E_coli_benchmark/REF"

mkdir -p "$DATA_DIR"
mkdir -p "$REF_DIR"

echo "[1/2] Downloading Reference Genome (U00096.3 / MG1655)"
cd "$REF_DIR"

# Download del genoma FASTA (.fna) e dell'annotazione (.gff) da NCBI RefSeq
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz

echo "Decompressing reference files..."
gunzip -f GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip -f GCF_000005845.2_ASM584v2_genomic.gff.gz

mv GCF_000005845.2_ASM584v2_genomic.fna e_coli_k12_mg1655.fna
mv GCF_000005845.2_ASM584v2_genomic.gff e_coli_k12_mg1655.gff

echo "Converting GFF3 to GTF using gffread..."
gffread e_coli_k12_mg1655.gff -T -o e_coli_k12_mg1655.gtf

cd "$PROJECT_ROOT"

echo "[2/2] Downloading RNA-seq FASTQ files from NCBI AWS Mirror"
SRA_IDS=("SRR35290504" "SRR35290506" "SRR35290507" "SRR35290509" "SRR35290510" "SRR35290512")

# Entriamo nella cartella DATA dedicata al benchmark
cd "$DATA_DIR"

for ID in "${SRA_IDS[@]}"; do
    echo "Fetching $ID directly from NCBI Cloud Storage..."
    wget -c "https://sra-pub-run-odp.s3.amazonaws.com/sra/${ID}/${ID}"
    
    echo "Extracting reads from local SRA file $ID..."
    if command -v fasterq-dump &> /dev/null; then
        fasterq-dump "$ID" --split-files
    elif command -v fastq-dump &> /dev/null; then
        fastq-dump --split-files "$ID"
    else
        echo "ERROR: Neither fasterq-dump nor fastq-dump was found in your PATH."
        echo "Please ensure sra-tools is installed in your active environment."
        exit 1
    fi
    
    echo "Compressing and renaming to standard format..."
    gzip -f "${ID}_1.fastq"
    gzip -f "${ID}_2.fastq"
    mv "${ID}_1.fastq.gz" "${ID}_R1.fastq.gz"
    mv "${ID}_2.fastq.gz" "${ID}_R2.fastq.gz"
    
    # Libera spazio rimuovendo il file SRA nativo
    rm -f "$ID"
done

cd "$PROJECT_ROOT"

echo "========================================================="
echo " SUCCESS: All full benchmark data and reference staged!"
echo "========================================================="
