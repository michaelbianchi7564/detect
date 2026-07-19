#!/bin/bash

# Esci immediatamente se un comando fallisce
set -e

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_ROOT"

# 1. Controlla che sia stato passato almeno il nome della run
if [ -z "$1" ]; then
    echo "Usage: ./run_pipeline.sh <run_name> [cores]"
    echo "Example: ./run_pipeline.sh E_coli_test 8"
    exit 1
fi

RUN=$1
CORES=${2:-4}  # Usa il secondo argomento, altrimenti usa 4 core di default

# 2. Verifica che esista il file di configurazione specifico per questa esecuzione
LOCAL_CONFIG="runs/${RUN}/config.yaml"
if [ ! -f "$LOCAL_CONFIG" ]; then
    echo "ERROR: Configuration file not found at: $LOCAL_CONFIG"
    echo "-------------------------------------------------------------------------------"
    echo "To fix this, please create a config.yaml inside your run directory."
    echo "To use the global template as a starting point please copy that file running:"
    echo "cp config/config.yaml runs/${RUN}/config.yaml"
    exit 1
fi

# 3. Crea in anticipo la cartella dei risultati (se non esiste) e fai il backup del config
# In questo modo i risultati conterranno per sempre i parametri esatti di questo specifico run
RESULTS_DIR="runs/${RUN}/results"
mkdir -p "$RESULTS_DIR"
cp "$LOCAL_CONFIG" "$RESULTS_DIR/run_config_backup.yaml"
echo "INFO: Config file backed up to $RESULTS_DIR/run_config_backup.yaml"

# 4. Lancia Snakemake combinando il config globale (default) e quello locale (specifico)
echo "=== Starting Snakemake for run: ${RUN} using ${CORES} cores ==="
snakemake \
    -s workflow/RNA_seq.smk \
    --configfile config/config.yaml \
    --configfile "$LOCAL_CONFIG" \
    --config run_name="${RUN}" \
    --cores "$CORES" \
    --use-conda
