#!/bin/bash

set -e
# ============================================================================================
# INTRO
# ============================================================================================
show_banner() {
	echo -e "\033[0;36m"

	cat << 'EOF'
========================================================================
  ____  _____ _____ _____ ____ _____ 
 |  _ \| ____|_   _| ____/ ___|_   _|
 | | | |  _|   | | |  _|| |     | |  
 | |_| | |___  | | | |__| |___  | |  
 |____/|_____| |_| |_____\____| |_|  

 DETECT: Differential Expression & Transcriptomic Analysis Pipeline
 Version:  v1.0.0
 Author:   Michael Bianchi
 GitHub:   https://github.com/michaelbianchi7564/detect
========================================================================
EOF
	
	echo -e '\033[0m'
}

show_banner

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_ROOT"

# 1. Check if at least run name is provided
if [ -z "$1" ]; then
    echo "Usage: ./run_pipeline.sh <run_name> [cores]"
    echo "Example: ./run_pipeline.sh E_coli_test 8"
    exit 1
fi

RUN=$1
CORES=${2:-4}  # Use the second arg, or use 4 cores by default

# 2. Check if a specific config file exists
LOCAL_CONFIG="runs/${RUN}/config.yaml"
if [ ! -f "$LOCAL_CONFIG" ]; then
    echo "ERROR: Configuration file not found at: $LOCAL_CONFIG"
    echo "-------------------------------------------------------------------------------"
    echo "To fix this, please create a config.yaml inside your run directory."
    echo "To use the global template as a starting point please copy that file running:"
    echo "cp config/config.yaml runs/${RUN}/config.yaml"
    exit 1
fi

# Check the match between run_name in the local config and the run folder
CONFIG_RUN_NAME=$(grep -E '^run_name:' "$LOCAL_CONFIG" | sed -E 's/^run_name:[[:space:]]*"?([^"# ]*)"?.*/\1/')
if [ -n "$CONFIG_RUN_NAME" ] && [ "$CONFIG_RUN_NAME" != "$RUN" ]; then
	echo "ERROR: run_name in $LOCAL_CONFIG ('$CONFIG_RUN_NAME') does not match the run folder ('$RUN')."
	exit 1
fi

# 3. Early create a results folder (if doesn't exist) and a backup of the config file
RESULTS_DIR="runs/${RUN}/results"
mkdir -p "$RESULTS_DIR"
cp "$LOCAL_CONFIG" "$RESULTS_DIR/run_config_backup.yaml"
echo "INFO: Config file backed up to $RESULTS_DIR/run_config_backup.yaml"

# 4. Run snakemake combining global and local config files
SNAKE_ARGS=(
	-s workflow/RNA_seq.smk
	--configfile config/config.yaml
	--configfile "$LOCAL_CONFIG"
	--config run_name="${RUN}"
)

LOG_FILE="$(mktemp)"
trap 'rm -f "$LOG_FILE"' EXIT

run_snakemake() {
	snakemake "${SNAKE_ARGS[@]}" \
		--cores "$CORES" \
		--use-conda \
		--rerun-incomplete 2>&1 | tee "$LOG_FILE"
	return "${PIPESTATUS[0]}"
}

echo "=== Checking directory locks for run: ${RUN} ==="

if ! run_snakemake; then
	if grep -qE "Directory cannot be locked|LockException" "$LOG_FILE"; then
		echo "WARNING: Directory locked (presumable previous interruption). Running --unlock..."
		snakemake "${SNAKE_ARGS[@]}" --unlock
		echo "-- Unlock Complete, rerunning pipeline ---"
		run_snakemake
	else
		echo "Error: Snakemake failed for different reason then lock"
		exit 1
	fi
fi
