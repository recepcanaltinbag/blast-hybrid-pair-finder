#!/usr/bin/env bash
# IS26 Composite Transposon Pipeline — run all steps in order
# Usage: bash run_pipeline.sh
# Tip: run inside tmux for overnight BLAST step

set -euo pipefail
cd "$(dirname "$0")"

PYTHON=python3
LOG_DIR="results/logs"
mkdir -p "$LOG_DIR"

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }

run_step() {
    local step=$1
    local script=$2
    echo ""
    echo "════════════════════════════════════════════════════"
    echo "  $(timestamp)  Starting: $step"
    echo "════════════════════════════════════════════════════"
    $PYTHON "scripts/$script" 2>&1 | tee "$LOG_DIR/${script%.py}.log"
    echo "  $(timestamp)  Finished: $step"
}

run_step "STEP 1 — Extract FASTA"        01_extract_fasta.py
run_step "STEP 2 — NCBI BLAST (slow!)"   02_run_blast.py
run_step "STEP 3 — Parse BLAST results"  03_parse_blast.py
run_step "STEP 4 — SRA metadata"         04_sra_metadata.py
run_step "STEP 5 — Merge & export"       05_merge_results.py
run_step "STEP 6 — Visualisations"       06_visualize.py

echo ""
echo "Pipeline complete. Results in: results/"
