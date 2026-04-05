"""
HybridSeqFinder — NCBI BLAST & SRA Hybrid Pair Discovery Pipeline
Central configuration — edit before running.
"""

import os

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR      = os.path.dirname(os.path.abspath(__file__))
GBK_DIR       = os.path.join(BASE_DIR, "ıs26")          # Turkish 'ı' — as-is
SUMMARY_CSV   = os.path.join(BASE_DIR, "merged_summary.csv")
RESULTS_DIR   = os.path.join(BASE_DIR, "results")
BLAST_XML_DIR = os.path.join(RESULTS_DIR, "blast_results")
FIGURES_DIR   = os.path.join(RESULTS_DIR, "figures")

INPUT_FASTA   = os.path.join(RESULTS_DIR, "novel_is26.fasta")
BLAST_HITS    = os.path.join(RESULTS_DIR, "blast_hits_filtered.tsv")
SRA_META      = os.path.join(RESULTS_DIR, "sra_metadata.tsv")
FINAL_TSV     = os.path.join(RESULTS_DIR, "final_results.tsv")
FINAL_EXCEL   = os.path.join(RESULTS_DIR, "final_results.xlsx")

# ── NCBI / Entrez ──────────────────────────────────────────────────────────
# REQUIRED: NCBI requires an email for Entrez queries.
# API key is optional but raises limit from 3 → 10 req/s.
# Get one free at: https://www.ncbi.nlm.nih.gov/account/
ENTREZ_EMAIL   = "recocanextra@gmail.com"   # <-- CHANGE THIS
ENTREZ_API_KEY = ""                          # leave empty if none

# ── BLAST parameters ───────────────────────────────────────────────────────
BLAST_PROGRAM   = "blastn"
BLAST_DATABASE  = "nt"
BLAST_HITLIST   = 500          # max hits per query
BLAST_EVALUE    = "0.001"
BLAST_ENTREZ_Q  = "Bacteria[Organism]"   # restrict to bacteria

# Minimum thresholds for keeping a hit
MIN_QUERY_COV   = 50.0   # %
MIN_IDENTITY    = 80.0   # %

# Seconds to wait between BLAST submissions (no API key → max 3 req/s)
BLAST_DELAY     = 20     # conservative; NCBI jobs take time anyway

# ── SRA filter ─────────────────────────────────────────────────────────────
LONG_READ_PLATFORMS = {
    "OXFORD_NANOPORE",
    "PACBIO_SMRT",
}
SHORT_READ_PLATFORMS = {
    "ILLUMINA",
    "ION_TORRENT",
}

# Minimum long+short pairs to report per query sequence
MIN_PAIRS_PER_QUERY = 5

# ── Additional output paths ────────────────────────────────────────────────
PAIRED_RUNS  = os.path.join(RESULTS_DIR, "paired_runs.tsv")
PAIRED_EXCEL = os.path.join(RESULTS_DIR, "paired_runs.xlsx")
