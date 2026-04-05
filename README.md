# HybridSeqFinder — NCBI BLAST → Hybrid SRA Pair Discovery Pipeline

> Submit any nucleotide sequences to NCBI, find which genomes contain them, and automatically discover **hybrid sequencing datasets** (long-read + short-read pairs) on those genomes — ready for high-quality assembly.

---

## What does it do?

Given a set of nucleotide sequences in FASTA format, this pipeline:

1. **BLASTs** each sequence against NCBI `nt` (or any NCBI nucleotide database)
2. **Retrieves** every NCBI nucleotide accession that contains a significant match
3. **Resolves** each accession to its BioSample via Entrez
4. **Scans** SRA for all sequencing runs on that BioSample
5. **Identifies hybrid pairs** — BioSamples with both long-read (Nanopore/PacBio) AND short-read (Illumina) data deposited
6. **Merges** everything into a single results table and generates publication-quality figures

The output is a curated list of SRA accession pairs (one long-read + one short-read run per BioSample) that can be directly downloaded for hybrid genome assembly.

---

## Use cases

- Find all publicly available hybrid sequencing datasets that contain your gene, operon, or genomic element of interest
- Discover which organisms carry your sequence AND have hybrid assembly-ready data in SRA
- Screen NCBI for hybrid assembly candidates before planning a wet-lab experiment
- Epidemiological surveillance: track a mobile genetic element across organisms and countries

### Example: IS26 Composite Transposons

The default configuration targets **Novel IS26 composite transposons** — multi-kilobase genomic islands carrying IS26 insertion sequences flanking antimicrobial resistance (AMR) gene clusters. These elements are difficult to fully resolve with short reads alone; hybrid assemblies are essential for accurate structural characterisation.

`scripts/01_extract_fasta.py` shows how to extract sequences of interest from GenBank assemblies. **Replace this step with your own FASTA file** to use the pipeline for any sequence set.

---

## Pipeline

```
Your sequences (FASTA)
        │
        ▼
  01_extract_fasta.py    ← optional: extract from GenBank; or bring your own FASTA
        │
        ▼
  02_run_blast.py        → results/blast_results/*.xml.gz
        │                   (remote blastn vs NCBI nt, one XML per query)
        ▼
  03_parse_blast.py      → results/blast_hits_filtered.tsv
        │                   (apply identity + coverage thresholds)
        ▼
  04_sra_metadata.py     → results/sra_metadata.tsv
        │                   results/paired_runs.tsv
        │                   (Entrez chain: accession → BioSample → SRA runs → pairs)
        ▼
  05_merge_results.py    → results/final_results.tsv
        │                   results/final_results.xlsx
        ▼
  06_visualize.py        → results/figures/
                            (6 publication-quality figures)
```

---

## Repository Structure

```
.
├── config.py                  # All parameters: paths, BLAST settings, thresholds
├── run_pipeline.sh            # Run all steps sequentially (tmux recommended)
├── scripts/
│   ├── 01_extract_fasta.py    # Example: extract sequences from GenBank assemblies
│   ├── 02_run_blast.py        # Remote NCBI BLAST (blastn vs nt)
│   ├── 03_parse_blast.py      # Parse XML, filter by identity + query coverage
│   ├── 04_sra_metadata.py     # Entrez: accession → BioSample → SRA → hybrid pairs
│   ├── 05_merge_results.py    # Join all tables → TSV + Excel
│   └── 06_visualize.py        # 6 publication-quality figures
└── results/                   # Auto-created during run
    ├── blast_hits_filtered.tsv
    ├── sra_metadata.tsv
    ├── paired_runs.tsv
    ├── final_results.tsv
    ├── final_results.xlsx
    ├── blast_results/         # Per-query BLAST XML (gzipped, auto-skipped on re-run)
    ├── figures/
    └── logs/
```

---

## Quick Start

```bash
# 1. Install dependencies
pip install biopython pandas numpy matplotlib seaborn plotly openpyxl tqdm folium scipy

# 2. Set your NCBI email in config.py (required by NCBI)
#    Optionally add an API key for 10x faster Entrez queries

# 3. Put your sequences in results/novel_is26.fasta
#    (or run 01_extract_fasta.py to generate them)

# 4. Run (use tmux — BLAST step takes hours)
bash run_pipeline.sh
```

---

## Configuration (`config.py`)

```python
ENTREZ_EMAIL   = "your@email.com"     # Required by NCBI
ENTREZ_API_KEY = ""                   # Optional — get free at ncbi.nlm.nih.gov/account

# BLAST
BLAST_DATABASE  = "nt"                # NCBI nucleotide database
BLAST_HITLIST   = 500                 # Max hits per query
BLAST_EVALUE    = "0.001"
BLAST_ENTREZ_Q  = "Bacteria[Organism]"  # Restrict search (change as needed)

# Quality filters
MIN_QUERY_COV  = 50.0    # % query coverage  (lower for long/mosaic sequences)
MIN_IDENTITY   = 80.0    # % nucleotide identity

# Hybrid pairing
LONG_READ_PLATFORMS  = {"OXFORD_NANOPORE", "PACBIO_SMRT"}
SHORT_READ_PLATFORMS = {"ILLUMINA", "ION_TORRENT"}
```

---

## Output

### `paired_runs.tsv` — the core output

Each row is one hybrid pair: a BioSample that has both a long-read and a short-read SRA run containing your sequence.

| Column | Description |
|---|---|
| `hit_accession` | NCBI nucleotide accession matching your query |
| `biosample` | BioSample accession (SAMN…) |
| `long_read_run` | Long-read SRA run (SRR/ERR/DRR) |
| `long_read_platform` | `OXFORD_NANOPORE` or `PACBIO_SMRT` |
| `long_read_bases_Mbp` | Data volume (Mbp) |
| `short_read_run` | Short-read SRA run |
| `short_read_platform` | `ILLUMINA` or `ION_TORRENT` |
| `organism` | Source organism |
| `country` | Country of isolation |
| `isolation_source` | Sample type (e.g. blood, urine, stool) |
| `host` | Host organism |
| `collection_date` | Date of sample collection |

### `final_results.xlsx`

Two-sheet Excel workbook:
- **Sheet 1** — All BLAST hits with full metadata
- **Sheet 2** — Hybrid pairs enriched with query metadata

### Figures (`results/figures/`)

| File | Description |
|---|---|
| `fig1_geo_map.html` + `.png` | Global distribution of BLAST hits (interactive + static) |
| `fig2_is26_organism_heatmap.png` | Query type × organism clustered heatmap |
| `fig3_amr_heatmap.png` | Gene × organism co-occurrence (hierarchical clustering) |
| `fig4_sequencing_landscape.png` | Platform breakdown + data volume (dual panel) |
| `fig5_blast_quality.png` | Identity vs coverage with marginal histograms |
| `fig6_hybrid_overview.png` | Hybrid pair count + sequencing volume per organism (lollipop) |

---

## Design Notes

**BioSample-level pairing:** Hybrid pairs are identified at the BioSample level — not BioProject. A BioSample represents a single biological specimen, so long and short reads from the same BioSample are directly comparable for hybrid assembly. BioProject scanning is intentionally excluded.

**Query coverage threshold:** Long or mosaic sequences rarely achieve 80%+ coverage in a single BLAST hit. A default of 50% is used to capture partial but biologically meaningful matches. Adjust `MIN_QUERY_COV` in `config.py` as needed.

**Resume support:** Steps 02 (BLAST) and 04 (SRA metadata) save results incrementally. Re-running after an interruption will continue from where it stopped — no redundant NCBI queries.

---

## Requirements

| Package | Version | Purpose |
|---|---|---|
| `biopython` | ≥1.79 | SeqIO, NCBIWWW, NCBIXML, Entrez |
| `pandas` | ≥1.3 | Data manipulation |
| `numpy` | ≥1.21 | Numerical operations |
| `matplotlib` | ≥3.4 | Static figures |
| `seaborn` | ≥0.11 | Statistical visualisations + clustermap |
| `scipy` | ≥1.7 | Hierarchical clustering |
| `plotly` | ≥5.0 | Interactive HTML figures |
| `openpyxl` | ≥3.0 | Excel export |
| `tqdm` | ≥4.0 | Progress bars |
| `folium` | ≥0.12 | Interactive geographic map |

Python 3.8+ required.

---

## Citation

If you use this pipeline in your research, please cite:

> Altınbağ RC. — HybridSeqFinder: NCBI BLAST to hybrid SRA pair discovery pipeline. 2025.
> https://github.com/recepcanaltinbag/is26-hybrid-sra-pipeline

---

## License

MIT License — see [LICENSE](LICENSE) for details.
