"""
STEP 1 — Extract query sequences from source GenBank files.

This is an example extractor for a specific input format (composite transposon CSV + GBK).
Replace this step with your own FASTA preparation if your input format differs.

Output: results/novel_is26.fasta   (path configurable via INPUT_FASTA in config.py)
  Header format: >SRAID|CycleID|element_type|resistance_genes
"""

import os
import sys
import glob
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import GBK_DIR, SUMMARY_CSV, INPUT_FASTA, RESULTS_DIR


def load_novel_is26(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    mask = (df["CompTN"] == "Novel") & (df["ISproducts"].str.contains("IS26", na=False))
    novel = df[mask].copy()
    print(f"[01] Query entries in CSV: {len(novel)}")
    return novel


def build_gbk_key(sra_id: str, kmer: str, cycle_id: str) -> str:
    """Reconstruct the GBK filename stem (without .gbk)."""
    return f"{sra_id}__{sra_id}_{kmer}__{cycle_id}"


def extract_sequences(novel_df: pd.DataFrame, gbk_dir: str, fasta_out: str) -> int:
    os.makedirs(os.path.dirname(fasta_out), exist_ok=True)

    written = 0
    missing = []

    with open(fasta_out, "w") as fh:
        for _, row in novel_df.iterrows():
            sra_id   = str(row["SRAID"]).strip()
            kmer     = str(row["kmer"]).strip()
            cycle_id = str(row["CycleID"]).strip()
            is_type  = str(row["ISproducts"]).strip().replace(" ", "_")
            ant      = str(row.get("Antproducts_names", "")).strip().replace(";", "|")

            stem    = build_gbk_key(sra_id, kmer, cycle_id)
            gbk_path = os.path.join(gbk_dir, stem + ".gbk")

            if not os.path.exists(gbk_path):
                missing.append(stem)
                continue

            record = SeqIO.read(gbk_path, "genbank")
            seq    = str(record.seq)
            if not seq or "N" * len(seq) == seq:
                continue

            # Sanitise header components
            ant_safe = ant if ant and ant != "nan" else "no_AMR"
            header   = f">{sra_id}|{cycle_id}|{is_type}|{ant_safe}"

            fh.write(header + "\n")
            # 60-char wrapped FASTA
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")
            written += 1

    print(f"[01] Written to FASTA : {written}")
    if missing:
        print(f"[01] WARNING — {len(missing)} GBK file(s) not found:")
        for m in missing:
            print(f"       {m}")

    return written


def main():
    novel_df = load_novel_is26(SUMMARY_CSV)
    n = extract_sequences(novel_df, GBK_DIR, INPUT_FASTA)
    if n == 0:
        print("[01] ERROR: No sequences written. Check paths.")
        sys.exit(1)
    print(f"[01] Done → {INPUT_FASTA}")


if __name__ == "__main__":
    main()
