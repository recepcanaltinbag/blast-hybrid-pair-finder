"""
STEP 3 — Parse BLAST XML results, apply quality filters, extract accessions.

query_id recovery strategy:
  NCBI BLAST replaces '|' in FASTA headers with 'No definition line'.
  We recover the real query_id by matching sequence length from the FASTA file
  to the XML filename stem (which was set from the FASTA record id).

Filters: Query coverage >= MIN_QUERY_COV, Identity >= MIN_IDENTITY

Output: results/blast_hits_filtered.tsv
"""

import os
import sys
import gzip
import glob
from typing import List, Dict

import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import BLAST_XML_DIR, BLAST_HITS, INPUT_FASTA, MIN_QUERY_COV, MIN_IDENTITY

G   = "\033[92m"; Y = "\033[93m"; R = "\033[91m"
C   = "\033[96m"; W = "\033[97m"; DIM = "\033[2m"; RST = "\033[0m"; BOLD = "\033[1m"

def hr(char="─", width=62):
    return DIM + char * width + RST


def parse_organism(hit_def: str) -> str:
    parts = hit_def.split("[")
    if len(parts) >= 2:
        return parts[-1].rstrip("]").strip()
    return " ".join(hit_def.split()[:2])


def build_xml_to_query_map(fasta_path: str, xml_dir: str) -> Dict[str, str]:
    """
    Map each XML filename stem → real query_id from FASTA.

    XML files are named after the FASTA record id with '|' replaced by '_'.
    We match by comparing the sanitised FASTA id to the XML stem.
    """
    mapping = {}

    # Load FASTA ids
    records = list(SeqIO.parse(fasta_path, "fasta"))
    # Build lookup: sanitised_id → original_id
    sanitised_to_real = {}
    for rec in records:
        safe = rec.id.replace("|", "_").replace("/", "_")
        sanitised_to_real[safe] = rec.id

    # Match XML filenames
    for xml_path in glob.glob(os.path.join(xml_dir, "*.xml.gz")):
        stem = os.path.basename(xml_path).replace(".xml.gz", "")
        if stem in sanitised_to_real:
            mapping[xml_path] = sanitised_to_real[stem]
        else:
            # Fallback: try prefix match on SRAID|CycleID
            matched = None
            for safe, real in sanitised_to_real.items():
                if stem.startswith(safe[:30]):
                    matched = real
                    break
            mapping[xml_path] = matched or stem   # keep stem if no match

    return mapping


def parse_xml_file(raw: str, query_id: str,
                   min_qcov: float, min_ident: float) -> List[Dict]:
    import io
    rows = []
    try:
        blast_records = list(NCBIXML.parse(io.StringIO(raw)))
    except Exception as e:
        tqdm.write(f"  {R}[!] Could not parse XML for {query_id}: {e}{RST}")
        return rows

    for blast_record in blast_records:
        query_len = blast_record.query_length

        for alignment in blast_record.alignments:
            accession = alignment.accession
            hit_title = alignment.title[:120]
            organism  = parse_organism(alignment.hit_def)

            for hsp in alignment.hsps:
                if query_len == 0:
                    continue
                identity_pct  = (hsp.identities / hsp.align_length) * 100
                query_cov_pct = (hsp.align_length / query_len) * 100

                if identity_pct < min_ident or query_cov_pct < min_qcov:
                    continue

                rows.append({
                    "query_id"      : query_id,
                    "query_len"     : query_len,
                    "hit_accession" : accession,
                    "hit_title"     : hit_title,
                    "organism"      : organism,
                    "evalue"        : hsp.expect,
                    "bitscore"      : hsp.bits,
                    "identity_pct"  : round(identity_pct, 2),
                    "query_cov_pct" : round(query_cov_pct, 2),
                    "align_len"     : hsp.align_length,
                    "hit_start"     : hsp.sbjct_start,
                    "hit_end"       : hsp.sbjct_end,
                })
    return rows


def main():
    xml_files = sorted(glob.glob(os.path.join(BLAST_XML_DIR, "*.xml.gz")))

    print()
    print(hr("═"))
    print(f"  {BOLD}{C}STEP 03 — BLAST RESULT PARSER{RST}")
    print(hr("═"))
    print(f"  {W}XML files found  :{RST} {len(xml_files)}")
    print(f"  {W}Identity filter  :{RST} >= {MIN_IDENTITY}%")
    print(f"  {W}Coverage filter  :{RST} >= {MIN_QUERY_COV}%")
    print(hr())

    if not xml_files:
        print(f"  {R}ERROR: No XML files in {BLAST_XML_DIR}{RST}")
        sys.exit(1)

    if not os.path.exists(INPUT_FASTA):
        print(f"  {R}ERROR: FASTA not found: {INPUT_FASTA}{RST}")
        sys.exit(1)

    print(f"  {DIM}Building XML → query_id mapping from FASTA...{RST}")
    xml_to_query = build_xml_to_query_map(INPUT_FASTA, BLAST_XML_DIR)
    matched = sum(1 for v in xml_to_query.values() if "|" in v)
    print(f"  {G}Matched:{RST} {matched}/{len(xml_files)} XML files to real query IDs\n")

    all_rows  = []

    bar_fmt = f"  {C}{{l_bar}}{{bar}}{RST}{DIM}| {{n_fmt}}/{{total_fmt}} [{{elapsed}}<{{remaining}}]{RST}"

    with tqdm(xml_files, bar_format=bar_fmt, ncols=68, unit="xml") as pbar:
        for xml_path in pbar:
            query_id = xml_to_query.get(xml_path, os.path.basename(xml_path))
            label    = query_id.split("|")[0][:35] if "|" in query_id else query_id[:35]
            pbar.set_description(f"  {C}Parsing{RST} {label:<35}")

            with gzip.open(xml_path, "rt", encoding="utf-8") as fh:
                raw = fh.read()
            total_hits = raw.count("<Hit>")

            rows = parse_xml_file(raw, query_id, MIN_QUERY_COV, MIN_IDENTITY)
            all_rows.extend(rows)

            passing = len(rows)
            pct     = (passing / total_hits * 100) if total_hits else 0
            status  = G if passing > 0 else Y

            tqdm.write(
                f"  {status}[{'✓' if passing else '○'}]{RST} "
                f"{label:<38} "
                f"{DIM}total:{RST}{total_hits:>4}  "
                f"{DIM}pass:{RST}{status}{passing:>4}{RST}  "
                f"{DIM}({pct:>5.1f}%){RST}"
            )

    if not all_rows:
        print(f"\n  {Y}WARNING: No hits passed the filters.{RST}")
        sys.exit(0)

    df = (pd.DataFrame(all_rows)
            .sort_values("bitscore", ascending=False)
            .drop_duplicates(subset=["query_id", "hit_accession"])
            .reset_index(drop=True))

    os.makedirs(os.path.dirname(BLAST_HITS), exist_ok=True)
    df.to_csv(BLAST_HITS, sep="\t", index=False)

    ident_mean = df["identity_pct"].mean()
    ident_med  = df["identity_pct"].median()
    cov_mean   = df["query_cov_pct"].mean()
    top_orgs   = df["organism"].value_counts().head(5)

    print()
    print(hr("═"))
    print(f"  {BOLD}{G}PARSE COMPLETE{RST}")
    print(hr("═"))
    print(f"  {W}Total passing hits      :{RST} {G}{len(df)}{RST}")
    print(f"  {W}Unique queries with hits:{RST} {df['query_id'].nunique()} / {len(xml_files)}")
    print(f"  {W}Unique accessions       :{RST} {df['hit_accession'].nunique()}")
    print(f"  {W}Unique organisms        :{RST} {df['organism'].nunique()}")
    print(hr())
    print(f"  {W}Identity   — mean:{RST} {ident_mean:.1f}%  {W}median:{RST} {ident_med:.1f}%")
    print(f"  {W}Query cov  — mean:{RST} {cov_mean:.1f}%")
    print(hr())
    print(f"  {W}Top 5 organisms:{RST}")
    for org, cnt in top_orgs.items():
        bar = "█" * min(int(cnt / top_orgs.max() * 20), 20)
        print(f"    {G}{bar:<20}{RST} {cnt:>4}  {DIM}{org}{RST}")
    print(hr("═"))
    print(f"\n  {G}→ Saved:{RST} {BLAST_HITS}")
    print(f"  {G}→ Next :{RST} python3 scripts/04_sra_metadata.py\n")


if __name__ == "__main__":
    main()
