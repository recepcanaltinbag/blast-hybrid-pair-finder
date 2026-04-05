"""
STEP 2 — Submit each query sequence to NCBI remote BLAST.

Strategy:
  - One query at a time (no API key → max 3 req/s)
  - Results saved as XML per query → skips already-done ones on re-run
  - Exponential backoff on network errors
  - Live progress bar + per-query stats

Usage:
    python3 scripts/02_run_blast.py

Runtime estimate (no API key, 68 seqs): ~3-6 hours depending on NCBI load.
Run overnight or in a tmux session.
"""

import os
import sys
import time
import random
import gzip
from datetime import timedelta

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio import Entrez
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    ENTREZ_EMAIL, ENTREZ_API_KEY,
    INPUT_FASTA, BLAST_XML_DIR,
    BLAST_PROGRAM, BLAST_DATABASE,
    BLAST_HITLIST, BLAST_EVALUE,
    BLAST_ENTREZ_Q, BLAST_DELAY,
)

Entrez.email   = ENTREZ_EMAIL
if ENTREZ_API_KEY:
    Entrez.api_key = ENTREZ_API_KEY

# ── Terminal colours ───────────────────────────────────────────────────────
G  = "\033[92m"   # green
Y  = "\033[93m"   # yellow
R  = "\033[91m"   # red
C  = "\033[96m"   # cyan
B  = "\033[94m"   # blue
W  = "\033[97m"   # white
DIM= "\033[2m"
RST= "\033[0m"
BOLD="\033[1m"


def hr(char="─", width=62):
    return DIM + char * width + RST


def fmt_dur(seconds: float) -> str:
    return str(timedelta(seconds=int(seconds)))


# ── Helpers ────────────────────────────────────────────────────────────────

def xml_path_for(record_id: str) -> str:
    safe = record_id.replace("|", "_").replace("/", "_")
    return os.path.join(BLAST_XML_DIR, safe + ".xml.gz")


def already_done(record_id: str) -> bool:
    return os.path.exists(xml_path_for(record_id))


def blast_one(record, pbar, attempt: int = 1):
    """Submit a single SeqRecord to BLAST. Returns raw XML string or None."""
    try:
        result_handle = NCBIWWW.qblast(
            program      = BLAST_PROGRAM,
            database     = BLAST_DATABASE,
            sequence     = str(record.seq),
            hitlist_size = BLAST_HITLIST,
            expect       = BLAST_EVALUE,
            entrez_query = BLAST_ENTREZ_Q,
            format_type  = "XML",
        )
        xml_data = result_handle.read()
        result_handle.close()
        return xml_data
    except Exception as e:
        wait = min(60 * attempt, 300)
        pbar.write(f"  {Y}[!] Attempt {attempt} failed: {e}{RST}")
        pbar.write(f"  {Y}    Retrying in {wait}s ...{RST}")
        time.sleep(wait)
        if attempt < 5:
            return blast_one(record, pbar, attempt + 1)
        pbar.write(f"  {R}[✗] Giving up after 5 attempts.{RST}")
        return None


def save_xml(record_id: str, xml_data: str):
    path = xml_path_for(record_id)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(xml_data)


def validate_xml(xml_data: str) -> bool:
    return "<BlastOutput>" in xml_data


def count_hits(xml_data: str) -> int:
    """Quick count of <Hit> tags in XML without full parse."""
    return xml_data.count("<Hit>")


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    if not os.path.exists(INPUT_FASTA):
        print(f"{R}[02] ERROR: FASTA not found: {INPUT_FASTA}{RST}")
        print(f"     Run 01_extract_fasta.py first.")
        sys.exit(1)

    os.makedirs(BLAST_XML_DIR, exist_ok=True)

    records     = list(SeqIO.parse(INPUT_FASTA, "fasta"))
    total       = len(records)
    done_before = sum(1 for r in records if already_done(r.id))
    remaining   = total - done_before

    # ── Header banner ──────────────────────────────────────────────────────
    print()
    print(hr("═"))
    print(f"  {BOLD}{C}HYBRID SEQ FINDER — NCBI BLAST{RST}")
    print(hr("═"))
    print(f"  {W}Database   :{RST} {BLAST_DATABASE}")
    print(f"  {W}Program    :{RST} {BLAST_PROGRAM}")
    print(f"  {W}Hitlist    :{RST} {BLAST_HITLIST} hits/query")
    print(f"  {W}Filter     :{RST} {BLAST_ENTREZ_Q}")
    print(f"  {W}Delay      :{RST} ~{BLAST_DELAY}s between queries (NCBI rate limit)")
    print(hr())
    print(f"  {G}Total queries  :{RST} {total}")
    print(f"  {G}Already cached :{RST} {done_before}")
    print(f"  {Y}Remaining      :{RST} {remaining}")
    est_min = remaining * (BLAST_DELAY + 90) / 60
    print(f"  {Y}Estimated time :{RST} ~{est_min:.0f} min  ({fmt_dur(est_min * 60)})")
    print(hr("═"))
    print()

    if remaining == 0:
        print(f"  {G}✓ All queries already done. Proceed to step 03.{RST}\n")
        return

    submitted  = 0
    failed     = 0
    total_hits = 0
    times      = []
    session_start = time.time()

    bar_fmt = (
        f"  {C}{{l_bar}}{{bar}}{RST}"
        f"{DIM}| {{n_fmt}}/{{total_fmt}} [{C}{{elapsed}}{RST}{DIM} < {{remaining}}{RST}"
        f"{DIM}, {{rate_fmt}}]{RST}"
    )

    with tqdm(
        total=remaining,
        bar_format=bar_fmt,
        ncols=70,
        unit="seq",
        dynamic_ncols=True,
    ) as pbar:

        for i, record in enumerate(records, 1):
            seq_id    = record.id
            seq_label = seq_id.split("|")[0]   # just SRA ID for display

            if already_done(seq_id):
                pbar.set_postfix_str(f"{DIM}SKIP {seq_label}{RST}", refresh=True)
                continue

            # ── Pre-query status ─────────────────────────────────────────
            pbar.set_description(f"  {B}BLAST{RST} {seq_label[:22]:<22}")
            pbar.set_postfix_str(f"{Y}submitting...{RST}", refresh=True)

            t0       = time.time()
            xml_data = blast_one(record, pbar)
            elapsed  = time.time() - t0

            if xml_data is None:
                failed += 1
                pbar.write(f"  {R}[✗] {seq_id[:55]} — FAILED{RST}")
                pbar.update(1)
                continue

            if not validate_xml(xml_data):
                pbar.write(f"  {Y}[!] {seq_id[:55]} — unexpected XML format{RST}")

            save_xml(seq_id, xml_data)
            submitted  += 1
            times.append(elapsed)
            n_hits      = count_hits(xml_data)
            total_hits += n_hits

            # ── Per-query result line ─────────────────────────────────────
            hits_str = f"{G}{n_hits} hits{RST}" if n_hits > 0 else f"{Y}0 hits{RST}"
            pbar.write(
                f"  {G}[✓]{RST} {C}{i:>3}/{total}{RST}  "
                f"{seq_label:<14}  {hits_str:<20}  "
                f"{DIM}{len(record.seq):>6} bp  {elapsed:>5.0f}s{RST}"
            )
            pbar.set_postfix_str(f"{G}✓ {n_hits} hits{RST}", refresh=True)
            pbar.update(1)

            # ── Running stats every 5 queries ─────────────────────────────
            if submitted % 5 == 0 and times:
                avg_t     = sum(times) / len(times)
                left      = remaining - submitted
                eta_s     = left * (avg_t + BLAST_DELAY)
                session_t = time.time() - session_start
                pbar.write(
                    f"\n  {hr()}\n"
                    f"  {W}Progress  {RST}{submitted}/{remaining}  "
                    f"({submitted/remaining*100:.0f}%)  "
                    f"| {W}Avg/query{RST} {avg_t:.0f}s  "
                    f"| {W}ETA{RST} {fmt_dur(eta_s)}  "
                    f"| {W}Session{RST} {fmt_dur(session_t)}\n"
                    f"  {hr()}\n"
                )

            # Polite delay
            if submitted + failed < remaining:
                jitter = random.uniform(0, 5)
                time.sleep(BLAST_DELAY + jitter)

    # ── Final summary ──────────────────────────────────────────────────────
    session_total = time.time() - session_start
    avg_hits = total_hits / submitted if submitted else 0

    print()
    print(hr("═"))
    print(f"  {BOLD}{G}BLAST COMPLETE{RST}")
    print(hr("═"))
    print(f"  {W}Submitted   :{RST} {G}{submitted}{RST}")
    print(f"  {W}Failed      :{RST} {(R + str(failed) + RST) if failed else (G + '0' + RST)}")
    print(f"  {W}Total hits  :{RST} {total_hits}")
    print(f"  {W}Avg hits/q  :{RST} {avg_hits:.1f}")
    print(f"  {W}Session time:{RST} {fmt_dur(session_total)}")
    print(f"  {W}XML saved to:{RST} {BLAST_XML_DIR}")
    print(hr("═"))
    print(f"\n  {G}→ Next step:{RST}  python3 scripts/03_parse_blast.py\n")


if __name__ == "__main__":
    main()
