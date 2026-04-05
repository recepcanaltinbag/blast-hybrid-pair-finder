"""
STEP 4 — SRA metadata + long/short read pairing.

For each BLAST hit accession:
  1. accession → nucleotide UID (esearch)
  2. nuc UID   → BioSample UID (elink)
  3. BioSample UID → SAMN accession (efetch)
  4. SAMN → all SRA runs (esearch + efetch runinfo)
  5. Check long+short pairs on that BioSample

Outputs:
  results/sra_metadata.tsv
  results/paired_runs.tsv
"""

import os, sys, time, xml.etree.ElementTree as ET
from collections import defaultdict
from typing import List, Dict

import pandas as pd
from Bio import Entrez
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    ENTREZ_EMAIL, ENTREZ_API_KEY,
    BLAST_HITS, SRA_META, PAIRED_RUNS,
    LONG_READ_PLATFORMS, SHORT_READ_PLATFORMS,
    MIN_PAIRS_PER_QUERY,
)

Entrez.email = ENTREZ_EMAIL
if ENTREZ_API_KEY:
    Entrez.api_key = ENTREZ_API_KEY

_DELAY = 0.35 if ENTREZ_API_KEY else 0.4

G   = "\033[92m"; Y = "\033[93m"; R = "\033[91m"
C   = "\033[96m"; W = "\033[97m"; DIM = "\033[2m"; RST = "\033[0m"; BOLD = "\033[1m"

def hr(char="─", width=62):
    return DIM + char * width + RST


# ── Entrez helper ──────────────────────────────────────────────────────────

def ecall(fn, *args, retries=4, **kwargs):
    for attempt in range(1, retries + 1):
        try:
            h    = fn(*args, **kwargs)
            data = h.read()
            h.close()
            time.sleep(_DELAY)
            return data
        except Exception as e:
            wait = 2 ** attempt
            tqdm.write(f"  {Y}[!] {fn.__name__} attempt {attempt}: {e} — wait {wait}s{RST}")
            time.sleep(wait)
    return None


# ── Step functions ─────────────────────────────────────────────────────────

def acc_to_nuc_uid(acc: str) -> str:
    data = ecall(Entrez.esearch, db="nucleotide",
                 term=f"{acc}[Accession]", retmax=3)
    if not data:
        return ""
    root = ET.fromstring(data)
    ids  = [el.text for el in root.iter("Id") if el.text]
    return ids[0] if ids else ""


def nuc_uid_to_samn(nuc_uid: str) -> str:
    # elink: nucleotide → biosample
    data = ecall(Entrez.elink, dbfrom="nucleotide",
                 db="biosample", id=nuc_uid)
    if not data:
        return ""
    root    = ET.fromstring(data)
    bs_uids = [el.text for el in root.findall(".//LinkSetDb/Link/Id") if el.text]
    if not bs_uids:
        return ""
    # efetch to get SAMN accession
    fetch = ecall(Entrez.efetch, db="biosample",
                  id=",".join(bs_uids[:5]), retmode="xml")
    if not fetch:
        return ""
    try:
        r = ET.fromstring(fetch)
        for bs in r.iter("BioSample"):
            acc = bs.get("accession", "")
            if acc.startswith(("SAMN", "SAMD", "SAME")):
                return acc
    except ET.ParseError:
        pass
    return ""


def samn_to_runs(samn: str) -> List[Dict]:
    # esearch SRA by BioSample
    data = ecall(Entrez.esearch, db="sra",
                 term=f"{samn}[BioSample]", retmax=200)
    if not data:
        return []
    root    = ET.fromstring(data)
    sra_ids = [el.text for el in root.iter("Id") if el.text]
    if not sra_ids:
        return []
    # efetch runinfo
    fetch = ecall(Entrez.efetch, db="sra", id=",".join(sra_ids),
                  rettype="runinfo", retmode="xml")
    if not fetch:
        return []
    runs = []
    try:
        root2 = ET.fromstring(fetch)
        for row in root2.iter("Row"):
            r       = {c.tag: (c.text or "").strip() for c in row}
            acc     = r.get("Run", "")
            if not acc:
                continue
            raw     = r.get("Platform", "") + " " + r.get("Model", "")
            plat    = normalise_platform(raw)
            bases   = r.get("bases", "")
            try:
                bases_mb = round(int(bases) / 1e6, 2)
            except (ValueError, TypeError):
                bases_mb = ""
            runs.append({
                "sra_run"          : acc,
                "platform"         : plat,
                "platform_raw"     : r.get("Platform", ""),
                "model"            : r.get("Model", ""),
                "library_strategy" : r.get("LibraryStrategy", ""),
                "library_layout"   : r.get("LibraryLayout", ""),
                "spots"            : r.get("spots", ""),
                "bases_Mbp"        : bases_mb,
                "biosample"        : r.get("BioSample", "") or samn,
                "bioproject"       : r.get("BioProject", ""),
                "organism_sra"     : r.get("ScientificName", ""),
                "taxid"            : r.get("TaxID", ""),
                "center"           : r.get("CenterName", ""),
                "sra_study"        : r.get("SRAStudy", ""),
            })
    except ET.ParseError:
        pass
    return runs


def normalise_platform(raw: str) -> str:
    raw = raw.upper()
    if any(k in raw for k in ("NANOPORE","MINION","GRIDION","PROMETHION")):
        return "OXFORD_NANOPORE"
    if any(k in raw for k in ("PACBIO","SEQUEL","RSII","HIFI")):
        return "PACBIO_SMRT"
    if any(k in raw for k in ("ILLUMINA","NEXTSEQ","NOVASEQ","HISEQ","MISEQ")):
        return "ILLUMINA"
    if "ION" in raw and "TORRENT" in raw:
        return "ION_TORRENT"
    return raw.strip() or "UNKNOWN"


def get_biosample_meta(samn: str) -> dict:
    data = ecall(Entrez.efetch, db="biosample", id=samn, retmode="xml")
    if not data:
        return {}
    meta = {}
    try:
        root = ET.fromstring(data)
        for attr in root.iter("Attribute"):
            name = (attr.get("harmonized_name") or attr.get("attribute_name","")).lower()
            meta[name] = (attr.text or "").strip()
    except ET.ParseError:
        pass
    return meta


def make_pairs(runs: List[Dict]) -> List[Dict]:
    by_bs = defaultdict(lambda: {"long": [], "short": []})
    for r in runs:
        bs = r.get("biosample") or "__"
        if r["platform"] in LONG_READ_PLATFORMS:
            by_bs[bs]["long"].append(r)
        elif r["platform"] in SHORT_READ_PLATFORMS:
            by_bs[bs]["short"].append(r)

    pairs = []
    for bs, bk in by_bs.items():
        if not bk["long"] or not bk["short"]:
            continue
        def best(lst):
            try:    return max(lst, key=lambda x: float(x.get("bases_Mbp") or 0))
            except: return lst[0]
        lr, sr = best(bk["long"]), best(bk["short"])
        pairs.append({
            "biosample"            : bs if bs != "__" else "",
            "bioproject"           : lr.get("bioproject") or sr.get("bioproject",""),
            "organism_sra"         : lr.get("organism_sra") or sr.get("organism_sra",""),
            "long_read_run"        : lr["sra_run"],
            "long_read_platform"   : lr["platform"],
            "long_read_strategy"   : lr.get("library_strategy",""),
            "long_read_bases_Mbp"  : lr.get("bases_Mbp",""),
            "long_read_spots"      : lr.get("spots",""),
            "short_read_run"       : sr["sra_run"],
            "short_read_platform"  : sr["platform"],
            "short_read_strategy"  : sr.get("library_strategy",""),
            "short_read_layout"    : sr.get("library_layout",""),
            "short_read_bases_Mbp" : sr.get("bases_Mbp",""),
            "short_read_spots"     : sr.get("spots",""),
            "n_long_runs"          : len(bk["long"]),
            "n_short_runs"         : len(bk["short"]),
        })
    return pairs


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    if not os.path.exists(BLAST_HITS):
        print(f"{R}ERROR: {BLAST_HITS} not found.{RST}"); sys.exit(1)

    hits_df    = pd.read_csv(BLAST_HITS, sep="\t")
    accessions = hits_df["hit_accession"].dropna().unique().tolist()

    # Resume
    done_acc  = set()
    all_runs  = []
    all_pairs = []
    if os.path.exists(SRA_META) and os.path.getsize(SRA_META) > 0:
        try:
            prev     = pd.read_csv(SRA_META, sep="\t")
            done_acc = set(prev["hit_accession"].dropna().unique())
            all_runs = prev.to_dict("records")
        except Exception:
            pass
    if os.path.exists(PAIRED_RUNS) and os.path.getsize(PAIRED_RUNS) > 0:
        try:
            all_pairs = pd.read_csv(PAIRED_RUNS, sep="\t").to_dict("records")
        except Exception:
            pass

    remaining = [a for a in accessions if a not in done_acc]

    print()
    print(hr("═"))
    print(f"  {BOLD}{C}STEP 04 — SRA METADATA + PAIRING{RST}")
    print(hr("═"))
    print(f"  {W}Total accessions :{RST} {len(accessions)}")
    print(f"  {W}Already cached   :{RST} {len(done_acc)}")
    print(f"  {Y}Remaining        :{RST} {len(remaining)}")
    print(f"  {W}Strategy         :{RST} BioSample-direct (no BioProject scan)")
    print(hr("═"))
    print()

    stats = {"pairs": 0, "long": 0, "short": 0, "no_samn": 0, "no_sra": 0}

    bar_fmt = (f"  {C}{{l_bar}}{{bar}}{RST}"
               f"{DIM}| {{n_fmt}}/{{total_fmt}} "
               f"[{{elapsed}}<{{remaining}}, {{rate_fmt}}]{RST}")

    with tqdm(remaining, bar_format=bar_fmt, ncols=70, unit="acc") as pbar:
        for acc in pbar:
            org = hits_df.loc[hits_df["hit_accession"] == acc, "organism"].iloc[0]
            pbar.set_description(f"  {C}{acc[:20]:<20}{RST}")

            # 1. acc → nuc UID
            nuc_uid = acc_to_nuc_uid(acc)
            if not nuc_uid:
                tqdm.write(f"  {DIM}[–] {acc:<22} no nuc UID{RST}")
                all_runs.append({"hit_accession": acc, "organism_blast": org})
                stats["no_samn"] += 1
                continue

            # 2. nuc UID → SAMN
            samn = nuc_uid_to_samn(nuc_uid)
            if not samn:
                tqdm.write(f"  {DIM}[–] {acc:<22} no BioSample{RST}")
                all_runs.append({"hit_accession": acc, "organism_blast": org})
                stats["no_samn"] += 1
                continue

            # 3. SAMN → SRA runs
            runs = samn_to_runs(samn)
            for r in runs:
                r["hit_accession"]  = acc
                r["organism_blast"] = org
                r["biosample"]      = r.get("biosample") or samn
            all_runs.extend(runs)

            if not runs:
                tqdm.write(f"  {DIM}[–] {acc:<22} {samn} — no SRA runs{RST}")
                all_runs.append({"hit_accession": acc, "organism_blast": org,
                                 "biosample": samn})
                stats["no_sra"] += 1
                continue

            # 4. pairs?
            pairs = make_pairs(runs)
            lr_n  = sum(1 for r in runs if r["platform"] in LONG_READ_PLATFORMS)
            sr_n  = sum(1 for r in runs if r["platform"] in SHORT_READ_PLATFORMS)

            if pairs:
                meta = get_biosample_meta(samn)
                for p in pairs:
                    p["hit_accession"]    = acc
                    p["organism_blast"]   = org
                    p["isolation_source"] = meta.get("isolation_source","")
                    p["host"]             = meta.get("host","")
                    p["collection_date"]  = meta.get("collection_date","")
                    p["country"]          = meta.get("geo_loc_name", meta.get("country",""))
                    p["strain"]           = meta.get("strain","")
                all_pairs.extend(pairs)
                stats["pairs"] += len(pairs)
                tqdm.write(
                    f"  {G}[✓]{RST} {acc:<22} {DIM}{samn}{RST}  "
                    f"{C}long:{lr_n}{RST}  short:{sr_n}  "
                    f"{G}{len(pairs)} pair(s){RST}"
                )
            else:
                reason = "no long" if lr_n == 0 else "no short"
                tqdm.write(
                    f"  {Y}[○]{RST} {acc:<22} {DIM}{samn}{RST}  "
                    f"long:{lr_n}  short:{sr_n}  {Y}{reason}{RST}"
                )

            stats["long"]  += lr_n
            stats["short"] += sr_n

            pbar.set_postfix_str(
                f"pairs:{stats['pairs']} long:{stats['long']}",
                refresh=True,
            )

            # Incremental save every 25
            if (len(done_acc) + pbar.n) % 25 == 0:
                pd.DataFrame(all_runs).to_csv(SRA_META,     sep="\t", index=False)
                pd.DataFrame(all_pairs).to_csv(PAIRED_RUNS, sep="\t", index=False)

    # Final save
    pd.DataFrame(all_runs).to_csv(SRA_META,     sep="\t", index=False)
    pd.DataFrame(all_pairs).to_csv(PAIRED_RUNS, sep="\t", index=False)

    pairs_df = pd.DataFrame(all_pairs)
    print()
    print(hr("═"))
    print(f"  {BOLD}{G}SRA METADATA COMPLETE{RST}")
    print(hr("═"))
    print(f"  {W}Processed        :{RST} {len(remaining)}")
    print(f"  {W}No BioSample     :{RST} {Y}{stats['no_samn']}{RST}")
    print(f"  {W}No SRA runs      :{RST} {Y}{stats['no_sra']}{RST}")
    print(f"  {W}Long-read runs   :{RST} {G}{stats['long']}{RST}")
    print(f"  {W}Short-read runs  :{RST} {stats['short']}")
    print(f"  {W}Total pairs      :{RST} {G}{stats['pairs']}{RST}")

    if not pairs_df.empty and "hit_accession" in pairs_df.columns:
        if "long_read_platform" in pairs_df.columns:
            print(hr())
            print(f"  {W}Long-read platforms:{RST}")
            vc = pairs_df["long_read_platform"].value_counts()
            for plat, cnt in vc.items():
                b = "█" * min(int(cnt / vc.max() * 20), 20)
                print(f"    {C}{b:<20}{RST} {cnt:>4}  {DIM}{plat}{RST}")

        print(hr())
        print(f"  {W}Top organisms with pairs:{RST}")
        orgs = pairs_df["organism_sra"].fillna(pairs_df["organism_blast"]).value_counts().head(8)
        for org, cnt in orgs.items():
            b = "█" * min(int(cnt / orgs.max() * 20), 20)
            print(f"    {G}{b:<20}{RST} {cnt:>3}  {DIM}{str(org)[:35]}{RST}")

    print(hr("═"))
    print(f"\n  {G}→ Saved:{RST} {SRA_META}")
    print(f"  {G}→ Saved:{RST} {PAIRED_RUNS}")
    print(f"  {G}→ Next :{RST} python3 scripts/05_merge_results.py\n")


if __name__ == "__main__":
    main()
