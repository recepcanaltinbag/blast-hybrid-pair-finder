"""
STEP 6 — Professional visualisations for HybridSeqFinder analysis.

Fig 1 — Geographic distribution   : World map (folium interactive + static PNG)
Fig 2 — Query type × organism      : Clustered co-occurrence heatmap
Fig 3 — AMR gene landscape        : Resistance gene × organism (hierarchical clustermap)
Fig 4 — Sequencing landscape      : Platform breakdown + data volume (dual panel)
Fig 5 — BLAST hit quality         : Identity vs coverage + marginal histograms
Fig 6 — Hybrid assembly overview  : Pair count + sequencing volume per organism (lollipop)
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import FINAL_TSV, PAIRED_RUNS, FIGURES_DIR, MIN_PAIRS_PER_QUERY

# ── Colour palette (publication-style) ────────────────────────────────────
DARK_BLUE  = "#1B4F72"
MED_BLUE   = "#2E86C1"
TEAL       = "#148F77"
ORANGE     = "#E67E22"
RED        = "#C0392B"
PURPLE     = "#7D3C98"
GREY       = "#95A5A6"
BG         = "#FAFBFC"
GRID_COLOR = "#E5E8EC"
ACCENT     = DARK_BLUE
ACCENT2    = MED_BLUE

# Nature-style categorical palette
CAT_COLORS = [
    "#1B4F72", "#E67E22", "#148F77", "#C0392B", "#7D3C98",
    "#F39C12", "#1ABC9C", "#8E44AD", "#2980B9", "#E74C3C",
    "#27AE60", "#D35400", "#16A085", "#2C3E50", "#A93226",
]

PLATFORM_COLORS = {
    "OXFORD_NANOPORE": MED_BLUE,
    "PACBIO_SMRT"    : ORANGE,
    "ILLUMINA"       : TEAL,
    "ION_TORRENT"    : PURPLE,
    "UNKNOWN"        : GREY,
}

# Global matplotlib style
sns.set_style("whitegrid")
plt.rcParams.update({
    "font.family"     : "sans-serif",
    "font.size"       : 10,
    "axes.titlesize"  : 12,
    "axes.titleweight": "bold",
    "axes.labelsize"  : 10,
    "xtick.labelsize" : 8,
    "ytick.labelsize" : 8,
})

G   = "\033[92m"; Y = "\033[93m"; R = "\033[91m"
C   = "\033[96m"; W = "\033[97m"; DIM = "\033[2m"; RST = "\033[0m"; BOLD = "\033[1m"

def hr(char="─", width=62):
    return DIM + char * width + RST


def savefig(fig, name, dpi=200):
    os.makedirs(FIGURES_DIR, exist_ok=True)
    path = os.path.join(FIGURES_DIR, name)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"    {G}→{RST} {path}")
    return path


def short_org(s):
    """Return genus + species only."""
    parts = str(s).split()
    return " ".join(parts[:2]) if len(parts) >= 2 else str(s)


def load_data():
    if not os.path.exists(FINAL_TSV):
        print(f"  {R}ERROR: {FINAL_TSV} not found. Run 05_merge_results.py first.{RST}")
        sys.exit(1)
    df = pd.read_csv(FINAL_TSV, sep="\t", low_memory=False)
    print(f"  {G}Loaded{RST} {len(df)} rows, {len(df.columns)} columns")
    return df


# ── Country centroids ──────────────────────────────────────────────────────
COORDS = {
    "China": (35.86, 104.19), "USA": (37.09, -95.71), "United States": (37.09, -95.71),
    "India": (20.59, 78.96), "Germany": (51.17, 10.45), "UK": (55.38, -3.44),
    "United Kingdom": (55.38, -3.44), "France": (46.23, 2.21),
    "Australia": (-25.27, 133.78), "Brazil": (-14.24, -51.93),
    "Japan": (36.20, 138.25), "South Korea": (35.91, 127.77),
    "Italy": (41.87, 12.57), "Spain": (40.46, -3.75),
    "Canada": (56.13, -106.35), "Netherlands": (52.13, 5.29),
    "Sweden": (60.13, 18.64), "Denmark": (56.26, 9.50),
    "Switzerland": (46.82, 8.23), "Belgium": (50.50, 4.47),
    "Turkey": (38.96, 35.24), "Pakistan": (30.38, 69.35),
    "Bangladesh": (23.68, 90.36), "Thailand": (15.87, 100.99),
    "Vietnam": (14.06, 108.28), "Egypt": (26.82, 30.80),
    "Kenya": (-0.02, 37.91), "Nigeria": (9.08, 8.68),
    "South Africa": (-30.56, 22.94), "Argentina": (-38.42, -63.62),
    "Mexico": (23.63, -102.55), "Portugal": (39.40, -8.22),
    "Greece": (39.07, 21.82), "Poland": (51.92, 19.15),
    "Czech Republic": (49.82, 15.47), "Hungary": (47.16, 19.50),
    "Romania": (45.94, 24.97), "Russia": (61.52, 105.32),
    "Norway": (60.47, 8.47), "Finland": (61.92, 25.75),
    "Ireland": (53.41, -8.24), "New Zealand": (-40.90, 174.89),
    "Israel": (31.05, 34.85), "Saudi Arabia": (23.89, 45.08),
    "Iran": (32.43, 53.69), "Singapore": (1.35, 103.82),
    "Malaysia": (4.21, 101.97), "Philippines": (12.88, 121.77),
    "Indonesia": (-0.79, 113.92), "Taiwan": (23.70, 121.00),
    "Austria": (47.52, 14.55), "Colombia": (4.57, -74.30),
    "Peru": (-9.19, -75.02), "Chile": (-35.67, -71.54),
    "Ethiopia": (9.15, 40.49), "Ghana": (7.95, -1.02),
    "Uganda": (1.37, 32.29), "Tanzania": (-6.37, 34.89),
    "Croatia": (45.10, 15.20), "Serbia": (44.02, 21.01),
    "Morocco": (31.79, -7.09), "Tunisia": (33.89, 9.54),
    "Jordan": (30.59, 36.24), "Iraq": (33.22, 43.68),
    "Kuwait": (29.31, 47.49), "Qatar": (25.35, 51.18),
    "UAE": (23.42, 53.85), "Sri Lanka": (7.87, 80.77),
    "Nepal": (28.39, 84.12), "Myanmar": (16.87, 96.19),
    "Cambodia": (12.57, 104.99), "Ukraine": (48.38, 31.17),
    "Kazakhstan": (48.02, 66.92), "Georgia": (42.32, 43.36),
}


# ── Fig 1: Geographic distribution ────────────────────────────────────────

def fig_geo_map(df):
    # ── Static PNG (always produced) ──────────────────────────────────────
    if "country" not in df.columns:
        print(f"    {Y}[skip] No 'country' column.{RST}")
        return

    raw = df.dropna(subset=["country"]).copy()
    raw["country"] = raw["country"].astype(str).str.split(":").str[0].str.strip()
    counts = raw["country"].value_counts().reset_index()
    counts.columns = ["country", "n"]

    # Unique organisms per country
    org_map = {}
    if "organism" in raw.columns:
        org_map = (raw.dropna(subset=["organism"])
                      .groupby("country")["organism"]
                      .apply(lambda x: x.astype(str).apply(short_org).nunique())
                      .to_dict())

    fig, ax = plt.subplots(figsize=(15, 7))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor("#D6EAF8")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.axhline(0, color="#ABEBC6", linewidth=0.5, alpha=0.7)

    max_n = max(counts["n"].max(), 1)
    mapped = 0
    for _, row in counts.iterrows():
        coord = COORDS.get(row["country"])
        if coord is None:
            continue
        lat, lon = coord
        n        = row["n"]
        size     = 40 + (n / max_n) * 450
        alpha    = 0.5 + (n / max_n) * 0.45
        ax.scatter(lon, lat, s=size, c=MED_BLUE, alpha=alpha,
                   edgecolors=DARK_BLUE, linewidths=1.0, zorder=5)
        ax.annotate(str(n), (lon, lat), fontsize=6.5, ha="center",
                    va="center", color="white", fontweight="bold", zorder=6)
        mapped += 1

    ax.set_xlabel("Longitude", fontsize=9, color=ACCENT)
    ax.set_ylabel("Latitude",  fontsize=9, color=ACCENT)
    ax.set_title(
        "Global Distribution of Query Sequence BLAST Hits\n"
        f"(n = {len(counts)} countries, bubble size ∝ NCBI nucleotide hit count)",
        fontsize=12, fontweight="bold", color=ACCENT, pad=10,
    )
    ax.grid(color=GRID_COLOR, linewidth=0.4, alpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Size legend
    for leg_n, leg_lbl in [(1, "1 hit"), (10, "10 hits"), (50, "50 hits")]:
        ax.scatter([], [], s=40 + (leg_n / max_n) * 450,
                   c=MED_BLUE, alpha=0.7, edgecolors=DARK_BLUE, linewidths=1, label=leg_lbl)
    ax.legend(title="NCBI hits", fontsize=8, title_fontsize=9,
              framealpha=0.9, loc="lower left")

    plt.tight_layout()
    savefig(fig, "fig1_geo_map.png")

    # ── Interactive folium HTML ────────────────────────────────────────────
    try:
        import folium
    except ImportError:
        try:
            import subprocess
            subprocess.check_call([sys.executable, "-m", "pip", "install", "folium", "-q"])
            import folium
        except Exception as e:
            print(f"    {Y}[skip folium] {e}{RST}")
            return

    m = folium.Map(location=[25, 15], zoom_start=2, tiles="CartoDB positron")
    for _, row in counts.iterrows():
        coord = COORDS.get(row["country"])
        if coord is None:
            continue
        n       = row["n"]
        n_org   = org_map.get(row["country"], 0)
        radius  = 6 + (n / max_n) * 22
        folium.CircleMarker(
            location=coord, radius=radius,
            color=DARK_BLUE, weight=1.5,
            fill=True, fill_color=MED_BLUE,
            fill_opacity=0.5 + (n / max_n) * 0.45,
            popup=folium.Popup(
                f"<div style='font-family:sans-serif'>"
                f"<b style='color:{DARK_BLUE}'>{row['country']}</b><br>"
                f"BLAST hits: <b>{n}</b><br>"
                f"Unique organisms: <b>{n_org}</b></div>",
                max_width=220
            ),
            tooltip=f"{row['country']}: {n} hits",
        ).add_to(m)

    legend = (
        '<div style="position:fixed;bottom:30px;left:30px;z-index:1000;'
        'background:white;padding:12px 16px;border-radius:8px;'
        'box-shadow:0 2px 6px rgba(0,0,0,0.2);font-family:sans-serif;font-size:12px">'
        f'<b style="color:{DARK_BLUE}">Query Sequence BLAST Hits</b><br>'
        'Circle size ∝ NCBI nucleotide hits<br>'
        '<span style="color:#888;font-size:10px">Source: NCBI SRA + Nucleotide</span>'
        '</div>'
    )
    m.get_root().html.add_child(folium.Element(legend))
    html_path = os.path.join(FIGURES_DIR, "fig1_geo_map.html")
    m.save(html_path)
    print(f"    {G}→{RST} {html_path}")


# ── Fig 2: Query type × organism clustered heatmap ─────────────────────────

def fig_is26_organism_heatmap(df):
    """
    Heatmap: rows = top organisms, columns = query element types.
    """
    needed = ["organism", "is26_type"]
    if any(c not in df.columns for c in needed):
        print(f"    {Y}[skip] Missing columns {[c for c in needed if c not in df.columns]}{RST}")
        return

    sub = df.dropna(subset=needed).copy()
    sub["organism"] = sub["organism"].astype(str).apply(short_org)
    sub["is26_type"] = sub["is26_type"].astype(str).str.replace("_", " ", 1)

    pivot = (sub.groupby(["organism", "is26_type"])
               .size()
               .unstack(fill_value=0))

    top_org = pivot.sum(axis=1).nlargest(20).index
    pivot   = pivot.loc[pivot.index.isin(top_org)]

    if pivot.empty or pivot.shape[0] < 2:
        print(f"    {Y}[skip] Not enough data for query×organism heatmap.{RST}")
        return

    pivot_log = np.log1p(pivot)
    fig_h = max(7, len(pivot) * 0.46)
    fig_w = max(10, len(pivot.columns) * 1.1)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    sns.heatmap(
        pivot_log, ax=ax,
        cmap="Blues",
        linewidths=0.4, linecolor=GRID_COLOR,
        annot=pivot.values,
        fmt="d",
        annot_kws={"size": 8, "color": "black"},
        cbar_kws={"label": "log₁₊(count)", "shrink": 0.65},
    )

    ax.set_title(
        "Query Element Type Distribution across Host Organisms\n"
        "(count annotated; colour = log₁₊ scale)",
        fontsize=12, fontweight="bold", color=ACCENT, pad=14,
    )
    ax.set_xlabel("Query Element Type", fontsize=10, color=ACCENT, labelpad=8)
    ax.set_ylabel("Organism", fontsize=10, color=ACCENT, labelpad=8)
    plt.xticks(rotation=38, ha="right", fontsize=8)
    plt.yticks(fontsize=8, style="italic", rotation=0)
    plt.tight_layout()
    savefig(fig, "fig2_is26_organism_heatmap.png")


# ── Fig 3: AMR gene landscape — hierarchically clustered ──────────────────

def fig_resistance_heatmap(df):
    needed = ["resistance_genes", "organism"]
    if any(c not in df.columns for c in needed):
        print(f"    {Y}[skip] Missing {[c for c in needed if c not in df.columns]}{RST}")
        return

    sub = df.dropna(subset=needed).copy()
    sub = sub[sub["resistance_genes"].astype(str) != "no_AMR"]
    sub["organism"] = sub["organism"].astype(str).apply(short_org)

    rows = []
    for _, row in sub.iterrows():
        for gene in str(row["resistance_genes"]).split(";"):
            gene = gene.strip()
            if gene and gene not in ("nan", "no_AMR", ""):
                rows.append({"organism": row["organism"], "gene": gene})

    if not rows:
        print(f"    {Y}[skip] No AMR gene data.{RST}")
        return

    long  = pd.DataFrame(rows)
    pivot = (long.groupby(["organism", "gene"])
                 .size()
                 .unstack(fill_value=0))

    top_org  = pivot.sum(axis=1).nlargest(20).index
    top_gene = pivot.sum(axis=0).nlargest(25).index
    pivot    = pivot.loc[pivot.index.isin(top_org),
                         pivot.columns.isin(top_gene)]

    if pivot.empty or pivot.shape[0] < 2 or pivot.shape[1] < 2:
        print(f"    {Y}[skip] Pivot too small for heatmap.{RST}")
        return

    pivot_log = np.log1p(pivot)
    annotate  = pivot.shape[1] <= 18

    # Try scipy clustermap (hierarchical clustering)
    try:
        g = sns.clustermap(
            pivot_log,
            cmap="YlOrRd",
            linewidths=0.3,
            linecolor=GRID_COLOR,
            annot=(pivot.values if annotate else False),
            fmt="d",
            annot_kws={"size": 7},
            figsize=(max(14, pivot.shape[1] * 0.70),
                     max(8,  pivot.shape[0] * 0.50)),
            cbar_kws={"label": "log₁₊(count)", "shrink": 0.45},
            yticklabels=True,
            xticklabels=True,
            method="ward",
            metric="euclidean",
        )
        g.fig.patch.set_facecolor(BG)
        g.ax_heatmap.set_title(
            "Resistance Gene × Organism Co-occurrence\n"
            "(hierarchically clustered by Ward linkage, log₁₊ scale)",
            fontsize=12, fontweight="bold", color=ACCENT, pad=10,
        )
        g.ax_heatmap.set_xlabel("Resistance Gene", fontsize=10, color=ACCENT)
        g.ax_heatmap.set_ylabel("Organism",        fontsize=10, color=ACCENT)
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=7)
        plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0,  fontsize=7, style="italic")

        os.makedirs(FIGURES_DIR, exist_ok=True)
        path = os.path.join(FIGURES_DIR, "fig3_amr_heatmap.png")
        g.fig.savefig(path, dpi=200, bbox_inches="tight", facecolor=BG)
        plt.close(g.fig)
        print(f"    {G}→{RST} {path}")

    except Exception as e:
        print(f"    {Y}[clustermap fallback: {e}]{RST}")
        fig, ax = plt.subplots(figsize=(max(14, pivot.shape[1] * 0.70),
                                        max(8,  pivot.shape[0] * 0.50)))
        fig.patch.set_facecolor(BG)
        sns.heatmap(
            pivot_log, ax=ax, cmap="YlOrRd",
            linewidths=0.3, linecolor=GRID_COLOR,
            annot=(pivot.values if annotate else False),
            fmt="d", annot_kws={"size": 7},
            cbar_kws={"label": "log₁₊(count)", "shrink": 0.6},
        )
        ax.set_title("Resistance Gene × Organism Co-occurrence\n(Query Sequences)",
                     fontsize=12, fontweight="bold", color=ACCENT, pad=12)
        ax.set_xlabel("Resistance Gene", fontsize=10, color=ACCENT)
        ax.set_ylabel("Organism",        fontsize=10, color=ACCENT)
        plt.xticks(rotation=45, ha="right", fontsize=7)
        plt.yticks(fontsize=7, style="italic", rotation=0)
        plt.tight_layout()
        savefig(fig, "fig3_amr_heatmap.png")


# ── Fig 4: Sequencing landscape — dual panel ──────────────────────────────

def fig_sequencing_landscape(df, pairs_path):
    """
    Panel A: long-read platform pair counts (bar).
    Panel B: data volume per platform in Gbp, long vs short (grouped bar).
    """
    pairs = pd.DataFrame()
    if os.path.exists(pairs_path):
        try:
            pairs = pd.read_csv(pairs_path, sep="\t")
        except Exception:
            pass

    if pairs.empty or "long_read_platform" not in pairs.columns:
        print(f"    {Y}[skip] No paired_runs data.{RST}")
        return

    for col in ["long_read_bases_Mbp", "short_read_bases_Mbp"]:
        if col in pairs.columns:
            pairs[col] = pd.to_numeric(pairs[col], errors="coerce").fillna(0)

    plat_counts = pairs["long_read_platform"].value_counts().reset_index()
    plat_counts.columns = ["platform", "n_pairs"]

    vol = (pairs.groupby("long_read_platform")
               .agg(
                   long_Gbp =("long_read_bases_Mbp",  lambda x: x.sum() / 1000),
                   short_Gbp=("short_read_bases_Mbp", lambda x: x.sum() / 1000),
               )
               .reset_index()
               .rename(columns={"long_read_platform": "platform"}))

    fig, axes = plt.subplots(1, 2, figsize=(13, max(4, len(plat_counts) * 1.2)))
    fig.patch.set_facecolor(BG)

    # ── Panel A: pair counts ────────────────────────────────────────────
    ax = axes[0]
    ax.set_facecolor(BG)
    colors_a = [PLATFORM_COLORS.get(p, GREY) for p in plat_counts["platform"]]
    bars = ax.barh(plat_counts["platform"], plat_counts["n_pairs"],
                   color=colors_a, edgecolor="white", height=0.55)
    for bar, val in zip(bars, plat_counts["n_pairs"]):
        ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", ha="left",
                fontsize=11, fontweight="bold", color=ACCENT)
    ax.set_xlabel("Number of long + short read pairs", fontsize=10, color=ACCENT)
    ax.set_title("A.  Hybrid Pair Count\nby Long-Read Platform",
                 fontsize=11, fontweight="bold", color=ACCENT)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color=GRID_COLOR, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.set_yticklabels([p.replace("_", " ").title() for p in plat_counts["platform"]],
                       fontsize=9)

    # ── Panel B: data volume ────────────────────────────────────────────
    ax2 = axes[1]
    ax2.set_facecolor(BG)
    x  = np.arange(len(vol))
    w  = 0.38
    c1 = [PLATFORM_COLORS.get(p, GREY) for p in vol["platform"]]

    bars_l = ax2.bar(x - w / 2, vol["long_Gbp"],  w,
                     color=c1, edgecolor="white", alpha=0.9, label="Long-read")
    bars_s = ax2.bar(x + w / 2, vol["short_Gbp"], w,
                     color=TEAL, edgecolor="white", alpha=0.75, label="Short-read")

    for bar in list(bars_l) + list(bars_s):
        h = bar.get_height()
        if h > 0.01:
            ax2.text(bar.get_x() + bar.get_width() / 2, h + 0.5,
                     f"{h:.1f}", ha="center", va="bottom", fontsize=7.5, color=ACCENT)

    ax2.set_xticks(x)
    ax2.set_xticklabels([p.replace("_", "\n") for p in vol["platform"]], fontsize=9)
    ax2.set_ylabel("Data Volume (Gbp)", fontsize=10, color=ACCENT)
    ax2.set_title("B.  Total Sequencing Volume\n(Long vs Short Read, Gbp)",
                  fontsize=11, fontweight="bold", color=ACCENT)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(axis="y", color=GRID_COLOR, linewidth=0.7)
    ax2.set_axisbelow(True)
    ax2.legend(fontsize=9, framealpha=0.9)

    fig.suptitle(
        "Sequencing Landscape of Query-linked BioSamples",
        fontsize=13, fontweight="bold", color=ACCENT, y=1.03,
    )
    plt.tight_layout()
    savefig(fig, "fig4_sequencing_landscape.png")


# ── Fig 5: BLAST hit quality + marginal histograms ────────────────────────

def fig_blast_quality(df):
    """
    Main scatter: identity_pct vs query_cov_pct.
    Point size ∝ query_len (shows length → coverage relationship).
    Marginal histograms via GridSpec.
    """
    needed = ["identity_pct", "query_cov_pct"]
    if any(c not in df.columns for c in needed):
        print(f"    {Y}[skip] Missing {needed}{RST}")
        return

    sub = df.dropna(subset=needed).copy()
    if sub.empty:
        print(f"    {Y}[skip] No BLAST quality data.{RST}")
        return

    sub["query_len"] = pd.to_numeric(sub.get("query_len", 0), errors="coerce").fillna(0)
    sub["org_short"] = sub["organism"].fillna("Unknown").astype(str).apply(short_org)

    top_orgs  = sub["org_short"].value_counts().head(8).index.tolist()
    color_map = {o: CAT_COLORS[i % len(CAT_COLORS)]
                 for i, o in enumerate(sorted(top_orgs))}
    color_map["Other"] = GREY

    sub["org_label"] = sub["org_short"].apply(lambda x: x if x in top_orgs else "Other")
    point_colors     = [color_map[o] for o in sub["org_label"]]
    point_sizes      = np.clip(sub["query_len"] / 180, 12, 220)

    fig = plt.figure(figsize=(11, 8), facecolor=BG)
    gs  = gridspec.GridSpec(
        2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
        hspace=0.06, wspace=0.06,
    )
    ax_main  = fig.add_subplot(gs[1, 0])
    ax_top   = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)

    for ax in (ax_main, ax_top, ax_right):
        ax.set_facecolor(BG)

    # Main scatter
    ax_main.scatter(
        sub["query_cov_pct"], sub["identity_pct"],
        c=point_colors, s=point_sizes,
        alpha=0.68, edgecolors="white", linewidths=0.4, zorder=3,
    )

    # Threshold lines
    ax_main.axvline(50, color=ORANGE, linestyle="--", linewidth=1.3,
                    alpha=0.85, zorder=2, label="50% cov. threshold")
    ax_main.axhline(80, color=RED,    linestyle="--", linewidth=1.3,
                    alpha=0.85, zorder=2, label="80% identity threshold")
    ax_main.fill_between([50, 106], 80, 106,
                         color=TEAL, alpha=0.07, zorder=1, label="Passing zone")

    ax_main.set_xlabel("Query Coverage (%)", fontsize=11, color=ACCENT, labelpad=6)
    ax_main.set_ylabel("Identity (%)",       fontsize=11, color=ACCENT, labelpad=6)
    ax_main.set_xlim(0, 106)
    ax_main.set_ylim(60, 106)
    ax_main.grid(color=GRID_COLOR, linewidth=0.5)
    ax_main.spines["top"].set_visible(False)
    ax_main.spines["right"].set_visible(False)

    # Organism legend
    org_handles = [mpatches.Patch(color=color_map[o], label=o)
                   for o in sorted(top_orgs)] + \
                  [mpatches.Patch(color=GREY, label="Other")]
    # Threshold legend
    thr_handles = [
        plt.Line2D([0], [0], color=ORANGE, linestyle="--", linewidth=1.3, label="50% cov. cutoff"),
        plt.Line2D([0], [0], color=RED,    linestyle="--", linewidth=1.3, label="80% identity cutoff"),
        mpatches.Patch(color=TEAL, alpha=0.3, label="Passing zone"),
    ]
    leg1 = ax_main.legend(handles=org_handles, title="Organism",
                          fontsize=7, title_fontsize=8,
                          framealpha=0.9, loc="lower right",
                          bbox_to_anchor=(1.0, 0.0))
    ax_main.add_artist(leg1)
    ax_main.legend(handles=thr_handles, fontsize=8, framealpha=0.9,
                   loc="upper left")

    # Size annotation
    ax_main.text(0.01, 0.01, "Point size ∝ query length",
                 transform=ax_main.transAxes, fontsize=7.5,
                 color=GREY, va="bottom")

    # Top histogram (coverage)
    ax_top.hist(sub["query_cov_pct"], bins=35, color=MED_BLUE,
                alpha=0.85, edgecolor="white", linewidth=0.4)
    ax_top.axvline(50, color=ORANGE, linestyle="--", linewidth=1, alpha=0.7)
    ax_top.set_ylabel("Count", fontsize=8, color=ACCENT)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.grid(color=GRID_COLOR, linewidth=0.4)
    plt.setp(ax_top.get_xticklabels(), visible=False)
    ax_top.set_title(
        "BLAST Hit Quality — Query Sequences\n"
        "(point size ∝ query length)",
        fontsize=11, fontweight="bold", color=ACCENT, pad=8,
    )

    # Right histogram (identity)
    ax_right.hist(sub["identity_pct"], bins=35, color=TEAL,
                  alpha=0.85, edgecolor="white", linewidth=0.4,
                  orientation="horizontal")
    ax_right.axhline(80, color=RED, linestyle="--", linewidth=1, alpha=0.7)
    ax_right.set_xlabel("Count", fontsize=8, color=ACCENT)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.grid(color=GRID_COLOR, linewidth=0.4)
    plt.setp(ax_right.get_yticklabels(), visible=False)

    savefig(fig, "fig5_blast_quality.png")


# ── Fig 6: Hybrid assembly overview — lollipop dual panel ─────────────────

def fig_hybrid_overview(pairs_path):
    """
    Panel A: pair count per organism (lollipop).
    Panel B: long + short read data volume per organism (lollipop).
    """
    if not os.path.exists(pairs_path):
        print(f"    {Y}[skip] paired_runs.tsv not found.{RST}")
        return

    try:
        pairs = pd.read_csv(pairs_path, sep="\t")
    except Exception as e:
        print(f"    {Y}[skip] {e}{RST}")
        return

    if pairs.empty:
        print(f"    {Y}[skip] No paired runs.{RST}")
        return

    # Resolve organism name
    if "organism_sra" in pairs.columns:
        org_raw = pairs["organism_sra"].fillna(
            pairs.get("organism_blast", pd.Series(["Unknown"] * len(pairs))))
    elif "organism_blast" in pairs.columns:
        org_raw = pairs["organism_blast"]
    else:
        org_raw = pd.Series(["Unknown"] * len(pairs))
    pairs["org_short"] = org_raw.fillna("Unknown").astype(str).apply(short_org)

    for col in ["long_read_bases_Mbp", "short_read_bases_Mbp"]:
        if col in pairs.columns:
            pairs[col] = pd.to_numeric(pairs[col], errors="coerce").fillna(0)

    agg = (pairs.groupby("org_short")
               .agg(
                   n_pairs  =("long_read_run",        "count"),
                   long_Gbp =("long_read_bases_Mbp",  lambda x: x.sum() / 1000),
                   short_Gbp=("short_read_bases_Mbp", lambda x: x.sum() / 1000),
               )
               .reset_index()
               .sort_values("n_pairs"))

    top = agg.tail(20).reset_index(drop=True)

    if top.empty:
        print(f"    {Y}[skip] No organism data.{RST}")
        return

    y = np.arange(len(top))
    fig, axes = plt.subplots(1, 2, figsize=(14, max(6, len(top) * 0.48)),
                             sharey=True)
    fig.patch.set_facecolor(BG)

    # ── Panel A: pair count lollipop ──────────────────────────────────
    ax = axes[0]
    ax.set_facecolor(BG)
    ax.hlines(y, 0, top["n_pairs"], color=GRID_COLOR, linewidth=1.8, zorder=1)
    ax.scatter(top["n_pairs"], y, color=MED_BLUE, s=90, zorder=3,
               edgecolors=DARK_BLUE, linewidths=1.0)
    for i, val in enumerate(top["n_pairs"]):
        ax.text(val + 0.5, i, str(val), va="center", ha="left",
                fontsize=7.5, color=ACCENT)
    ax.axvline(MIN_PAIRS_PER_QUERY, color=RED, linestyle="--",
               linewidth=1.2, alpha=0.75,
               label=f"≥{MIN_PAIRS_PER_QUERY} pairs target")
    ax.set_yticks(y)
    ax.set_yticklabels(top["org_short"], fontsize=8, style="italic")
    ax.set_xlabel("Number of hybrid pairs", fontsize=10, color=ACCENT)
    ax.set_title("A.  Hybrid Pair Count per Organism",
                 fontsize=11, fontweight="bold", color=ACCENT)
    ax.legend(fontsize=8, framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color=GRID_COLOR, linewidth=0.5)
    ax.set_axisbelow(True)

    # ── Panel B: data volume lollipop ─────────────────────────────────
    ax2 = axes[1]
    ax2.set_facecolor(BG)
    offset = 0.20

    # Long read
    ax2.hlines(y + offset, 0, top["long_Gbp"],
               color=GRID_COLOR, linewidth=1.4, zorder=1)
    ax2.scatter(top["long_Gbp"], y + offset,
                color=MED_BLUE, s=75, zorder=3,
                edgecolors=DARK_BLUE, linewidths=0.9,
                label="Long-read (Nanopore/PacBio)")

    # Short read
    ax2.hlines(y - offset, 0, top["short_Gbp"],
               color=GRID_COLOR, linewidth=1.4, zorder=1)
    ax2.scatter(top["short_Gbp"], y - offset,
                color=TEAL, s=75, zorder=3,
                edgecolors="#0E6655", linewidths=0.9,
                label="Short-read (Illumina)")

    ax2.set_xlabel("Data Volume (Gbp)", fontsize=10, color=ACCENT)
    ax2.set_title("B.  Sequencing Volume per Organism\n(Long vs Short Read)",
                  fontsize=11, fontweight="bold", color=ACCENT)
    ax2.legend(fontsize=8, framealpha=0.9, loc="lower right")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(axis="x", color=GRID_COLOR, linewidth=0.5)
    ax2.set_axisbelow(True)

    fig.suptitle(
        "Hybrid Assembly Potential — Query Sequences\n"
        f"(top {len(top)} organisms by hybrid pair count)",
        fontsize=13, fontweight="bold", color=ACCENT, y=1.01,
    )
    plt.tight_layout()
    savefig(fig, "fig6_hybrid_overview.png")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    print()
    print(hr("═"))
    print(f"  {BOLD}{C}STEP 06 — VISUALISATIONS{RST}")
    print(hr("═"))

    df = load_data()
    print()

    figures = [
        ("Fig 1 — Geographic distribution",   lambda: fig_geo_map(df)),
        ("Fig 2 — Query type × organism",       lambda: fig_is26_organism_heatmap(df)),
        ("Fig 3 — AMR gene landscape",         lambda: fig_resistance_heatmap(df)),
        ("Fig 4 — Sequencing landscape",       lambda: fig_sequencing_landscape(df, PAIRED_RUNS)),
        ("Fig 5 — BLAST hit quality",          lambda: fig_blast_quality(df)),
        ("Fig 6 — Hybrid assembly overview",   lambda: fig_hybrid_overview(PAIRED_RUNS)),
    ]

    ok = 0
    for label, fn in figures:
        print(f"  {C}{label}{RST}")
        try:
            fn()
            ok += 1
        except Exception as e:
            print(f"    {R}[ERROR] {e}{RST}")
        print()

    print(hr("═"))
    print(f"  {BOLD}{G}VISUALISATIONS COMPLETE{RST}")
    print(hr("═"))
    print(f"  {W}Figures produced :{RST} {ok} / {len(figures)}")
    print(f"  {W}Output directory :{RST} {FIGURES_DIR}")
    print(hr("═"))
    print()


if __name__ == "__main__":
    main()
