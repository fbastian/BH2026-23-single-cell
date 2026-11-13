"""Step 02c - QC v1 thresholds (MASLD final-QC recipe, 2026-10-29, ported to mouse) + ridge figures BEFORE removal.

Strata: sample x pre-QC L0 (02b); strata < 100 cells take the sample-level fit. MAD = unscaled (as md_02).
  log1p n_counts, log1p n_genes : median +/- 5 MAD; immune: no upper bound (MASLD D1)
  pct_mito (13 protein-coding mt-): median + 3 MAD, floored at nuclei 5 % / cells 10 % / epithelial cells 20 % (D5)
  pct_hb : > 20 % removed (D2)
Comparison rule (authors, 3a_SCT.R): per sample, all lineages pooled, 3 scaled MAD (1.4826) on log genes / log UMI
both sides, pct_mito upper (their mito = all ^mt- genes; here the 13 protein-coding, the same set in mouse).
Doublets are not in this step (scDblFinder + Scrublet consensus comes next).

Input : <OUT>/02_qc_v1/cells_preqc.parquet
Output: <OUT>/02_qc_v1/thresholds_v1.tsv, cells_flags_v1.parquet, qc_summary_sample_L0.tsv
        <OUT>/02_qc_v1/figures/QC_ridge_<cell|nucleus>_<metric>.{pdf,png}
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/masld_vault/cross_species/guilliams2022"
D = f"{OUT}/02_qc_v1"
FIG = f"{D}/figures"
os.makedirs(FIG, exist_ok=True)
LIN = ["immune", "epithelial", "endothelial", "stromal"]
PAL = {"immune": "#56B4E9", "epithelial": "#009E73", "endothelial": "#CC79A7", "stromal": "#E69F00"}
MIN_CELLS = 100
c = pd.read_parquet(f"{D}/cells_preqc.parquet")
c["L0"] = c["L0"].astype(str)


def mad(x):
    x = np.asarray(x, float)
    med = np.median(x)
    return med, np.median(np.abs(x - med))


# ---- v1 thresholds per sample x L0 ------------------------------------------------------------------------------
rows = []
for (s, L), g in c.groupby(["sample", "L0"], observed=True):
    fit = g if len(g) >= MIN_CELLS else c[c["sample"] == s]
    susp = g["suspension"].iloc[0]
    r = dict(sample=s, L0=L, suspension=susp, n=len(g), fit="stratum" if len(g) >= MIN_CELLS else "sample")
    for met in ["log1p_n_counts", "log1p_n_genes"]:
        med, md = mad(fit[met])
        r[f"{met}_lo"] = med - 5 * md
        r[f"{met}_hi"] = np.inf if L == "immune" else med + 5 * md
    med, md = mad(fit["pct_mito"])
    floor = 5.0 if susp == "nucleus" else (20.0 if L == "epithelial" else 10.0)
    r["pct_mito_hi"] = max(med + 3 * md, floor)
    r["pct_hb_hi"] = 20.0
    rows.append(r)
T = pd.DataFrame(rows)
T.to_csv(f"{D}/thresholds_v1.tsv", sep="\t", index=False)

# ---- authors' rule per sample -----------------------------------------------------------------------------------
A = []
for s, g in c.groupby("sample", observed=True):
    r = dict(sample=s)
    for met in ["log1p_n_counts", "log1p_n_genes"]:
        med, md = mad(g[met])
        r[f"{met}_lo"], r[f"{met}_hi"] = med - 3 * 1.4826 * md, med + 3 * 1.4826 * md
    med, md = mad(g["pct_mito"])
    r["pct_mito_hi"] = med + 3 * 1.4826 * md
    A.append(r)
A = pd.DataFrame(A).set_index("sample")

# ---- flags --------------------------------------------------------------------------------------------------------
t = c[["sample", "L0"]].merge(T, on=["sample", "L0"], how="left").set_index(c.index)
F = pd.DataFrame(index=c.index)
F["v1_counts"] = (c["log1p_n_counts"] < t["log1p_n_counts_lo"]) | (c["log1p_n_counts"] > t["log1p_n_counts_hi"])
F["v1_genes"] = (c["log1p_n_genes"] < t["log1p_n_genes_lo"]) | (c["log1p_n_genes"] > t["log1p_n_genes_hi"])
F["v1_mito"] = c["pct_mito"] > t["pct_mito_hi"]
F["v1_hb"] = c["pct_hb"] > t["pct_hb_hi"]
F["v1_fail"] = F[["v1_counts", "v1_genes", "v1_mito", "v1_hb"]].any(axis=1)
a = A.loc[c["sample"]].set_index(c.index)
F["auth_fail"] = ((c["log1p_n_counts"] < a["log1p_n_counts_lo"]) | (c["log1p_n_counts"] > a["log1p_n_counts_hi"])
                  | (c["log1p_n_genes"] < a["log1p_n_genes_lo"]) | (c["log1p_n_genes"] > a["log1p_n_genes_hi"])
                  | (c["pct_mito"] > a["pct_mito_hi"]))
F.join(c).to_parquet(f"{D}/cells_flags_v1.parquet")
S = F.join(c[["sample", "L0", "digestion", "suspension"]]).groupby(["sample", "L0"], observed=True).agg(
    digestion=("digestion", "first"), n=("v1_fail", "size"), v1_counts=("v1_counts", "sum"),
    v1_genes=("v1_genes", "sum"), v1_mito=("v1_mito", "sum"), v1_hb=("v1_hb", "sum"), v1_fail=("v1_fail", "sum"),
    authors_fail=("auth_fail", "sum"))
S.to_csv(f"{D}/qc_summary_sample_L0.tsv", sep="\t")
tot = F.join(c[["suspension", "L0"]]).groupby(["suspension", "L0"])[["v1_fail", "auth_fail"]].agg(["sum", "size"])
print(tot.to_string())
print(f"v1 removes {F['v1_fail'].sum()} of {len(F)} ({100 * F['v1_fail'].mean():.1f} %); "
      f"authors' rule {F['auth_fail'].sum()} ({100 * F['auth_fail'].mean():.1f} %)", flush=True)

# ---- ridge figures ------------------------------------------------------------------------------------------------
HALO = [pe.withStroke(linewidth=3, foreground="white")]
METRICS = {"log1p_n_counts": "log1p #UMI", "log1p_n_genes": "log1p #genes", "pct_mito": "% mito (13 mt- genes)",
           "pct_hb": "% haemoglobin", "pct_top20": "% top-20 genes", "complexity": "log10 genes / log10 UMI"}
V1KEY = {"log1p_n_counts": ("log1p_n_counts_lo", "log1p_n_counts_hi"),
         "log1p_n_genes": ("log1p_n_genes_lo", "log1p_n_genes_hi"),
         "pct_mito": (None, "pct_mito_hi"), "pct_hb": (None, "pct_hb_hi")}
order = c[["sample", "digestion"]].drop_duplicates().sort_values(["digestion", "sample"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6, "pdf.fonttype": 42})
for susp in ["cell", "nucleus"]:
    samples = order[order["sample"].isin(c.loc[c["suspension"] == susp, "sample"])]
    for met, xl in METRICS.items():
        sub = c[c["suspension"] == susp]
        lo_x, hi_x = np.nanpercentile(sub[met], [0.2, 99.8])
        pad = 0.05 * (hi_x - lo_x)
        xs = np.linspace(lo_x - pad, hi_x + pad, 300)
        n = len(samples)
        fig, axes = plt.subplots(1, 4, figsize=(7.2, 0.24 * n + 1.0), sharey=True)
        for j, L in enumerate(LIN):
            ax = axes[j]
            for i, (s, dg) in enumerate(samples[["sample", "digestion"]].itertuples(index=False)):
                y0 = n - 1 - i
                v = sub.loc[(sub["sample"] == s) & (sub["L0"] == L), met].dropna().values
                if len(v) >= 10 and np.std(v) > 0:
                    k = gaussian_kde(v)(xs)
                    k = k / k.max() * 0.9
                    ax.fill_between(xs, y0, y0 + k, color=PAL[L], alpha=0.85 if len(v) >= MIN_CELLS else 0.35,
                                    lw=0, zorder=n - i)
                    ax.plot(xs, y0 + k, color="white", lw=0.4, zorder=n - i + 0.1)
                elif len(v):
                    ax.plot(v, np.full(len(v), y0 + 0.1), "|", color=PAL[L], ms=3)
                if met in V1KEY and len(v):
                    tr = T[(T["sample"] == s) & (T["L0"] == L)].iloc[0]
                    for key in V1KEY[met]:
                        if key and np.isfinite(tr[key]):
                            ax.plot([tr[key]] * 2, [y0, y0 + 0.95], color="#D55E00", ls="--", lw=0.9,
                                    zorder=500, path_effects=HALO)
                if met in ("log1p_n_counts", "log1p_n_genes", "pct_mito") and len(v):
                    ar = A.loc[s]
                    for key in ([f"{met}_lo", f"{met}_hi"] if met != "pct_mito" else ["pct_mito_hi"]):
                        ax.plot([ar[key]] * 2, [y0, y0 + 0.95], color="#555555", ls=":", lw=0.8, zorder=499)
                if j == 0:
                    ax.text(xs[0], y0 + 0.15, "", fontsize=5)
            ax.set_xlim(xs[0], xs[-1])
            nL = int((sub["L0"] == L).sum())
            fl = int(F.loc[sub.index][sub["L0"] == L]["v1_fail"].sum()) if nL else 0
            ax.set_title(f"{L}\n{nL:,} cells, v1 fails {fl:,}", fontsize=6)
            ax.set_xlabel(xl)
            ax.spines[["top", "right", "left"]].set_visible(False)
        axes[0].set_yticks(np.arange(n) + 0.3)
        axes[0].set_yticklabels([f"{s} ({dg})" for s, dg in samples[["sample", "digestion"]].itertuples(index=False)][::-1],
                                fontsize=5)
        handles = [Line2D([], [], color="#D55E00", ls="--", lw=1, label="QC v1 (MASLD recipe, sample x L0)"),
                   Line2D([], [], color="#555555", ls=":", lw=1, label="authors: 3 MAD per sample")]
        fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.01))
        fig.suptitle(f"Guilliams 2022 Fig 1B samples - {susp}s - {xl} - before removal", fontsize=7)
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        for ext in ("pdf", "png"):
            fig.savefig(f"{FIG}/QC_ridge_{susp}_{met}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print("figure", susp, met, flush=True)
