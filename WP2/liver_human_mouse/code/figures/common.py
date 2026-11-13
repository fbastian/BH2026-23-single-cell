"""Shared figure settings (Nature spec: 183 mm wide, <= 170 mm tall, Arial 5-7 pt, panel letters 8 pt bold, TrueType text,
points rasterised at 450 dpi, Wong palette) + lineage maps for mouse (annotation_aor) and HLiCA (author_cell_type)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb
plt.rcParams.update({"font.family": "Arial", "font.size": 6, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
                     "axes.titlesize": 6, "axes.labelsize": 6, "xtick.labelsize": 5, "ytick.labelsize": 5, "legend.fontsize": 5,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 450})
MM = 1 / 25.4
D = "/mnt/masld_vault/cross_species/figures"
WONG = {"Hepatocyte": "#E69F00", "Cholangiocyte": "#000000", "Endothelial": "#0072B2", "Mesenchymal": "#009E73",
        "Myeloid": "#D55E00", "DC": "#CC79A7", "B/plasma": "#F0E442", "T/NK/ILC": "#56B4E9", "Other": "#BBBBBB"}
LIN_ORDER = list(WONG)
MOUSE_LIN = {"Hepatocyte": "Hepatocyte", "Cholangiocyte": "Cholangiocyte", "Endothelial": "Endothelial", "Fibroblast": "Mesenchymal",
             "Mesothelial": "Mesenchymal", "KC": "Myeloid", "Monocyte": "Myeloid", "Neutrophil": "Myeloid", "Basophil": "Myeloid",
             "cDC1": "DC", "cDC2": "DC", "Mig. cDCs": "DC", "Transitional DC": "DC", "pDCs": "DC", "B cell": "B/plasma",
             "T cell": "T/NK/ILC", "NK cell": "T/NK/ILC", "ILC1": "T/NK/ILC", "unresolved": "Other"}
MOUSE_ORDER = list(MOUSE_LIN)
HUMAN_LIN = {**{k: "Hepatocyte" for k in ["Periportal Hepatocyte", "Pericentral Hepatocyte", "SERPINE1+ Hepatocyte", "UGT+ Hepatocyte",
                                          "Mito+ Hepatocyte", "Ribosomal+ Hepatocyte"]},
             **{k: "Cholangiocyte" for k in ["LAMC2+ Cholangiocyte", "ApoLipo", "Keratin", "CXCL8+ Keratin", "Mucus Secreting"]},
             **{k: "Endothelial" for k in ["Periportal LSEC", "LSEC", "Central Venous LSEC", "Portal Vein", "Hepatic Artery",
                                           "Vascular Endothelial Cell", "Lymphatic"]},
             **{k: "Mesenchymal" for k in ["Portal Fibroblast", "Hepatic Stellate Cell", "CUX2+ Hepatic Stellate Cell",
                                           "Vascular Smooth Muscle Cell", "NRXN1+ Stromal Cell"]},
             **{k: "Myeloid" for k in ["Kupffer Cells", "TREM2+ Macrophages", "Classical Monocytes", "Non-classical Monocytes",
                                       "Activated Monocytes", "BAG3+ Monocytes", "MAMLD1+ Trans Monocytes", "Neutrophils"]},
             **{k: "DC" for k in ["Type 1 cDCs", "Type 2 cDCs", "DCmac", "migDCs", "pDC"]},
             **{k: "B/plasma" for k in ["B Cell", "IgA B cells", "IgG B cells", "Plasma B Cell"]},
             **{k: "T/NK/ILC" for k in ["Helper T Cell", "Regulatory T Cell", "CD8 T Cell", "Exhausted CD8 T Cell", "MAIT T Cell",
                                        "Bright NK Cell", "Dim NK Cell"]},
             "Cycling": "Other"}
HUMAN_ORDER = list(HUMAN_LIN)


def shades(order, lin):
    """Per-type colours: the lineage's Wong colour, shaded light->dark within the lineage."""
    out = {}
    for L in LIN_ORDER:
        ts = [t for t in order if lin[t] == L]
        base = np.array(to_rgb(WONG[L]))
        for i, t in enumerate(ts):
            f = 0 if len(ts) == 1 else (i / (len(ts) - 1) - 0.5) * 0.7  # -0.35..0.35
            c = base + (1 - base) * (-f) if f < 0 else base * (1 - f)
            if L == "Cholangiocyte":
                c = np.array([0.15 + 0.12 * i] * 3)
            out[t] = tuple(np.clip(c, 0, 1))
    return out


def ax_mm(fig, x, y, w, h, W=183, H=170):
    """axes from mm, origin top-left."""
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def letter(fig, x, y, s, W=183, H=170):
    fig.text(x / W, 1 - y / H, s, fontsize=8, fontweight="bold", va="top", ha="left")


def umap_ax(ax, U, colors, s=0.15, order=None):
    o = np.random.default_rng(0).permutation(len(U)) if order is None else order
    ax.scatter(U[o, 0], U[o, 1], c=np.asarray(colors, dtype=object)[o].tolist(), s=s, lw=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    # arrows a fixed 7 mm long whatever the panel size (user 2026-11-12: arrows too small in small panels)
    pos = ax.get_position(); fw, fh = ax.figure.get_size_inches()
    fx, fy = 7 / (pos.width * fw * 25.4), 7 / (pos.height * fh * 25.4)
    x0, y0 = U[:, 0].min(), U[:, 1].min(); dx, dy = np.ptp(U[:, 0]) * fx, np.ptp(U[:, 1]) * fy
    ax.annotate("", (x0 + dx, y0), (x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.5, color="k", mutation_scale=4))
    ax.annotate("", (x0, y0 + dy), (x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.5, color="k", mutation_scale=4))
    ax.text(x0, y0 - dy * 0.15, "UMAP1", fontsize=5, va="top"); ax.text(x0 - dx * 0.15, y0, "UMAP2", fontsize=5, rotation=90, ha="right")


# HLiCA authors' palette (github quon-titative-biology/HLiCA src/Annotations/Liver_Cell_Atlas/Analysis_03/colour_palette.R,
# myColors_celltype; names assigned positionally as in the R code, so LSEC / Vascular Endothelial Cell follow R, not the comments)
_HC = ["#994676", "#ce1256", "#d9667f", "#ce1256", "#b80783", "#b01629", "#d9667f", "#914e84", "#b01629", "#d9667f", "#d9667f", "#ce1256",
       "#FFA500", "#FF8C00", "#FFB347", "#FF7F50", "#FFAA33",
       "#8B4513", "#A0522D", "#CD853F", "#B5651D", "#DEB887", "#8B3E2F", "#A52A2A",
       "#F0E68C", "#FFD700", "#FFFACD", "#F5DEB3", "#f09d02", "#f09d02",
       "#663399", "#800080", "#9370DB", "#9851b0",
       "#006400", "#339500", "#4CAF00", "#66C700", "#b0db93", "#5f9c4e", "#4c9638",
       "#6495ED", "#1d91c0", "#74a9cf", "#9ecae1", "#0570b0", "#08519c", "#eb150c", "#000066"]
_HN = ["Kupffer Cells", "Classical Monocytes", "Type 2 cDCs", "Non-classical Monocytes", "TREM2+ Macrophages", "MAMLD1+ Trans Monocytes", "DCmac",
       "Neutrophils", "Activated Monocytes", "migDCs", "Type 1 cDCs", "BAG3+ Monocytes",
       "ApoLipo", "Keratin", "LAMC2+ Cholangiocyte", "CXCL8+ Keratin", "Mucus Secreting",
       "Periportal LSEC", "Portal Vein", "Central Venous LSEC", "Lymphatic", "Hepatic Artery", "Vascular Endothelial Cell", "LSEC",
       "Hepatic Stellate Cell", "Vascular Smooth Muscle Cell", "Portal Fibroblast", "CUX2+ Hepatic Stellate Cell", "NRXN1+ Hepatic Stellate Cell",
       "NRXN1+ Stromal Cell", "Plasma B Cell", "B Cell", "IgA B cells", "IgG B cells",
       "Bright NK Cell", "Dim NK Cell", "Helper T Cell", "CD8 T Cell", "MAIT T Cell", "Exhausted CD8 T Cell", "Regulatory T Cell",
       "Ribosomal+ Hepatocyte", "UGT+ Hepatocyte", "Mito+ Hepatocyte", "Pericentral Hepatocyte", "Periportal Hepatocyte", "SERPINE1+ Hepatocyte",
       "pDC", "Cycling"]
HLICA_COLORS = dict(zip(_HN, _HC))
# HLiCA lineages (myColors_celltype_alpha)
HLICA_LIN_COLORS = {"Hepatocyte": "#6495ED", "Cholangiocyte": "#FFA500", "Endothelia": "#CD853F", "Mesenchyme": "#FFD700",
                    "Myeloid Cell": "#ce1256", "Lymphocyte": "#339500", "pDC": "#eb150c", "Cycling": "#000066"}
HLICA_LIN = {}
for t in _HN:
    L = HUMAN_LIN.get(t)
    HLICA_LIN[t] = {"Hepatocyte": "Hepatocyte", "Cholangiocyte": "Cholangiocyte", "Endothelial": "Endothelia", "Mesenchymal": "Mesenchyme",
                    "Myeloid": "Myeloid Cell", "DC": "Myeloid Cell", "B/plasma": "Lymphocyte", "T/NK/ILC": "Lymphocyte"}.get(L, t)
HLICA_LIN["pDC"] = "pDC"; HLICA_LIN["Cycling"] = "Cycling"
HLICA_ORDER = [t for L in HLICA_LIN_COLORS for t in HUMAN_ORDER if HLICA_LIN.get(t) == L]
