"""Table 1 data (2026-11-12): every SAMap match mouse type -> HLiCA type with mapping score >= 0.4 (the 40 edges of Fig 3e),
with the number of gene pairs GenePairFinder returned (align_thr 0.1) and the top 6 (human / mouse; one name when both agree).
Output <D>/tables/table1_samap_matches.tsv (rendered by tab1_kable.R)."""
import sys; sys.path.insert(0, "/home/ubuntu/liver_cross_species/figures/code")
import pandas as pd
from common import D, MOUSE_ORDER, MOUSE_LIN, HUMAN_ORDER, LIN_ORDER
O = "/mnt/masld_vault/cross_species/samap_run/out"
M = pd.read_csv(f"{O}/mapping_table.tsv", sep="\t", index_col=0)
gp = pd.read_csv(f"{O}/gene_pairs_thr0.1.tsv", sep="\t")
rows = []
for m in [t for t in MOUSE_ORDER if "mm_" + t in M.index]:
    for h in [t for t in HUMAN_ORDER if "hs_" + t in M.columns]:
        s = M.loc["mm_" + m, "hs_" + h]
        if s < 0.4: continue
        c = f"hs_{h};mm_{m}"
        pr = [x.split(";") for x in gp[c].dropna().astype(str)] if c in gp.columns else []
        top = ", ".join(a[3:] if a[3:].upper() == b[3:].upper() else f"{a[3:]}/{b[3:]}" for a, b in pr[:6])
        rows.append({"Lineage": MOUSE_LIN[m], "Mouse type": m, "HLiCA type": h, "SAMap score": round(float(s), 2), "Gene pairs": len(pr), "Top gene pairs (human/mouse)": top})
t = pd.DataFrame(rows)
t["lo"] = t.Lineage.map({k: i for i, k in enumerate(LIN_ORDER)}); t["mo"] = t["Mouse type"].map({k: i for i, k in enumerate(MOUSE_ORDER)})
t = t.sort_values(["lo", "mo", "SAMap score"], ascending=[True, True, False]).drop(columns=["lo", "mo"])
t.to_csv(f"{D}/tables/table1_samap_matches.tsv", sep="\t", index=False)
print(t.to_string(index=False)[:3000]); print(len(t))
