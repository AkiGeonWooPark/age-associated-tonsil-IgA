#!/usr/bin/env python3
"""
Treg - GC B cell spatial co-occurrence (Squidpy)

`adata` = Xenium niche object.
"""
import numpy as np
import pandas as pd
import squidpy as sq
from scipy.stats import wilcoxon

adata = adata[adata.obs["compartment"] == "boundary"].copy()

# Co-occurrence label: Treg / RUNX3+ GC B / RUNX3- GC B / other
def make_label(r):
    if r["anno"] == "Treg":
        return "Treg"
    if r["anno"] == "GC B cell":
        return "GCB_RUNX3pos" if r["runx3_pos"] == 1 else "GCB_RUNX3neg"
    return "other"
adata.obs["co_label"] = adata.obs.apply(make_label, axis=1).astype("category")

# ---- Per-follicle Delaunay graph + co-occurrence ----------------------------
rows = []
for fid, sub in adata.obs.groupby("follicle"):
    a = adata[sub.index].copy()
    sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True)
    sq.gr.co_occurrence(a, cluster_key="co_label")
    occ = a.uns["co_label_co_occurrence"]["occ"]          # (clusters, clusters, bins)
    intervals = a.uns["co_label_co_occurrence"]["interval"]
    mids = (intervals[:-1] + intervals[1:]) / 2
    cats = list(a.obs["co_label"].cat.categories)
    treg_i = cats.index("Treg")
    for tgt in ["GCB_RUNX3pos", "GCB_RUNX3neg"]:
        if tgt not in cats:
            continue
        curve = occ[cats.index(tgt), treg_i, :]           # co-occurrence vs Treg
        m = mids < 30
        rows.append({"follicle": fid, "runx3": tgt,
                     "age_group": a.obs["age_group"][0],
                     "AUC": np.trapz(curve[m], mids[m])})

df = pd.DataFrame(rows)

# ---- Paired RUNX3+ vs RUNX3- (0-30 um AUC) ----------------------------------
wide = df.pivot_table(index="follicle", columns="runx3", values="AUC").dropna()
stat, p = wilcoxon(wide["GCB_RUNX3pos"], wide["GCB_RUNX3neg"])
print(f"co-occurrence AUC RUNX3+ vs RUNX3- (n={len(wide)} follicles): p={p:.3g}")
df.to_csv("RUNX3_Treg_cooccurrence_follicle.csv", index=False)
