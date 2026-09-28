#!/usr/bin/env python3
"""
Treg -> B cell spatial communication (COMMOT)
`adata` = Xenium niche object.
"""
import numpy as np
import pandas as pd
import commot as ct
import omnipath as op
from scipy.stats import mannwhitneyu

# ---- LR pairs: single-gene OmniPath pairs in the panel, TGFB1 only -----------
lr = op.interactions.import_intercell_network(
    transmitter_params={"categories": "ligand"},
    receiver_params={"categories": "receptor"})
panel = set(adata.var_names)
m = (~lr["genesymbol_intercell_source"].str.contains("_", na=False)
     & lr["genesymbol_intercell_source"].isin(panel)
     & lr["genesymbol_intercell_target"].isin(panel))
df_ligrec = (lr[m][["genesymbol_intercell_source", "genesymbol_intercell_target"]]
             .drop_duplicates()
             .rename(columns={"genesymbol_intercell_source": "ligand",
                              "genesymbol_intercell_target": "receptor"}))
df_ligrec["pathway"] = "OmniPath"
df_ligrec = df_ligrec[df_ligrec["ligand"] == "TGFB1"].reset_index(drop=True)

# ---- Run COMMOT per follicle, extract Treg -> GC B signal --------------------
adata = adata[adata.obs["compartment"] == "boundary"].copy()
rows = []
for fid, sub in adata.obs.groupby("follicle"):
    a = adata[sub.index].copy()
    if (a.obs["anno"] == "Treg").sum() < 3 or (a.obs["anno"] == "GC B cell").sum() < 3:
        continue
    ct.tl.spatial_communication(a, database_name="LR", df_ligrec=df_ligrec,
                                dis_thr=100, heteromeric=False, pathway_sum=False)
    key = "commot-LR-TGFB1-TGFBR1"
    if key not in a.obsp:
        continue
    treg = np.where(a.obs["anno"].values == "Treg")[0]
    gcb  = np.where(a.obs["anno"].values == "GC B cell")[0]
    recv = np.asarray(a.obsp[key][np.ix_(treg, gcb)].sum(axis=0)).ravel()  # per-GCB
    rows.append({"follicle": fid, "age_group": a.obs["age_group"][0],
                 "mean_recv": recv.mean()})

df = pd.DataFrame(rows)
young = df[df["age_group"] == "Younger"]["mean_recv"]
old   = df[df["age_group"] == "Older"]["mean_recv"]
stat, p = mannwhitneyu(old, young, alternative="greater")
print(f"Treg->GCB TGFB1-TGFBR1: younger={young.mean():.4f}, older={old.mean():.4f}, p={p:.3g}")
df.to_csv("commot_follicle_TGFB1_TGFBR1.csv", index=False)
