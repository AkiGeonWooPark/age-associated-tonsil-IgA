#!/usr/bin/env python3
"""
RUNX3 transcription-factor activity inference (decoupleR)

"""
import decoupler as dc
from scipy.stats import wilcoxon

# ---- 1. TF activity with CollecTRI + ULM ------------------------------------
collectri = dc.op.collectri(organism="human")
dc.mt.ulm(data=adata, net=collectri, layer="lognorm")
score = dc.pp.get_obsm(adata=adata, key="score_ulm")     # cells x TF activity
adata.obs["runx3_activity"] = score[:, "RUNX3"].X.ravel()

# ---- 2. Paired IgA vs IgG within expanded clones ----------------------------
df = adata.obs[["clone", "anno", "isotype", "runx3_activity"]]
df = df[df["isotype"].isin(["IgA", "IgG"])]
paired = (df.groupby(["clone", "anno", "isotype"])["runx3_activity"]
            .mean().unstack("isotype").dropna(subset=["IgA", "IgG"]))

stat, p = wilcoxon(paired["IgA"], paired["IgG"])
print(f"paired IgA vs IgG RUNX3 activity: n={len(paired)}, p={p:.3g}")
