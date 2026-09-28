#!/usr/bin/env python3
"""
Xenium cellular-niche analysis and Treg-enriched niche
`adata` = annotated Xenium object.
"""
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans

CELLTYPE_COL, CORE_COL = "anno", "core"
K_NEIGH, N_NICHE = 20, 10

# ---- 1. Per-core KNN neighbourhood composition ------------------------------
ct = adata.obs[CELLTYPE_COL].astype("category")
celltypes = list(ct.cat.categories)
onehot = pd.get_dummies(ct).reindex(columns=celltypes, fill_value=0).values.astype(np.float32)

coords = adata.obsm["spatial"]
cores = adata.obs[CORE_COL].astype(str).values
comp = np.zeros((adata.n_obs, len(celltypes)), dtype=np.float32)
for c in np.unique(cores):
    idx = np.where(cores == c)[0]
    k = min(K_NEIGH + 1, len(idx))                    # +1 to include self
    _, knn = NearestNeighbors(n_neighbors=k).fit(coords[idx]).kneighbors(coords[idx])
    comp[idx] = onehot[idx][knn].sum(axis=1)
comp_frac = comp / np.clip(comp.sum(1, keepdims=True), 1, None)

# ---- 2. k-means niches (N0..N9) ---------------------------------------------
niche = KMeans(n_clusters=N_NICHE, random_state=0, n_init=10).fit_predict(comp_frac)
adata.obs["niche"] = pd.Categorical([f"N{int(x)}" for x in niche])


# ---- 3. Niche ligand-receptor analysis (stLearn, summarized) ----------------
# st.tl.cci.run on a 250x250 grid per core gives LR activity (significant spots).
# Per-core Fisher's exact test compares the Treg niche vs other niches; evidence
# is combined across cores (Fisher's method), with geometric-mean OR and BH FDR
# < 0.05. Age effect: pool significant/non-significant spots in N10 by age group
# and test each Treg-niche LR pair with a two-sided Fisher's exact test (BH).
