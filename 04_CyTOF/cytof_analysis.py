#!/usr/bin/env python3
"""
Mass cytometry (CyTOF) analysis (Scanpy)

`adata` holds live CD45+ cells (marker x cell), gated in Cytobank and already
downsampled to 10,000 cells per sample. Two panels (B: 42 markers, T: 40 markers)
were analyzed the same way; the B-cell panel is shown. 
"""
import numpy as np
import scanpy as sc

# ---- 1. arcsinh transform (cofactor = 5) ------------------------------------
adata.X = np.arcsinh(adata.X / 5.0)

# ---- 2. Embedding and clustering --------------------------------------------
sc.pp.neighbors(adata, n_neighbors=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden")

# ---- 3. Annotate major populations; drop non-immune / granulocytes / doublets
# Clusters annotated from canonical marker expression -> obs["lineage"].
adata_immune = adata[~adata.obs["lineage"].isin(
    ["Non-immune", "Granulocyte", "Doublet"])].copy()

# ---- 4. Extract B-cell compartment and re-cluster ---------------------------
adata_B = adata_immune[adata_immune.obs["lineage"] == "B_cell"].copy()
sc.pp.scale(adata_B)
sc.pp.pca(adata_B)
sc.pp.neighbors(adata_B, n_neighbors=20)
sc.tl.umap(adata_B)
sc.tl.leiden(adata_B, resolution=1.0, key_added="b_anno")
# B subsets annotated from IgD/IgM/IgA/CD38/CD27 -> naive / GC / MBC / plasma.
# IgA+ fraction per subset and Treg frequency compared between age groups
# (two-sided Mann-Whitney U test); see figure scripts.
adata_B.write("CyTOF_Bcells.h5ad")
