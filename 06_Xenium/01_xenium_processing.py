#!/usr/bin/env python3
"""
Xenium in situ spatial transcriptomics processing (stLearn / Scanpy)

`adata` holds the Xenium cell-by-gene matrix (430-gene panel, 15 tonsil cores)
with spatial coordinates. 
"""
import numpy as np
import scanpy as sc
import stlearn as st

# ---- 1. QC ------------------------------------------------------------------
st.pp.filter_genes(adata, min_counts=10)   # genes in >=10 cells
st.pp.filter_cells(adata, min_counts=10)   # cells with >=10 transcripts

# ---- 2. Normalization + Freeman-Tukey variance stabilization ----------------
st.pp.normalize_total(adata)
adata.X = np.sqrt(adata.X.toarray()) + np.sqrt(adata.X.toarray() + 1)

# ---- 3. Embedding and clustering --------------------------------------------
st.em.run_pca(adata, n_comps=50, random_state=0)
st.pp.neighbors(adata, n_neighbors=25, use_rep="X_pca", random_state=0)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden")
sc.tl.umap(adata)

# ---- 4. Manual lineage annotation -------------------------------------------
# Major lineages annotated from canonical markers: B cells, T/NK cells, plasma
# cells, epithelial cells, FRCs, and myeloid cells.
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
adata.obs["lineage"] = adata.obs["leiden"].replace({...})

# ---- 5. Per-lineage re-clustering for detailed annotation -------------------
# Each major lineage was extracted and independently re-clustered (same flow) to
# yield the fine cell-type labels used downstream (`anno`).
adata_B = adata[adata.obs["lineage"] == "B_cell"].copy()
# ... re-cluster + annotate, then write back into the full object.

adata.write("adata_xenium_annotated.h5ad")
