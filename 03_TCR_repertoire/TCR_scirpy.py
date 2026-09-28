#!/usr/bin/env python3
"""
Single-cell TCR repertoire analysis (scirpy)

`tcr_all` is the concatenated AIRR AnnData (Cell Ranger VDJ).
"""
import numpy as np
import anndata as ad
import mudata as mu
import scipy.sparse as sp
import scirpy as ir

# ---- 1. Build MuData (pair AIRR with matched GEX annotation) ----------------

gex = ad.AnnData(X=sp.csr_matrix((len(tcr_all), 1)), obs=tcr_all.obs[[]].copy())
mdata = mu.MuData({"gex": gex, "airr": tcr_all})

# ---- 2. Chain indexing + QC -------------------------------------------------
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)
mdata = mdata[mdata.obs["airr:chain_pairing"] != "multichain"].copy()
mdata = mdata[~np.isin(mdata.obs["airr:chain_pairing"],
                       ["orphan VDJ", "orphan VJ"])].copy()

# ---- 3. TCRdist-based clonotype clusters ------------------------------------
# Distances from amino-acid CDR3 similarity across paired TRA/TRB chains.
ir.pp.ir_dist(mdata, metric="tcrdist", sequence="aa", cutoff=15)
ir.tl.define_clonotype_clusters(
    mdata, sequence="aa", metric="tcrdist",
    receptor_arms="all", dual_ir="any",
)  # -> obs["airr:cc_aa_tcrdist"]

# ---- 4. Clonal expansion and diversity --------------------------------------
# Clonal expansion = fraction of cells in expanded clonotype clusters.
ir.tl.clonal_expansion(mdata, target_col="cc_aa_tcrdist")
# Repertoire diversity summarized as the Shannon index per T-cell subset.
ir.tl.alpha_diversity(mdata, groupby="gex:anno",
                      target_col="cc_aa_tcrdist", metric="shannon")

mdata.obs[["airr:cc_aa_tcrdist"]].to_csv("TCR_clonotypes.csv")
