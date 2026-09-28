#!/usr/bin/env python3
"""
Xenium follicle annotation and RUNX3+ B-cell spatial distribution

`adata` = annotated Xenium object. 
"""
import numpy as np
import cv2
from sklearn.mixture import GaussianMixture
import stlearn as st

# ---- 1. Build a spatial grid per core ---------------------------------------
grid = st.tl.cci.grid(adata, n_row=250, n_col=250, use_label="anno")

# ---- 2. Extract follicle (GC) contours from the grid ------------------------
GC_LABELS = ["Mantle_B", "GC_B", "Tfh", "FDC", "Plasmablast"]

def extract_follicle_contours(grid, gc_labels, min_area=50):
    """Binary occupancy mask of follicle-associated grid bins -> smoothed
    contours via OpenCV. Contours were manually reviewed for coherence."""
    mask = grid.obs["anno"].isin(gc_labels).astype(np.uint8)
    img = mask.values.reshape(grid.uns["grid_shape"])       # 2D occupancy
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    img = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel)
    img = cv2.morphologyEx(img, cv2.MORPH_CLOSE, kernel)
    cts, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    return [c for c in cts if cv2.contourArea(c) >= min_area]

contours = extract_follicle_contours(grid, GC_LABELS)

# ---- 3. Assign nearest follicle + signed distance ---------------------------
def signed_distance_to_boundary(xy, contour):
    """+ outside / - inside the contour (grid units, later scaled to um)."""
    d = np.array([cv2.pointPolygonTest(contour, (float(x), float(y)), True)
                  for x, y in xy])
    return -d   # pointPolygonTest is + inside; flip so inside is negative
# ... loop over cores/contours, take nearest follicle per cell, store
#     adata.obs["follicle"] and adata.obs["signed_dist_um"].

# ---- 4. Spatial compartments -------------------------------------------------
def compartment(d_um):
    if d_um < -30:  return "core"
    if d_um <  40:  return "boundary"
    if d_um <  80:  return "outer_boundary"
    return "far_outside"
adata.obs["compartment"] = adata.obs["signed_dist_um"].map(compartment)

# ---- 5. RUNX3+ B cells (2-component GMM threshold) ---------------------------
bcells = adata[adata.obs["lineage"] == "B_cell"].copy()
x = np.asarray(bcells[:, "RUNX3"].X.todense()).ravel().reshape(-1, 1)
gm = GaussianMixture(n_components=2, random_state=0).fit(x)
xs = np.linspace(*np.sort(gm.means_.ravel()), 1000).reshape(-1, 1)
thr = xs[np.argmin(np.abs(np.diff(gm.predict_proba(xs), axis=1)))][0]
bcells.obs["runx3_pos"] = (x.ravel() > thr).astype(int)
