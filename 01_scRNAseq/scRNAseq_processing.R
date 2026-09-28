# =============================================================================
# Single-cell RNA-seq processing (Seurat v5)
#
# Gene-expression libraries from 48 tonsil and 16 mesenteric lymph node (mLN)
# samples were processed independently with Cell Ranger v7.1.0 against the
# GRCh38-2020-A reference. This script shows the common downstream flow:
#   QC -> normalization -> HVG -> PCA -> Harmony integration -> doublet removal
#   -> UMAP / graph-based clustering -> manual annotation -> lineage re-clustering
#
# =============================================================================

library(Seurat)          # v5
library(harmony)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)

# ---- 1. Quality control ------------------------------------------------------
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 10)

# ---- 2. Normalization, HVG, PCA, Harmony integration -------------------------
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- IntegrateLayers(obj, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = "harmony")

# ---- 3. Doublet removal (scDblFinder, per sample) ---------------------------
sce <- obj
sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
sce <- as.SingleCellExperiment(sce)
sce <- scDblFinder(sce, samples = "donor")
obj$doublet <- sce$scDblFinder.class
obj <- subset(obj, subset = doublet == "singlet")

# ---- 4. Embedding and graph-based clustering ---------------------------------
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)

# ---- 5. Manual annotation of broad lineages ---------------------------------
# Clusters were annotated from canonical markers (e.g. CD3D/NKG7 = T/NK,
# MS4A1/CD40 = B, LMO2/MKI67 = GC B, JCHAIN/XBP1 = plasma, CR2 = FDC).
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
obj <- RenameIdents(obj, ...)          # broad lineage labels
obj$lineage <- Idents(obj)

# ---- 6. Lineage extraction and refined re-clustering ------------------------
# B, T/NK and myeloid compartments were extracted and independently re-clustered
# (same Normalize -> HVG -> PCA -> Harmony -> UMAP -> Leiden flow) to resolve
# fine cell states, then merged back into a single object. This yielded 26
# annotated cell types in tonsil and 21 in mLN (stored in `anno`).

saveRDS(obj, "scRNA_annotated.rds")

# =============================================================================
# Cross-tissue cell-type similarity (tonsil vs mLN)
# Average log-normalized expression was computed per annotated cell type using
# genes shared between the two datasets, and pairwise Spearman correlations were
# calculated between tonsil and mLN cell types.
# =============================================================================
# shared <- intersect(rownames(tonsil), rownames(mln))
# t_mean <- AverageExpression(tonsil, features = shared, group.by = "anno")$RNA
# m_mean <- AverageExpression(mln,    features = shared, group.by = "anno")$RNA
# cor_mat <- cor(t_mean, m_mean, method = "spearman")
