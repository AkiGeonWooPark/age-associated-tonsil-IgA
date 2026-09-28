# =============================================================================
# Pseudobulk differential expression (DESeq2)
# Two comparisons are shown:
#   (1) age-associated : younger vs older, within each B-cell subset
#   (2) isotype-associated : IgA+ vs IgG+ cells within each B-cell subset
# Genes with adjusted p < 0.05 and |log2FC| > 0.5 were considered significant.
# =============================================================================

library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)

# ---- Aggregate to donor x cell-type x age-group pseudobulk -------------------
pseudo <- AggregateExpression(
  obj, assays = "RNA", return.seurat = TRUE,
  group.by = c("donor", "anno", "age_group")
)

# ---- Per-subset DESeq2: younger vs older ------------------------------------
# For each subset, genes were pre-filtered to those detected in >=10% of cells in
# either group, then tested with DESeq2 on the pseudobulk object.
anno_levels <- unique(na.omit(obj$anno))

deg_by_subset <- lapply(anno_levels, function(a) {
  sc_sub <- subset(obj, subset = anno == a)
  if (!all(c("Younger", "Older") %in% sc_sub$age_group)) return(NULL)

  expr  <- GetAssayData(sc_sub, assay = "RNA", layer = "counts")
  g     <- sc_sub$age_group
  pct_y <- Matrix::rowMeans(expr[, g == "Younger", drop = FALSE] > 0)
  pct_o <- Matrix::rowMeans(expr[, g == "Older",   drop = FALSE] > 0)
  genes_keep <- names(which(pct_y >= 0.1 | pct_o >= 0.1))
  if (length(genes_keep) < 2) return(NULL)

  obj_sub <- subset(pseudo, subset = anno == a)
  Idents(obj_sub) <- "age_group"
  if (!all(c("Younger", "Older") %in% Idents(obj_sub))) return(NULL)
  obj_sub <- subset(obj_sub, features = intersect(genes_keep, rownames(obj_sub)))

  deg <- FindMarkers(obj_sub, ident.1 = "Older", ident.2 = "Younger",
                     test.use = "DESeq2")
  deg$gene <- rownames(deg); deg$anno <- a
  deg
})
deg_age <- bind_rows(deg_by_subset)
write.xlsx(deg_age, "DEG_Older_vs_Younger_bySubset.xlsx")

# ---- Isotype comparison (IgA+ vs IgG+) within each subset -------------------
# Same pseudobulk/DESeq2 approach, grouping cells by `isotype` instead of
# `age_group`; tested per subset (see 02_BCR_repertoire for the isotype call).
