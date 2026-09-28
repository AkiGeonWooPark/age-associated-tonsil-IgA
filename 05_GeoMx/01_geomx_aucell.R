# =============================================================================
# GeoMx DSP pathway activity (AUCell via escape)
#
# GeoMx data were QC-filtered and Q3-normalized in the GeoMx DSP Analysis Suite.
# `mat` = Q3-normalized target-count matrix (genes x segments);
# `meta` = per-segment table with generic columns `segment`, `region`, `age_group`.
# =============================================================================

library(escape)
library(msigdbr)
library(dplyr)
library(tibble)

# ---- 1. Gene sets from MSigDB (C2) ------------------------------------------
msig <- msigdbr(species = "Homo sapiens", collection = "C2")
gs <- list(
  KEGG_IgA = msig %>% filter(gs_name == "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION") %>%
             pull(gene_symbol) %>% unique(),
  M18933   = msig %>% filter(gs_name == "BIOCARTA_TGFB_PATHWAY") %>%
             pull(gene_symbol) %>% unique()
)

# ---- 2. AUCell enrichment per segment ---------------------------------------
scores <- escape.matrix(mat, gene.sets = gs, method = "AUCell",
                        groups = 30, min.size = NULL)
scores <- as.data.frame(scores) %>% rownames_to_column("segment")

# ---- 3. Join metadata and compare age groups by compartment -----------------
scores <- left_join(scores, meta, by = "segment")

# CD20+ (follicular) and CD20- (extrafollicular) analysed separately;
# age-group difference tested with a two-sided Mann-Whitney U test.
cd20pos <- scores %>% filter(region == "CD20pos")
cd20neg <- scores %>% filter(region == "CD20neg")
# wilcox.test(KEGG_IgA ~ age_group, data = cd20pos)   # etc.
write.csv(scores, "geomx_aucell_scores.csv", row.names = FALSE)
