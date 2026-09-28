# =============================================================================
# GeoMx DSP ligand-receptor inference (BulkSignalR)
# `counts` = Q3-normalized target-count matrix (col 1 = TargetName, then segments);
# `meta`   = per-segment table with generic columns `segment`, `region`, `age_group`.
# =============================================================================

library(BulkSignalR)
library(dplyr)
library(openxlsx)

# ---- CD20+ segments split by age group --------------------------------------
cd20      <- meta %>% filter(region == "CD20pos")
seg_young <- cd20 %>% filter(age_group == "Younger") %>% pull(segment)
seg_old   <- cd20 %>% filter(age_group == "Older")   %>% pull(segment)

# ---- 1. Build data model + LR inference per age group ------------------------
run_bsr <- function(count_mat) {
  dm  <- BSRDataModel(counts = count_mat, normalize = FALSE, method = "UQ",
                      log.transformed = FALSE)
  dm  <- learnParameters(dm, quick = TRUE)
  inf <- BSRInference(dm, min.cor = 0.3, reference = "GOBP")
  list(dm = dm, inf = inf, lr = LRinter(inf))
}
bsr_young <- run_bsr(counts[, c("TargetName", seg_young)])
bsr_old   <- run_bsr(counts[, c("TargetName", seg_old)])

# ---- 2. TGFB1-mediated interactions, younger vs older -----------------------
# LRinter() gives L, R, pathway, LR.corr, and BH-adjusted qval per interaction.
tgfb1_young <- bsr_young$lr %>% filter(L == "TGFB1")
tgfb1_old   <- bsr_old$lr   %>% filter(L == "TGFB1")

comp <- full_join(
  tgfb1_young %>% select(L, R, pw.id, pw.name, LRcorr_young = LR.corr, qval_young = qval),
  tgfb1_old   %>% select(L, R, pw.id, LRcorr_old = LR.corr, qval_old = qval),
  by = c("L", "R", "pw.id")
) %>% mutate(delta = LRcorr_old - LRcorr_young) %>% arrange(desc(delta))

# Interactions significant in the older group (q < 0.05) and enriched vs younger
enriched_old <- comp %>% filter(qval_old < 0.05, is.na(qval_young) | delta > 0)
write.xlsx(comp, "BulkSignalR_TGFB1_Younger_vs_Older.xlsx")
