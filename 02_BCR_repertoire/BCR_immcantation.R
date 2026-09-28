# =============================================================================
# Single-cell BCR repertoire analysis (Immcantation)
#
# `data_all` is the combined AIRR rearrangement table (IgBLAST + Change-O output).
# =============================================================================

library(alakazam)
library(shazam)
library(scoper)
library(dowser)
library(tidyverse)

# ---- QC: keep productive, drop multi-heavy and light-only cells --------------
data_all <- data_all %>% filter(productive)
multi_heavy <- names(which(table(filter(data_all, locus == "IGH")$cell_id) > 1))
data_all <- filter(data_all, !cell_id %in% multi_heavy)
heavy <- data_all %>% filter(locus == "IGH") %>% pull(cell_id)
light <- data_all %>% filter(locus %in% c("IGK", "IGL")) %>% pull(cell_id)
data_all <- filter(data_all, !cell_id %in% setdiff(light, heavy))

# ---- 1. Link to matched transcriptomic annotation (by barcode) --------------
data_all <- left_join(data_all, meta, by = "cell_id") %>% filter(!is.na(anno))

# ---- 2. Per-donor clonal threshold (nearest-neighbour distance) -------------
samples <- unique(data_all$donor)
nearest <- lapply(samples, function(s)
  distToNearest(filter(data_all, locus == "IGH", donor == s),
                cellIdColumn = "cell_id", model = "ham"))
names(nearest) <- samples

thresholds <- lapply(nearest, function(x)
  findThreshold(x$dist_nearest, method = "gmm", model = "gamma-norm",
                cutoff = "user", spc = 0.995)@threshold)

# ---- 3. Clonal assignment (per donor) ---------------------------------------
clones <- lapply(samples, function(s)
  hierarchicalClones(filter(data_all, donor == s),
                     cell_id = "cell_id", threshold = thresholds[[s]],
                     only_heavy = TRUE, split_light = TRUE,
                     summarize_clones = FALSE, fields = "donor"))
names(clones) <- samples
# Clones with >=2 cells were considered expanded.

# ---- 4. Germline reconstruction and SHM frequency ---------------------------
references <- readIMGT(dir = "/path/to/human/vdj/")
clones <- lapply(clones, function(x) createGermlines(x, references, nproc = 1))

bcr_heavy <- bind_rows(clones) %>% filter(locus == "IGH")
bcr_mut <- observedMutations(bcr_heavy,
                             sequenceColumn = "sequence_alignment",
                             germlineColumn = "germline_alignment_d_mask",
                             regionDefinition = IMGT_V, frequency = TRUE,
                             combine = TRUE, nproc = 4)   # mu_freq = SHM

# Isotype grouping: IGHA/IGHG = class-switched, IGHM/IGHD = unswitched
bcr_mut <- bcr_mut %>% mutate(isotype = case_when(
  grepl("^IGHA", c_call) ~ "IgA",
  grepl("^IGHG", c_call) ~ "IgG",
  grepl("^IGHM|^IGHD", c_call) ~ "IgM_IgD",
  TRUE ~ "Other"))

# ---- 5. Clonal diversity (Shannon entropy per donor x subset) ---------------
shannon <- function(x) { p <- table(x) / length(x); -sum(p * log(p)) }
diversity_df <- bcr_mut %>%
  group_by(donor, age, anno) %>%
  summarise(shannon_H = shannon(clone_id), n = n(), .groups = "drop")
# Age association tested by Spearman correlation within each subset (BH-adjusted).

# ---- 6. Antigen-driven selection (BASELINe) ---------------------------------
# Clone-specific germlines are reconstructed, clones collapsed, and selection
# strength (Sigma) estimated over the V region, grouped by subset / isotype.
clones_g <- lapply(clones, function(x)
  createGermlines(x, references, v_call = "v_call_genotyped", nproc = 1))
comb <- do.call(rbind, clones_g)
comb <- observedMutations(comb, sequenceColumn = "sequence_alignment",
                          germlineColumn = "germline_alignment_d_mask",
                          regionDefinition = NULL, frequency = FALSE,
                          combine = TRUE, nproc = 4)
z <- collapseClones(comb, cloneColumn = "clone_id")
b <- calcBaseline(z, regionDefinition = IMGT_V)
g_subset  <- groupBaseline(b, groupBy = "anno")      # by cell subset
g_isotype <- groupBaseline(b, groupBy = "isotype")   # by isotype
# testBaseline() / plotBaselineDensity() used to compare selection between groups.
