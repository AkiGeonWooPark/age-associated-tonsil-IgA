# age-associated-tonsil-IgA

Analysis code for **"Age-associated follicular remodeling and selective IgA
maturation in the pediatric tonsil."**

This repository contains distilled, commented scripts that document the main
analysis flow for each data modality. They are intended to make the analytical
steps in the manuscript reproducible in outline; data loading, paths, and
plotting/figure code are simplified or omitted, and metadata column names are
generic placeholders (e.g. `anno`, `age_group`, `donor`). Each script assumes its
input object is already in memory and maps to a section of the Methods.

## Data

All raw and processed data are deposited in the NCBI Gene Expression Omnibus (GEO):

| Modality | Cohort | Accession |
|---|---|---|
| scRNA-seq, scBCR-seq | Tonsil | GSE311663 |
| scTCR-seq | Tonsil | GSE347858 |
| GeoMx DSP | Tonsil | GSE347694, GSE347727 |
| Xenium in situ | Tonsil | GSE347859 |
| scRNA/scBCR/scTCR-seq | Mesenteric lymph node | GSE347985 |

## Repository layout

| Directory | Analysis | Methods section |
|---|---|---|
| `01_scRNAseq/` | Seurat v5 processing (QC, Harmony, scDblFinder, annotation, cross-tissue correlation), pseudobulk DESeq2 DE, and decoupleR RUNX3 activity (CollecTRI + ULM, paired IgA vs IgG) | *Single-Cell RNA Sequencing Analysis*; *…differential expression analysis*; *RUNX3 transcription factor activity inference* |
| `02_BCR_repertoire/` | Immcantation BCR pipeline (clonal assignment, SHM, diversity, isotype, BASELINe) and multivariable regression of clonal IgA proportion (forest plot) | *Single-Cell BCR Repertoire Analysis*; *Multivariable analysis of clonal IgA proportion* |
| `03_TCR_repertoire/` | scirpy TCR pipeline: TCRdist clonotype clusters, clonal expansion, Shannon diversity | *Single-Cell TCR Repertoire Analysis* |
| `04_CyTOF/` | Scanpy CyTOF pipeline (arcsinh, Leiden, annotation) | *Mass Cytometry (CyTOF) Analysis* |
| `05_GeoMx/` | GeoMx AUCell pathway scores (escape) and BulkSignalR TGFB1 ligand–receptor inference | *GeoMx DSP data processing and analysis* |
| `06_Xenium/` | stLearn/Scanpy processing, follicle annotation & signed distance, KNN/k-means niches, Treg→B-cell COMMOT communication, Treg–GC B squidpy co-occurrence | *Xenium … Analysis*; *… follicle annotation*; *Xenium Niche Analysis*; *Treg–B cell spatial communication*; *Treg–B cell co-occurrence analysis* |

Each analysis was applied independently to the tonsil and mesenteric lymph node
cohorts where both were profiled; the scripts show the shared flow (tonsil paths
as example).

## Software

- **R** (v4.3) — Seurat v5, harmony, scDblFinder, DESeq2, the Immcantation
  suite (alakazam, shazam, scoper, dowser), escape + msigdbr (AUCell), and
  BulkSignalR.
- **Python** (v3.9) — scanpy, scirpy, decoupleR, stLearn, squidpy, COMMOT,
  omnipath, scikit-learn, OpenCV.

Statistical tests, multiple-testing correction, and figure generation follow the
*Statistical Analysis* section of the manuscript.
