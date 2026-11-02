# README

> **Project:** Radiomic Intratumoral Heterogeneity (RITH) Phenotype Analysis in Clear Cell Renal Cell Carcinoma (ccRCC)
> **Language:** R (≥ 4.2)
> **Author / Maintainer:** Qingyuan Zheng
> **Last updated:** 2026.11.2

---

## 1. Project Overview

This repository contains an end-to-end R workflow for identifying and characterising a **Radiomic Intratumoral Heterogeneity (RITH)** phenotype in ccRCC from CT-derived PyRadiomics features.

The workflow:

1. **Cleans** the raw radiomics feature table and removes samples with inconsistent heterogeneity.
2. **Filters** features by median absolute deviation (MAD > 1).
3. **Determines feature direction** (HHR / LHR) based on the correlation of each feature with a pseudo-heterogeneity-tendency score.
4. **Clusters** samples with **Similarity Network Fusion (SNF) + spectral clustering (K = 2)** into **High RITH** and **Low RITH** subgroups.
5. **Computes a Heterogeneity Index (HI)** as `Σ z(HHR) − Σ z(LHR)`.
6. **Validates** the phenotype by PCA, boxplot / violin / density plots, and statistical testing.
7. **Associates** HI with clinical variables (pT stage, ISUP grade), survival (PFI), morphological features, and bulk transcriptomics (DEG, GO / KEGG / GSEA).

The final outputs are: the cleaned feature table with a `RITH_Phenotype` label, the per-sample `Heterogeneity_Index`, and a series of figures / tables describing how HI relates to prognosis and biology.

---

## 2. Directory Structure

```
ccRCC/
├── 0.Data_Cleaning.R                    # Clean features, define direction, compute HI, drop HI<0 in High-RITH
├── 1.Phenotype_Clustering_k2.R          # k=2 SNF phenotype clustering (fixed seed 1234)
├── 2.PCA_Dimensionality_Reduction.R     # PCA on MAD-selected features
├── 3.HI_Boxplot.R                       # Violin / box / density / histogram of HI
├── 4.Feature_HI_Correlation.R           # Feature-wise Spearman correlation with HI
├── 5.TMEimmune26.geneSurvival.R         # Kaplan–Meier PFI survival analysis
├── 6.HI_Clinical_Variable_Correlation.R # HI vs pT stage / ISUP grade
├── 7.Pathological_MITH.R                # Morphological features (Morph_*) vs HI
├── 8.Biological_Analysis.R              # limma DEG, volcano, GO / KEGG / GSEA
└── README.md
```

The reference working directory used by the scripts is `/mnt/data4/表型研究2026/ccRCC/RHWU/`. All sub-scripts (`2`–`8`) are expected to live in subfolders such as `2.pca降维`, `3.HI箱线图`, etc.; their `setwd()` calls reflect that layout and **must be adjusted** if you run them elsewhere.

---

## 3. Dependencies

### 3.1 R packages

```r
# Data manipulation
install.packages(c("readr", "dplyr", "tidyr", "stringr"))

# Clustering / dimensionality reduction
install.packages(c("SNFtool", "Spectrum", "cluster", "NbClust", "factoextra"))

# Visualisation
install.packages(c("ggplot2", "pheatmap", "ggrepel", "gridExtra"))

# Survival
install.packages(c("survival", "survminer"))

# Statistics / clinical plots
install.packages(c("ggpubr", "rstatix"))

# Transcriptomics / enrichment
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "org.Hs.eg.db",
                       "clusterProfiler", "enrichplot",
                       "ComplexHeatmap", "circlize", "RColorBrewer"))
```

### 3.2 External data files

| File | Used by | Required columns |
|---|---|---|
| `RITH_Data_Cleaned.csv` | `0`–`4` | `Image`, PyRadiomics features (`original_*`, `log-sigma-*-mm-3D_*`, `wavelet-*`) |
| `survival.txt` | `5` | `RITH_Phenotype`, `PFI.time`, `PFI` |
| `HI.csv` | `6` | `Heterogeneity_Index`, `pT_stage`, `ISUP_grade` |
| `Feats_T.csv` | `7` | `id`, `RITH_Phenotype`, `Morph_*` |
| `RHWU_ccRCC_bulk_RNA_Seq_gene_counts.csv` | `8` | Ensembl IDs as rownames; sample names as columns |
| `c2.cp.kegg.Hs.symbols.gmt` | `8` | Standard MSigDB KEGG GMT (SYMBOL) |

> The RNA-seq sample names in `8.Biological_Analysis.R` (`"王风云11010926"`, `"任志海11029686"`, …) are **data identifiers** and must be kept as-is, unless you rename the input columns.

---

## 4. Workflow and Execution Order

The scripts are designed to run **sequentially**:

```
0.Data_Cleaning.R
        │   → RITH_Data_Cleaned.csv (samples with HI<0 in High-RITH removed)
        ▼
1.Phenotype_Clustering_k2.R
        │   → RITH_Phenotype_and_HI_Final.csv
        │   → RITH_Direction_Separated_Heatmap.pdf
        ▼
2.PCA_Dimensionality_Reduction.R ──► Standardized_PC_Scores.csv / Final_RITH_Results_with_PC.csv
3.HI_Boxplot.R                   ──► HI_Violin_Plot / HI_Box_Plot / HI_Density_Plot / HI_Histogram
4.Feature_HI_Correlation.R       ──► All_Features_HI_Correlation.csv / HI_Scatter_Method1_FullNames.pdf
5.TMEimmune26.geneSurvival.R     ──► sur.RITH_Phenotype.tiff  (KM plot, PFI)
6.HI_Clinical_Variable_Correlation.R ──► HI_correlation_plots.pdf / HI_correlation_results.csv
7.Pathological_MITH.R            ──► Morphological_Analysis_Results/ (01–08)
8.Biological_Analysis.R          ──► DEG_results_limma.csv / volcano*.pdf / GO*.pdf / KEGG*.pdf / GSEA*.pdf
```

Each script is standalone: it re-loads `RITH_Data_Cleaned.csv` (or its own input file), re-derives `RITH_Phenotype` and `Heterogeneity_Index` from scratch, and writes its own output files into the current working directory.

> **Note on `0` vs `1`:** `0.Data_Cleaning.R` performs the *cleaning* pass and removes High-RITH samples with `HI < 0`; `1.Phenotype_Clustering_k2.R` re-runs the same SNF pipeline on the cleaned table using `set.seed(1234)` to produce the **final, reproducible** phenotype heatmap. Run `0` first, then `1`.

---

## 5. Script Descriptions

### `0.Data_Cleaning.R`
- Loads the raw feature table (`RITH_Data_Cleaned.csv`), classifies each feature by filter type (`original`, `log-sigma-*-mm-3D`, `wavelet-*`) and class (`firstorder`, `glcm`, `glrlm`, `glszm`, `gldm`, `ngtdm`).
- Keeps features with **MAD > 1**.
- Computes a pseudo heterogeneity tendency score `feature_means = rowMeans(X_sig)`.
- Computes the Spearman correlation of each feature with `feature_means`; features with `r > 0` → **HHR**, `r < 0` → **LHR**.
- SNF affinity matrix (`K = 20`, `sigma = 0.5`), spectral clustering `K = 2`.
- Labels the cluster with the higher mean feature value as **High RITH** (red), the other as **Low RITH** (blue).
- Computes `HI = Σ z(HHR) − Σ z(LHR)`.
- Removes High-RITH samples with `HI < 0`, writes the cleaned table back to `RITH_Data_Cleaned.csv`.
- Draws a direction-separated heatmap (`RITH_Direction_Separated_Heatmap.png`).

### `1.Phenotype_Clustering_k2.R`
Same pipeline as `0`, with `set.seed(1234)` and PDF output (`RITH_Direction_Separated_Heatmap.pdf`). **This is the version of record for the phenotype label.**

### `2.PCA_Dimensionality_Reduction.R`
PCA on MAD-selected features; writes `Standardized_PC_Scores.csv` and the combined `Final_RITH_Results_with_PC.csv`. Figures:
`Standardized_PCA_Phenotype.pdf`, `Standardized_PCA_HI.png`, `Standardized_PCA_Phenotype_HI_filtered.pdf`, `Standardized_PCA_HI_filtered.pdf`, `Standardized_PCA_Combined.png`.

### `3.HI_Boxplot.R`
Four HI visualisations: violin, box, density, and facetted histogram; also outputs `HI_Statistics_Summary.csv` (mean, SD, median, IQR, quartiles per group) and Cohen's *d*.

### `4.Feature_HI_Correlation.R`
Spearman correlation of every feature with HI; identifies the top 5 positive and top 5 negative features; outputs `All_Features_HI_Correlation.csv`, `Selected_Features_Detailed_Info.csv`, `HI_Scatter_Method1_FullNames.pdf` (95 % central-range filtered).

### `5.TMEimmune26.geneSurvival.R`
Kaplan–Meier progression-free-interval (PFI) analysis stratified by RITH phenotype. Writes `sur.RITH_Phenotype.tiff`. Group labels come from `RITH_Phenotype` (no re-clustering).

### `6.HI_Clinical_Variable_Correlation.R`
Spearman correlation of HI with **pT stage** and **ISUP grade**, plus scatterplots and Kruskal–Wallis boxplots. Writes `HI_correlation_plots.pdf` / `.png` and `HI_correlation_results.csv`.

### `7.Pathological_MITH.R`
Associates HI with each `Morph_*` feature:
- Pearson correlation per morphological feature → `01_HI_vs_Morphological_Features_Correlation.pdf`
- Violin plot per feature with Welch t-test → `02_Morphological_Features_Violin_Plots_with_Stats.pdf`
- Top-10 correlated / top-10 differential features → `05_*.pdf`, `06_*.pdf`
- Text report → `08_Analysis_Report.txt`

### `8.Biological_Analysis.R`
Bulk RNA-seq analysis (Treatment vs Control):
- `limma` + `voom` DE analysis → `DEG_results_limma.csv`
- Volcano plots (`volcano_basic`, `volcano_labeled`, `volcano_ECM_highlight`)
- PCA of samples (`PCA_plot.pdf`)
- GO / KEGG enrichment (`GO_all.csv`, `KEGG_all.csv`, `GO_barplot.pdf`, `GO_dotplot.pdf`, `KEGG_barplot.pdf`, `KEGG_dotplot.pdf`)
- GSEA against `c2.cp.kegg.Hs.symbols.gmt` → `GSEA_KEGG_all.csv`, `GSEA_ECM_pathway.pdf`

---

## 6. Key Methodological Notes

| Step | Method / Parameter |
|---|---|
| Feature filtering | `mad(x) > 1` on raw features |
| Pseudo heterogeneity tendency | `rowMeans(X_sig)` |
| Direction rule | Spearman `r > 0` → HHR; `r < 0` → LHR |
| Similarity | Euclidean distance, `affinityMatrix(K = 20, sigma = 0.5)` |
| Clustering | Spectral clustering, `K = 2` |
| High / Low label | Cluster with higher mean feature value = High RITH |
| Heterogeneity Index (HI) | `Σ z(HHR) − Σ z(LHR)` (column-wise Z-score) |
| Sample exclusion | High-RITH samples with `HI < 0` removed |
| Survival endpoint | Progression-free interval (`PFI.time`, `PFI`) |

> **Rationale for removing `HI < 0` in High-RITH:** If the phenotype is defined by *high* heterogeneity, samples that were assigned to the High-RITH cluster but have a negative HI contradict the definition; they are excluded to keep the phenotype internally consistent.

---

## 7. Outputs

### 7.1 Phenotype / HI files
| File | Description |
|---|---|
| `RITH_Data_Cleaned.csv` | Cleaned feature table with `RITH_Phenotype` and `Heterogeneity_Index` |
| `RITH_Phenotype_and_HI_Final.csv` | `Image`, `Cluster_Label`, `RITH_Phenotype`, `Feature_Mean`, `Heterogeneity_Index` |
| `Feature_Correlation_with_Heterogeneity.csv` | Per-feature Spearman correlation and direction |
| `Feature_Direction_Definitions_CorrelationBased.csv` | Per-feature direction definition |
| `Standardized_PC_Scores.csv` | Per-sample standardized PC1/PC2 |
| `Final_RITH_Results_with_PC.csv` | Final combined table (phenotype + HI + PCs) |
| `HI_Statistics_Summary.csv` | Mean / SD / median / quartiles of HI per group |
| `All_Features_HI_Correlation.csv` | Full feature × HI Spearman correlation table |
| `Selected_Features_Detailed_Info.csv` | Top-5 positive / negative features with direction and strength |

### 7.2 Figures
| File | Source script | Purpose |
|---|---|---|
| `RITH_Direction_Separated_Heatmap.png` / `.pdf` | `0`, `1` | Sample-by-feature heatmap split by direction |
| `Standardized_PCA_Phenotype.pdf` / `.png` | `2` | PCA coloured by phenotype |
| `Standardized_PCA_HI.png` / `.pdf` | `2` | PCA coloured by HI |
| `Standardized_PCA_Combined.png` | `2` | Side-by-side PCA |
| `HI_Violin_Plot.pdf` / `.png` | `3` | HI violin plot |
| `HI_Box_Plot.pdf` / `.png` | `3` | HI boxplot |
| `HI_Density_Plot.pdf` / `.png` | `3` | HI density plot |
| `HI_Histogram.pdf` / `.png` | `3` | HI histogram (facetted) |
| `HI_Scatter_Method1_FullNames.pdf` | `4` | HI vs top correlated features |
| `sur.RITH_Phenotype.tiff` | `5` | KM curve (PFI) |
| `HI_correlation_plots.pdf` / `.png` | `6` | HI vs pT / ISUP (scatter + box) |
| `Morphological_Analysis_Results/01–08` | `7` | Morphology vs HI (correlation, violin, stats, report) |
| `volcano_basic`, `volcano_labeled`, `volcano_ECM_highlight` | `8` | DE volcano plots |
| `PCA_plot.pdf` | `8` | Sample PCA (Treatment vs Control) |
| `GO_barplot`, `GO_dotplot`, `KEGG_barplot`, `KEGG_dotplot` | `8` | Enrichment plots |
| `GSEA_ECM_pathway.pdf` | `8` | GSEA plot for selected pathway |

---

## 8. How to Run

```bash
# 1. Place all input files in the expected sub-folders
# 2. Set working directory inside each script (or override with setwd())
# 3. Run in order:
Rscript 0.Data_Cleaning.R
Rscript 1.Phenotype_Clustering_k2.R
Rscript 2.PCA_Dimensionality_Reduction.R
Rscript 3.HI_Boxplot.R
Rscript 4.Feature_HI_Correlation.R
Rscript 5.TMEimmune26.geneSurvival.R
Rscript 6.HI_Clinical_Variable_Correlation.R
Rscript 7.Pathological_MITH.R
Rscript 8.Biological_Analysis.R
```

---

## 9. Notes, Caveats and Reproducibility

- **Random seeds.** `0.Data_Cleaning.R` uses `set.seed(123)`; `1.Phenotype_Clustering_k2.R` uses `set.seed(1234)`. Use the seed from `1` as the *canonical* run for the published phenotype.
- **Working-directory safety.** The `setwd()` calls in `2`–`8` are absolute and match the original analysis environment. **Edit them** if you move the project.
- **Feature naming.** Feature names follow the PyRadiomics convention: `{filter}_{class}_{name}` (e.g. `original_glcm_Idmn`, `wavelet-HLL_firstorder_Skewness`, `log-sigma-3-0-mm-3D_glrlm_GrayLevelNonUniformity`).
- **Direction rule is dataset-dependent.** The HHR/LHR split is derived from the *sign of the correlation with the pseudo-heterogeneity tendency*. If a different feature-extraction pipeline is used, re-verify the HHR/LHR split before interpreting HI as a biological readout.
- **High-RITH `HI < 0` removal** should be documented in the final manuscript; it is not a standard step and materially affects the sample size.
- **Pathways.** `8.Biological_Analysis.R` sets `ecm_pathway <- "KEGG_STARCH_AND_SUCROSE_METABOLISM"` in the current version; change it to the pathway of interest (e.g. `"KEGG_ECM_RECEPTOR_INTERACTION"`) before drawing the GSEA plot.
- **RNA-seq sample identifiers** are Chinese names and are **not** anonymised; they should be replaced with study IDs before public release.

---

## 10. Citation / Contact

If you use this pipeline, please cite the corresponding manuscript:

> *(Author list). (Year). Radiomic intratumoral heterogeneity predicts progression-free survival and immune microenvironment features in clear cell renal cell carcinoma. *Journal*, volume(issue), pages. doi:…

For questions or issues, please open an issue in this repository or contact the maintainer.
