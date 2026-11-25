# Cross-Cohort Transcriptomic Analysis of NRF2 PTM Regulators in Hepatocellular Carcinoma (HCCDB v2.0)

This repository contains the curated gene list, R analysis and visualization script, cohort audit log, and statistical summary tables supporting the cohort analysis (Figure 5 and Section 2.10) in our review article:

> **Post-Translational Modifications of NRF2: Mechanistic Basis and Therapeutic Implications for Hepatocellular Carcinoma**  
> Yumeng Wu, Wenjing Zhao, Haibin Yin, Guoxing Li, Chunyan Gu, Feng Xiao, Aiguo Shen, and Shudong Zhu  
> *Submitted to Frontiers in Cell and Developmental Biology (Specialty Section: Signaling)*

---

## 1. Data Source

All bulk transcriptomic expression matrices (`mRNA_level3`) and curated sample annotation files (`sample`) were downloaded directly from the **Hepatocellular Carcinoma Database v2.0 (HCCDB v2.0)** official portal:

- **Portal URL:** http://lifeome.net/database/hccdb2
- **Reference:** Jiang, Z., Wu, Y., Miao, Y., Deng, K., Yang, F., Xu, S., et al. (2024). HCCDB v2.0: Decompose Expression Variations by Single-cell RNA-seq and Spatial Transcriptomics in HCC. *Genomics, Proteomics & Bioinformatics*, 22(1), qzae011. https://doi.org/10.1093/gpbjnl/qzae011

Raw HCCDB v2.0 archives are publicly distributed by the database maintainers and are not redistributed in this repository. To avoid cross-project batch effects and control-group heterogeneity, no external normal liver repositories (e.g., GTEx) were pooled with HCC cohorts, and all statistical tests were performed strictly within each cohort comparing `HCC` tumor samples against `Adjacent` non-tumor liver samples.

---

## 2. Software Environment & Dependencies

- **R version:** R 4.6.1 (Windows 11 x64)
- **Core R packages:** `data.table`, `dplyr`, `tidyr`, `ggplot2`, `ComplexHeatmap`, `circlize`, `grid`, and base `stats` (exact versions of all loaded packages and platform details are recorded in `results/sessionInfo.txt`).

---

## 3. Repository Structure & File Descriptions

| File Path | Description |
| :--- | :--- |
| `master_gene_list.csv` | Curated panel of **49 NRF2-associated genes** (`NFE2L2`, 47 core regulators across the six PTM modules and PARP1 transcriptional coactivation, plus the indirect acetylation regulator `PCDH20`), including display labels, HCCDB query symbols (e.g., `OGA` queried via legacy symbol `MGEA5`), PTM categories, tiers, and published proteomic annotations (`Table 2.1` and Gao et al., 2019). |
| `run_hccdb_revision.R` | Single-file, reproducible R script that unpacks downloaded HCCDB v2.0 archives, audits cohort sample labels and expression scales, applies explicit $\log_2(x)$ transformation to linear-scale `HCCDB3` (GSE25097), performs within-cohort two-sided Wilcoxon rank-sum tests with Benjamini–Hochberg (BH) FDR adjustment, exports all summary CSVs, and generates Figure 5. |
| `results/cohort_check.csv` | Audit log covering all **25 bulk HCCDB v2.0 datasets**: platform type, raw sample-type counts, matched `n_HCC` and `n_Adjacent`, expression range (`min`/`median`/`max`), `log2_transformed` flag, and inclusion/exclusion rationale (**18 cohorts included**: 4 RNA-seq and 14 microarray, 2,482 HCC vs. 1,829 adjacent samples; **7 excluded**: 6 tumor-only cohorts and `HCCDB26` containing only non-adjacent `Normal` controls). Non-target labels (`Healthy`, `Cirrhotic`, `Normal`) within included cohorts were strictly excluded from the `Adjacent` group. |
| `results/gene_recurrence_summary.csv` | Gene-level cross-cohort recurrence table across all **18 included cohorts** (`n_detected_all`, `n_up_all`, `n_down_all`) and separately across the **17 independent non-TCGA cohorts** excluding `HCCDB15` / TCGA-LIHC (`n_detected_nontcga`, `n_up_nontcga`, `n_down_nontcga`, effect-size sensitivity counts at `padj < 0.05 & |delta_mean| >= 0.5`, and `median_delta_nontcga`), alongside `HCCDB15` (TCGA-LIHC) statistics. |
| `results/per_cohort_stats.csv` | Complete gene-by-cohort statistical output table containing per-cohort sample sizes, group means/medians, `delta_mean`, `delta_median`, raw Wilcoxon `p_value`, within-cohort BH-adjusted `padj`, and detection status (`NA` preserved for unprofiled genes). |
| `results/sessionInfo.txt` | Full `sessionInfo()` capture documenting the exact R version, operating system, and versions of all attached and loaded packages used to generate the published results. |

---

## 4. How to Reproduce

1. Download the bulk `sample.zip` and `mRNA_level3.zip` archives from the HCCDB v2.0 Download page (http://lifeome.net/database/hccdb2) into a local read-only directory.
2. Update the raw data input path at the top of `run_hccdb_revision.R` to point to your local HCCDB v2.0 folder.
3. Run the script from the repository root:
   ```bash
   Rscript run_hccdb_revision.R
