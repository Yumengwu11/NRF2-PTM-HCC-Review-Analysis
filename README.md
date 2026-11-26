# Cross-Cohort Transcriptomic and Sensitivity Analysis of NRF2 PTM Regulators in HCC (HCCDB v2.0)

This repository contains the curated gene list, single reproducible R analysis and visualization script, clinical metadata audit logs, and statistical summary tables supporting Figure 5A–C and Section 2.10 in our review article:

> **Post-Translational Modifications of NRF2: Mechanistic Basis and Therapeutic Implications for Hepatocellular Carcinoma**  
> Yumeng Wu, Wenjing Zhao, Haibin Yin, Guoxing Li, Chunyan Gu, Feng Xiao, Aiguo Shen, and Shudong Zhu  
> *Submitted to Frontiers in Cell and Developmental Biology (Specialty Section: Signaling)*

---

## 1. Data Source

All bulk transcriptomic expression matrices (`mRNA_level3`), sample annotation files (`sample`), and patient clinical tables (`patient`) were downloaded directly from the **Hepatocellular Carcinoma Database v2.0 (HCCDB v2.0)** official portal:

- **Portal URL:** http://lifeome.net/database/hccdb2
- **Reference:** Jiang, Z., Wu, Y., Miao, Y., Deng, K., Yang, F., Xu, S., et al. (2024). HCCDB v2.0: Decompose Expression Variations by Single-cell RNA-seq and Spatial Transcriptomics in HCC. *Genomics, Proteomics & Bioinformatics*, 22(1), qzae011. https://doi.org/10.1093/gpbjnl/qzae011

Raw HCCDB v2.0 archives are publicly distributed by the database maintainers and are not redistributed here. No external normal liver datasets (such as GTEx) were pooled with tumor cohorts, and all comparisons were executed strictly within each cohort (`HCC` vs. `Adjacent` liver, excluding non-target labels such as `Healthy`, `Cirrhotic`, or `Normal`).

---

## 2. Software Environment & Dependencies

- **R version:** R 4.6.1 (Windows 11 x64)
- **Core R packages:** `data.table`, `dplyr`, `tidyr`, `ggplot2`, `ComplexHeatmap`, `circlize`, `grid`, `cowplot` / `patchwork` (exact versions of all attached and loaded packages are recorded in `results/sessionInfo.txt`).

---

## 3. Repository Structure & File Descriptions

| File Path | Description |
| :--- | :--- |
| `master_gene_list.csv` | Curated panel of **49 NRF2-associated genes** (`NFE2L2`, 47 core regulators across the six PTM modules and PARP1 transcriptional coactivation, plus indirect acetylation regulator `PCDH20`), including display labels, legacy query symbols (e.g., `OGA`/`MGEA5`), PTM categories, tiers, and published proteomic annotations (`Table 2.1` and Gao et al., 2019). |
| `run_hccdb_analysis.R` | Single, self-contained R script that audits all 25 HCCDB v2.0 bulk cohorts, applies explicit `log2(x)` transformation to linear-scale `HCCDB3` (GSE25097), resolves legacy gene aliases, executes within-cohort Tumor vs. Adjacent Wilcoxon/BH analyses (**Figure 5A**), tumor-intrinsic 10-gene NRF2 ARE signature Spearman correlations (**Figure 5B**), and patient-paired / stage-stratified / etiology-stratified sensitivity analyses (**Figure 5C**), and exports all CSV tables and the composite Figure 5. |
| `results/cohort_check.csv` | Audit log across all **25 bulk HCCDB v2.0 datasets** (18 cohorts included: 4 RNA-seq and 14 microarray, 2,482 HCC vs. 1,829 adjacent samples; 7 excluded for lacking adjacent controls). |
| `results/gene_recurrence_summary.csv` | Gene-level cross-cohort recurrence table across all **18 included cohorts** and separately across the **17 independent non-TCGA cohorts** (excluding `HCCDB15`/TCGA-LIHC), including effect-size sensitivity counts (`padj < 0.05 & |delta_mean| >= 0.5`). |
| `results/per_cohort_stats.csv` | Complete gene-by-cohort statistical output (sample sizes, group means/medians, `delta_mean`, `delta_median`, Wilcoxon `p_value`, within-cohort BH `padj`, and `NA` status for unprofiled genes). |
| `results/supp_R4_1_per_cohort_are_corr.csv` & `supp_R4_1_are_corr_summary.csv` | Within-cohort and cross-cohort Spearman correlation results between each of the 49 regulators and the 10-gene NRF2 ARE target signature (`NQO1`, `GCLC`, `GCLM`, `TXNRD1`, `SRXN1`, `GSR`, `G6PD`, `ME1`, `AKR1C1`, `AKR1C3`; excluding `SQSTM1` and `NFE2L2`) across HCC tumor samples (supporting **Figure 5B**). |
| `results/supp_R4_5_clinical_metadata_audit.csv` | Cohort-by-cohort audit of patient-matching eligibility (11 cohorts, 1,003 matched pairs), clinical stage annotations (7 cohorts), etiology annotations (4 cohorts), and MASLD/NASH availability (unannotated in 17/18 cohorts; `n = 1` in `HCCDB30`). |
| `results/supp_R4_5_paired_per_cohort.csv` & `supp_R4_5_paired_summary.csv` | Patient-paired two-sided Wilcoxon signed-rank test results across the 11 paired cohorts (supporting **Figure 5C**, left panel). |
| `results/supp_R4_5_stratified_stats.csv` | Stage-stratified (`Early` [I–II / BCLC 0–A] and `Late` [III–IV / BCLC B–C] vs. `Adjacent`, plus `Late` vs. `Early`) and etiology-stratified (`HBV`, `HCV`, `Non-viral` vs. `Adjacent`) test results (supporting **Figure 5C**, middle and right panels). |
| `results/sessionInfo.txt` | Full `sessionInfo()` output capturing the exact R version and package versions used. |

---

## 4. How to Reproduce

1. Download the bulk `sample.zip`, `mRNA_level3.zip`, and `patient.zip` archives from the HCCDB v2.0 Download page (http://lifeome.net/database/hccdb2) into a local read-only directory.
2. Update the raw data input path at the top of `run_hccdb_analysis.R` to point to your local HCCDB v2.0 folder.
3. Run the single script from the repository root:
   ```bash
   Rscript run_hccdb_analysis.R
