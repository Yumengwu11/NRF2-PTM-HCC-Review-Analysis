#!/usr/bin/env Rscript
# ==============================================================================
# run_hccdb_analysis.R
# Unified, self-contained R pipeline for the NRF2 Review Revision:
#   1. Independent within-cohort screening & Tumor vs. Adjacent differential
#      expression across 18 HCCDB cohorts (49 NRF2 PTM regulators).
#   2. R4-1: Tumor-only Spearman correlation with a 10-gene classical NRF2 ARE
#      activity signature (excluding SQSTM1 & NFE2L2).
#   3. R4-5: Clinical metadata audit, patient-paired Wilcoxon signed-rank
#      sensitivity analysis (11 cohorts, n=1,003 pairs), and stage (7 cohorts)
#      / etiology (4 cohorts) stratified sensitivity analyses.
#   4. Pure-R composite Figure 5A-C (PDF + 300 dpi PNG).
#
# Hard constraints:
#   - Raw HCCDB zip directory is strictly READ-ONLY; reuses tmp_hccdb/ if present.
#   - Uses only pre-installed packages in D:/ruanjian/R-4.6.1/my_library.
#   - No cross-cohort expression matrix merging; check.names = FALSE throughout.
# ==============================================================================

# 0. Library path & package checks (DO NOT auto-install)
custom_lib <- "D:/ruanjian/R-4.6.1/my_library"
if (dir.exists(custom_lib)) {
  .libPaths(c(custom_lib, .libPaths()))
}

required_pkgs <- c(
  "data.table", "dplyr", "ggplot2", "ComplexHeatmap",
  "circlize", "grid", "cowplot", "ggrepel", "scales"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    sprintf("Missing required R package(s) in %s: %s. Please confirm before installing.",
            custom_lib, paste(missing_pkgs, collapse = ", ")),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(cowplot)
  library(ggrepel)
  library(scales)
})

cat("======================================================================\n")
cat("NRF2 Review Revision - Unified HCCDB Analysis & Figure 5A-C Pipeline\n")
cat("R version:", R.version.string, "\n")
cat("Library paths:", paste(.libPaths(), collapse = " | "), "\n")
cat("======================================================================\n\n")

# ------------------------------------------------------------------------------
# Paths & Constants
# ------------------------------------------------------------------------------
raw_hccdb_dir  <- "D:/project_TACE-PGCC-cell/04_Bioinformatics/PUBLIC_DATA/HCCDB"
gene_list_file <- "master_gene_list.csv"
tmp_dir        <- "tmp_hccdb"
res_dir        <- "results"

if (!file.exists(gene_list_file)) {
  stop("Gene list file not found in working directory: ", gene_list_file, call. = FALSE)
}
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

rnaseq_cohorts <- c("HCCDB15", "HCCDB18", "HCCDB25", "HCCDB30")

# 7 legacy alias fallback mappings (official symbol tried first, then aliases in order)
legacy_alias_map <- list(
  "OGA"     = c("MGEA5"),
  "KAT8"    = c("MYST1", "MOF"),
  "NAA10"   = c("ARD1A", "ARD1"),
  "DCAF11"  = c("WDR23"),
  "SKP1"    = c("SKP1A"),
  "USP17L2" = c("DUB3"),
  "SYVN1"   = c("HRD1")
)

# 10 classical NRF2 ARE downstream target genes (strictly excluding SQSTM1 & NFE2L2)
are_target_genes <- c("NQO1", "GCLC", "GCLM", "TXNRD1", "SRXN1", "GSR", "G6PD", "ME1", "AKR1C1", "AKR1C3")

# 9 focal genes for Panel C & console summary
core9_genes <- c("SQSTM1", "SAE1", "PARP1", "KEAP1", "CUL3", "RBX1", "NFE2L2", "DCAF11", "FYN")

# 11 highlighted genes for Panel B callouts
panelB_label_genes <- c("SQSTM1", "SIRT7", "NFE2L2", "PARP1", "SAE1", "KEAP1",
                        "RBX1", "CUL3", "USP11", "DCAF11", "FYN")

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------
get_cohort_num <- function(cid) {
  as.integer(sub("^HCCDB", "", cid))
}

fmt_table <- function(x) {
  x_clean <- ifelse(is.na(x) | !nzchar(x), "NA", as.character(x))
  tb <- table(x_clean)
  paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = "; ")
}

format_type_counts <- function(type_vec) {
  type_clean <- ifelse(is.na(type_vec) | !nzchar(type_vec), "<NA>", type_vec)
  tb <- table(type_clean)
  ord_names <- c(
    intersect(c("HCC", "Adjacent"), names(tb)),
    sort(setdiff(names(tb), c("HCC", "Adjacent")))
  )
  tb <- tb[ord_names]
  paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = "; ")
}

parse_sample_full <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  sp <- lapply(lines, function(x) trimws(strsplit(x, "\t", fixed = TRUE)[[1]]))
  keys <- vapply(sp, function(x) if (length(x) > 0) x[1] else "", character(1))

  sid_idx <- which(toupper(keys) == "SAMPLE_ID")[1]
  if (is.na(sid_idx)) {
    dt <- fread(filepath, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
    return(dt)
  }
  sids <- sp[[sid_idx]][-1]
  n_s <- length(sids)
  res <- list(SAMPLE_ID = sids)
  for (i in seq_along(sp)) {
    if (i == sid_idx) next
    k <- keys[i]
    if (!nzchar(k)) next
    vals <- sp[[i]][-1]
    if (length(vals) < n_s) {
      vals <- c(vals, rep(NA_character_, n_s - length(vals)))
    } else if (length(vals) > n_s) {
      vals <- vals[seq_len(n_s)]
    }
    vals[vals %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
    if (k %in% names(res)) k <- paste0(k, "_", i)
    res[[k]] <- vals
  }
  as.data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
}

parse_patient_table <- function(filepath) {
  if (is.na(filepath) || !file.exists(filepath)) return(NULL)
  lines <- readLines(filepath, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 2L) return(NULL)

  tok2 <- trimws(strsplit(lines[2], "\t", fixed = TRUE)[[1]])
  is_colwise <- grepl("^HCCDB-", tok2[1])
  if (is_colwise) {
    df <- fread(
      filepath, sep = "\t", header = TRUE, check.names = FALSE,
      data.table = FALSE, na.strings = c("NA", "NaN", "", "null", "NULL")
    )
    for (cn in colnames(df)) {
      if (is.character(df[[cn]])) df[[cn]] <- trimws(df[[cn]])
    }
    return(df)
  }

  sp <- lapply(lines, function(x) trimws(strsplit(x, "\t", fixed = TRUE)[[1]]))
  keys <- vapply(sp, function(x) if (length(x) > 0) x[1] else "", character(1))
  pid_idx <- which(toupper(keys) == "PATIENT_ID")[1]
  if (is.na(pid_idx)) pid_idx <- 1L
  n_p <- length(sp[[pid_idx]]) - 1L
  res <- list()
  for (i in seq_along(sp)) {
    k <- keys[i]
    if (!nzchar(k)) next
    vals <- sp[[i]][-1]
    if (length(vals) < n_p) {
      vals <- c(vals, rep(NA_character_, n_p - length(vals)))
    } else if (length(vals) > n_p) {
      vals <- vals[seq_len(n_p)]
    }
    vals[vals %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
    if (k %in% names(res)) k <- paste0(k, "_", i)
    res[[k]] <- vals
  }
  as.data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
}

read_mrna_raw <- function(filepath) {
  dt <- fread(
    filepath, sep = "\t", header = TRUE, check.names = FALSE,
    data.table = FALSE, na.strings = c("NA", "NaN", "", "null", "NULL")
  )
  cnames <- colnames(dt)
  cnames_low <- tolower(trimws(cnames))
  sym_col <- which(cnames_low %in% c("symbol", "gene_symbol", "genesymbol", "gene"))[1]
  entrez_col <- which(grepl("^entrez", cnames_low))[1]
  if (is.na(sym_col)) {
    if (ncol(dt) >= 2 && is.numeric(dt[[1]]) && is.character(dt[[2]])) {
      entrez_col <- 1L; sym_col <- 2L
    } else {
      sym_col <- 1L
    }
  }
  meta_cols <- unique(na.omit(c(entrez_col, sym_col)))
  sample_cols <- setdiff(seq_len(ncol(dt)), meta_cols)

  symbols <- trimws(as.character(dt[[sym_col]]))
  mat <- as.matrix(dt[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- trimws(cnames[sample_cols])
  list(symbols = symbols, mat = mat)
}

match_gene_row <- function(disp_sym, query_sym, syms, all_row_means, alias_map = legacy_alias_map) {
  alias_syms <- if (disp_sym %in% names(alias_map)) alias_map[[disp_sym]] else character(0)
  cand_priority <- unique(c(disp_sym, query_sym, alias_syms))
  hit_idx <- integer(0)
  for (cand in cand_priority) {
    idx <- which(syms == cand)
    if (length(idx) > 0L) {
      hit_idx <- idx
      break
    }
  }
  if (length(hit_idx) == 0L) {
    return(list(row_idx = NA_integer_, matched_symbol = NA_character_, n_hits = 0L))
  }
  if (length(hit_idx) > 1L) {
    mvals <- all_row_means[hit_idx]
    best_pos <- if (all(is.na(mvals))) 1L else which.max(mvals)
    chosen <- hit_idx[best_pos]
  } else {
    chosen <- hit_idx[1L]
  }
  list(row_idx = chosen, matched_symbol = syms[chosen], n_hits = length(hit_idx))
}

map_ajcc_stage <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s <- sub("^STAGE\\s+", "", s)
  out <- rep(NA_character_, length(s))
  out[s %in% c("I", "IA", "IB", "IC", "1", "1A", "1B", "II", "IIA", "IIB", "IIC", "2", "2A", "2B")] <- "Early_I_II"
  out[s %in% c("III", "IIIA", "IIIB", "IIIC", "3", "3A", "3B", "3C", "IV", "IVA", "IVB", "4", "4A", "4B")] <- "Late_III_IV"
  out
}

map_bclc_stage <- function(x) {
  s <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(s))
  out[s %in% c("0", "A")] <- "Early_BCLC_0_A"
  out[s %in% c("B", "C")] <- "Late_BCLC_B_C"
  out
}

run_unpaired_panel <- function(mat, gene_matches, master_genes, cols_g1, cols_g2,
                               cid, domain, subgroup_label, comparison_label,
                               ref_type, mapping_note) {
  rows <- vector("list", nrow(master_genes))
  for (g in seq_len(nrow(master_genes))) {
    ridx <- gene_matches$row_idx[g]
    msym <- gene_matches$matched_symbol[g]
    if (is.na(ridx)) {
      rows[[g]] <- data.frame(
        cohort = cid, is_tcga = (cid == "HCCDB15"), strat_domain = domain,
        subgroup = subgroup_label, comparison = comparison_label, ref_type = ref_type,
        mapping_note = mapping_note,
        display_label = master_genes$display_label[g], matched_symbol = NA_character_,
        ptm_category = master_genes$ptm_category[g],
        n_g1 = length(cols_g1), n_g2 = length(cols_g2),
        mean_g1 = NA_real_, mean_g2 = NA_real_, delta_mean = NA_real_,
        median_g1 = NA_real_, median_g2 = NA_real_, delta_median = NA_real_,
        p_value = NA_real_, padj = NA_real_, status = "not_detected",
        stringsAsFactors = FALSE
      )
      next
    }
    x1 <- as.numeric(mat[ridx, cols_g1]); x1 <- x1[!is.na(x1)]
    x2 <- as.numeric(mat[ridx, cols_g2]); x2 <- x2[!is.na(x2)]
    m1 <- if (length(x1) > 0) mean(x1) else NA_real_
    m2 <- if (length(x2) > 0) mean(x2) else NA_real_
    md1 <- if (length(x1) > 0) median(x1) else NA_real_
    md2 <- if (length(x2) > 0) median(x2) else NA_real_
    dm <- if (length(x1) > 0 && length(x2) > 0) (m1 - m2) else NA_real_
    dmd <- if (length(x1) > 0 && length(x2) > 0) (md1 - md2) else NA_real_

    if (length(x1) < 3L || length(x2) < 3L) {
      st <- "insufficient_nonNA"; pval <- NA_real_
    } else if (length(unique(c(x1, x2))) <= 1L || stats::var(c(x1, x2)) == 0) {
      st <- "zero_variance"; pval <- NA_real_
    } else {
      wt <- suppressWarnings(stats::wilcox.test(x1, x2, alternative = "two.sided", exact = FALSE))
      st <- "tested"; pval <- as.numeric(wt$p.value)
    }

    rows[[g]] <- data.frame(
      cohort = cid, is_tcga = (cid == "HCCDB15"), strat_domain = domain,
      subgroup = subgroup_label, comparison = comparison_label, ref_type = ref_type,
      mapping_note = mapping_note,
      display_label = master_genes$display_label[g], matched_symbol = msym,
      ptm_category = master_genes$ptm_category[g],
      n_g1 = length(x1), n_g2 = length(x2),
      mean_g1 = m1, mean_g2 = m2, delta_mean = dm,
      median_g1 = md1, median_g2 = md2, delta_median = dmd,
      p_value = pval, padj = NA_real_, status = st,
      stringsAsFactors = FALSE
    )
  }
  df <- bind_rows(rows)
  df$padj <- stats::p.adjust(df$p_value, method = "BH")
  df$direction <- ifelse(
    is.na(df$padj), NA_character_,
    ifelse(df$padj < 0.05 & df$delta_mean > 0, "Up",
           ifelse(df$padj < 0.05 & df$delta_mean < 0, "Down", "NS"))
  )
  df
}

# ==============================================================================
# STEP 1: Reuse tmp_hccdb/ (or unzip if empty) & Discover All 25 Cohorts
# ==============================================================================
existing_tmp_files <- list.files(tmp_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)

if (length(existing_tmp_files) > 0L) {
  cat(sprintf("[Step 1] Reusing %d existing extracted .txt files in '%s/' (read-only source untouched).\n",
              length(existing_tmp_files), tmp_dir))
} else {
  if (!dir.exists(raw_hccdb_dir)) {
    stop("Read-only HCCDB directory not found: ", raw_hccdb_dir, call. = FALSE)
  }
  zip_files <- sort(list.files(raw_hccdb_dir, pattern = "\\.zip$", full.names = TRUE))
  cat(sprintf("[Step 1] Extracting %d .zip archives into '%s/'...\n", length(zip_files), tmp_dir))
  for (zf in zip_files) unzip(zf, exdir = tmp_dir, overwrite = TRUE)
  existing_tmp_files <- list.files(tmp_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
}

sample_files_vec <- existing_tmp_files[grepl("\\.sample\\.txt$", basename(existing_tmp_files))]
sample_tbl <- data.frame(
  cohort = sub("\\.sample\\.txt$", "", basename(sample_files_vec)),
  sample_path = sample_files_vec,
  stringsAsFactors = FALSE
)

mrna_candidates <- existing_tmp_files[!grepl("\\.(sample|patient)\\.txt$", basename(existing_tmp_files))]
mrna_tbl <- bind_rows(lapply(mrna_candidates, function(fp) {
  bn <- basename(fp)
  if (grepl("^HCCDB[0-9]+_mRNA_level3\\.txt$", bn)) {
    cid <- sub("_mRNA_level3\\.txt$", "", bn)
  } else {
    hdr <- readLines(fp, n = 1L, warn = FALSE)
    toks <- trimws(strsplit(hdr, "\t", fixed = TRUE)[[1]])
    sid_hit <- toks[grepl("^HCCDB-[0-9]+\\.", toks)][1]
    cid <- if (!is.na(sid_hit)) sub("^HCCDB-([0-9]+)\\..*$", "HCCDB\\1", sid_hit) else NA_character_
  }
  data.frame(cohort = cid, mrna_path = fp, stringsAsFactors = FALSE)
})) %>% filter(!is.na(cohort))

cohort_files <- full_join(sample_tbl, mrna_tbl, by = "cohort") %>%
  mutate(
    cohort_num = get_cohort_num(cohort),
    patient_path = file.path(tmp_dir, paste0(cohort, ".patient.txt"))
  ) %>%
  arrange(cohort_num)

# On-demand extraction of any missing .patient.txt if .patient.zip exists in read-only source
for (i in seq_len(nrow(cohort_files))) {
  if (!file.exists(cohort_files$patient_path[i])) {
    pzip <- file.path(raw_hccdb_dir, paste0(cohort_files$cohort[i], ".patient.zip"))
    if (file.exists(pzip)) {
      unzip(pzip, exdir = tmp_dir, overwrite = TRUE)
    } else {
      cohort_files$patient_path[i] <- NA_character_
    }
  }
}

master_genes <- read.csv(gene_list_file, stringsAsFactors = FALSE, check.names = FALSE)
master_genes$display_label <- trimws(master_genes$display_label)
master_genes$hccdb_query_symbol <- trimws(master_genes$hccdb_query_symbol)
cat(sprintf("  Mapped %d cohorts; loaded %d target genes from '%s'.\n\n",
            nrow(cohort_files), nrow(master_genes), gene_list_file))

# ==============================================================================
# STEP 2: Independent Within-Cohort Screening, DE, R4-1 (ARE), & R4-5 Analyses
# ==============================================================================
cat("[Step 2] Running independent within-cohort DE, R4-1 ARE correlation, and R4-5 sensitivity...\n")

cohort_check_list     <- list()
per_cohort_stats_list <- list()
r41_per_cohort_list   <- list()
audit_list            <- list()
paired_list           <- list()
strat_list            <- list()

for (i in seq_len(nrow(cohort_files))) {
  cid   <- cohort_files$cohort[i]
  spath <- cohort_files$sample_path[i]
  ppath <- cohort_files$patient_path[i]
  mpath <- cohort_files$mrna_path[i]
  plat  <- ifelse(cid %in% rnaseq_cohorts, "RNA-seq", "Microarray")

  s_df <- parse_sample_full(spath)
  raw_counts_str <- format_type_counts(s_df$TYPE)

  mrna_raw <- read_mrna_raw(mpath)
  mat  <- mrna_raw$mat
  syms <- mrna_raw$symbols
  mrna_sids <- colnames(mat)

  # Strict exact HCC and Adjacent matching
  s_mrna   <- s_df[s_df$SAMPLE_ID %in% mrna_sids, , drop = FALSE]
  hcc_cols <- mrna_sids[mrna_sids %in% s_df$SAMPLE_ID[s_df$TYPE == "HCC"]]
  adj_cols <- mrna_sids[mrna_sids %in% s_df$SAMPLE_ID[s_df$TYPE == "Adjacent"]]
  n_hcc <- length(hcc_cols)
  n_adj <- length(adj_cols)

  expr_min_raw <- suppressWarnings(min(mat, na.rm = TRUE))
  expr_med_raw <- suppressWarnings(median(mat, na.rm = TRUE))
  expr_max_raw <- suppressWarnings(max(mat, na.rm = TRUE))

  did_log2 <- FALSE
  scale_note <- ""
  if (cid == "HCCDB3") {
    if (is.finite(expr_min_raw) && expr_min_raw > 0) {
      mat <- log2(mat)
      scale_note <- "Linear scale (GSE25097): explicit log2(x) applied"
    } else {
      mat <- log2(mat + 1)
      scale_note <- "Linear scale (GSE25097, <=0 present): explicit log2(x+1) applied"
    }
    did_log2 <- TRUE
  } else if (is.finite(expr_max_raw) && is.finite(expr_min_raw) && expr_max_raw > 100 && expr_min_raw >= 0) {
    mat <- log2(mat + 1)
    did_log2 <- TRUE
    scale_note <- "Auto log2(x+1) applied (max > 100 & min >= 0)"
  }

  included <- (n_hcc >= 3L && n_adj >= 3L)
  if (included) {
    other_types <- setdiff(unique(s_df$TYPE), c("HCC", "Adjacent"))
    notes_vec <- character(0)
    if (length(other_types) > 0) {
      notes_vec <- c(notes_vec, sprintf("ignored non-target sample types: %s", paste(other_types, collapse = ", ")))
    }
    if (nzchar(scale_note)) notes_vec <- c(notes_vec, scale_note)
    excl_reason <- if (length(notes_vec) > 0) sprintf("Included (%s)", paste(notes_vec, collapse = "; ")) else "Included"
  } else {
    if (n_adj == 0L) {
      non_hcc <- setdiff(unique(s_df$TYPE), "HCC")
      excl_reason <- if (length(non_hcc) > 0) {
        sprintf("No Adjacent samples (has non-Adjacent control label(s): %s; not mixed)", paste(non_hcc, collapse = ", "))
      } else {
        "No Adjacent samples (HCC tumor-only cohort)"
      }
    } else {
      excl_reason <- sprintf("Insufficient samples (HCC=%d, Adj=%d)", n_hcc, n_adj)
    }
  }

  cohort_check_list[[cid]] <- data.frame(
    cohort = cid, platform = plat,
    sample_file = basename(spath), mrna_file = basename(mpath),
    raw_type_counts = raw_counts_str, n_mrna_samples_total = ncol(mat),
    n_HCC = n_hcc, n_Adjacent = n_adj,
    expr_min = round(expr_min_raw, 6), expr_median = round(expr_med_raw, 6), expr_max = round(expr_max_raw, 6),
    log2_transformed = did_log2, included = included, exclusion_reason = excl_reason,
    stringsAsFactors = FALSE
  )

  if (!included) next

  all_row_means <- rowMeans(mat, na.rm = TRUE)
  gm_list <- lapply(seq_len(nrow(master_genes)), function(g) {
    match_gene_row(master_genes$display_label[g], master_genes$hccdb_query_symbol[g], syms, all_row_means)
  })
  gene_matches <- data.frame(
    row_idx = vapply(gm_list, function(x) x$row_idx, integer(1)),
    matched_symbol = vapply(gm_list, function(x) x$matched_symbol, character(1)),
    stringsAsFactors = FALSE
  )

  # --- 2A. Main Tumor vs. Adjacent Wilcoxon + within-cohort BH ---
  gene_rows <- vector("list", nrow(master_genes))
  for (g in seq_len(nrow(master_genes))) {
    ridx <- gene_matches$row_idx[g]
    msym <- gene_matches$matched_symbol[g]
    if (is.na(ridx)) {
      gene_rows[[g]] <- data.frame(
        cohort = cid, is_tcga = (cid == "HCCDB15"), platform = plat,
        display_label = master_genes$display_label[g],
        hccdb_query_symbol = master_genes$hccdb_query_symbol[g],
        matched_symbol = NA_character_,
        ptm_category = master_genes$ptm_category[g], tier = master_genes$tier[g],
        n_HCC_cohort = n_hcc, n_Adjacent_cohort = n_adj,
        n_HCC_nonNA = 0L, n_Adjacent_nonNA = 0L,
        mean_HCC = NA_real_, mean_Adjacent = NA_real_, delta_mean = NA_real_,
        median_HCC = NA_real_, median_Adjacent = NA_real_, delta_median = NA_real_,
        p_value = NA_real_, padj = NA_real_, status = "not_detected",
        stringsAsFactors = FALSE
      )
      next
    }
    xh <- as.numeric(mat[ridx, hcc_cols]); xh <- xh[!is.na(xh)]
    xa <- as.numeric(mat[ridx, adj_cols]); xa <- xa[!is.na(xa)]
    mh <- if (length(xh) > 0) mean(xh) else NA_real_
    ma <- if (length(xa) > 0) mean(xa) else NA_real_
    mdh <- if (length(xh) > 0) median(xh) else NA_real_
    mda <- if (length(xa) > 0) median(xa) else NA_real_
    dm <- if (length(xh) > 0 && length(xa) > 0) (mh - ma) else NA_real_
    dmd <- if (length(xh) > 0 && length(xa) > 0) (mdh - mda) else NA_real_

    if (length(xh) < 3L || length(xa) < 3L) {
      st <- "insufficient_nonNA"; pv <- NA_real_
    } else if (length(unique(c(xh, xa))) <= 1L || stats::var(c(xh, xa)) == 0) {
      st <- "zero_variance"; pv <- NA_real_
    } else {
      wt <- suppressWarnings(stats::wilcox.test(xh, xa, alternative = "two.sided", exact = FALSE))
      st <- "tested"; pv <- as.numeric(wt$p.value)
    }

    gene_rows[[g]] <- data.frame(
      cohort = cid, is_tcga = (cid == "HCCDB15"), platform = plat,
      display_label = master_genes$display_label[g],
      hccdb_query_symbol = master_genes$hccdb_query_symbol[g],
      matched_symbol = msym,
      ptm_category = master_genes$ptm_category[g], tier = master_genes$tier[g],
      n_HCC_cohort = n_hcc, n_Adjacent_cohort = n_adj,
      n_HCC_nonNA = length(xh), n_Adjacent_nonNA = length(xa),
      mean_HCC = mh, mean_Adjacent = ma, delta_mean = dm,
      median_HCC = mdh, median_Adjacent = mda, delta_median = dmd,
      p_value = pv, padj = NA_real_, status = st,
      stringsAsFactors = FALSE
    )
  }
  cdf_de <- bind_rows(gene_rows)
  cdf_de$padj <- stats::p.adjust(cdf_de$p_value, method = "BH")
  cdf_de$direction <- ifelse(
    is.na(cdf_de$padj), NA_character_,
    ifelse(cdf_de$padj < 0.05 & cdf_de$delta_mean > 0, "Up",
           ifelse(cdf_de$padj < 0.05 & cdf_de$delta_mean < 0, "Down", "NS"))
  )
  per_cohort_stats_list[[cid]] <- cdf_de

  # --- 2B. R4-1: Tumor-only 10-gene NRF2 ARE score Spearman correlation ---
  are_idx <- vapply(are_target_genes, function(ag) {
    match_gene_row(ag, ag, syms, all_row_means, alias_map = list())$row_idx
  }, integer(1))
  are_found_names <- are_target_genes[!is.na(are_idx)]
  are_found_idx   <- are_idx[!is.na(are_idx)]
  valid_are_z <- list()
  for (k in seq_along(are_found_idx)) {
    v <- as.numeric(mat[are_found_idx[k], hcc_cols])
    sd_v <- stats::sd(v, na.rm = TRUE)
    if (is.finite(sd_v) && sd_v > 0) {
      valid_are_z[[are_found_names[k]]] <- (v - mean(v, na.rm = TRUE)) / sd_v
    }
  }
  n_are_used <- length(valid_are_z)
  if (n_are_used >= 5L && length(hcc_cols) >= 5L) {
    are_score <- colMeans(do.call(rbind, valid_are_z), na.rm = TRUE)
    names(are_score) <- hcc_cols
    corr_rows <- vector("list", nrow(master_genes))
    for (g in seq_len(nrow(master_genes))) {
      ridx <- gene_matches$row_idx[g]
      msym <- gene_matches$matched_symbol[g]
      if (is.na(ridx)) {
        corr_rows[[g]] <- data.frame(
          cohort = cid, is_tcga = (cid == "HCCDB15"), platform = plat,
          n_HCC_tumors = length(hcc_cols), n_ARE_genes_used = n_are_used,
          ARE_genes_used = paste(names(valid_are_z), collapse = ";"),
          display_label = master_genes$display_label[g], matched_symbol = NA_character_,
          ptm_category = master_genes$ptm_category[g], tier = master_genes$tier[g],
          n_valid_pairs = 0L, rho = NA_real_, p_value = NA_real_, padj = NA_real_,
          status = "not_detected", stringsAsFactors = FALSE
        )
        next
      }
      gv <- as.numeric(mat[ridx, hcc_cols])
      ok <- !is.na(gv) & !is.na(are_score)
      gv_ok <- gv[ok]; sc_ok <- are_score[ok]
      if (length(gv_ok) < 5L) {
        st <- "insufficient_nonNA"; rho_v <- NA_real_; pv <- NA_real_
      } else if (stats::sd(gv_ok) == 0 || stats::sd(sc_ok) == 0) {
        st <- "zero_variance"; rho_v <- NA_real_; pv <- NA_real_
      } else {
        ct <- suppressWarnings(stats::cor.test(gv_ok, sc_ok, method = "spearman", exact = FALSE))
        st <- "tested"; rho_v <- as.numeric(ct$estimate); pv <- as.numeric(ct$p.value)
      }
      corr_rows[[g]] <- data.frame(
        cohort = cid, is_tcga = (cid == "HCCDB15"), platform = plat,
        n_HCC_tumors = length(hcc_cols), n_ARE_genes_used = n_are_used,
        ARE_genes_used = paste(names(valid_are_z), collapse = ";"),
        display_label = master_genes$display_label[g], matched_symbol = msym,
        ptm_category = master_genes$ptm_category[g], tier = master_genes$tier[g],
        n_valid_pairs = length(gv_ok), rho = rho_v, p_value = pv, padj = NA_real_,
        status = st, stringsAsFactors = FALSE
      )
    }
    cdf_r41 <- bind_rows(corr_rows)
    cdf_r41$padj <- stats::p.adjust(cdf_r41$p_value, method = "BH")
    cdf_r41$corr_dir <- ifelse(
      is.na(cdf_r41$padj), NA_character_,
      ifelse(cdf_r41$padj < 0.05 & cdf_r41$rho > 0, "Pos",
             ifelse(cdf_r41$padj < 0.05 & cdf_r41$rho < 0, "Neg", "NS"))
    )
    r41_per_cohort_list[[cid]] <- cdf_r41
  }

  # --- 2C. R4-5: Clinical audit, patient-paired & stage/etiology stratification ---
  p_df <- parse_patient_table(ppath)
  pid_col_s <- intersect(c("PATIENT_ID", "PATIENTA_ID", "PATIENT"), colnames(s_mrna))[1]
  s_mrna$.PID <- as.character(s_mrna[[pid_col_s]])
  hcc_pids <- unique(s_mrna$.PID[s_mrna$TYPE == "HCC" & !is.na(s_mrna$.PID)])
  adj_pids <- unique(s_mrna$.PID[s_mrna$TYPE == "Adjacent" & !is.na(s_mrna$.PID)])
  matched_pids <- intersect(hcc_pids, adj_pids)
  n_pairs <- length(matched_pids)

  multi_hcc_pts <- character(0); multi_adj_pts <- character(0)
  if (n_pairs > 0L) {
    tb_h <- table(s_mrna$.PID[s_mrna$TYPE == "HCC" & s_mrna$.PID %in% matched_pids])
    tb_a <- table(s_mrna$.PID[s_mrna$TYPE == "Adjacent" & s_mrna$.PID %in% matched_pids])
    multi_hcc_pts <- names(tb_h[tb_h > 1L])
    multi_adj_pts <- names(tb_a[tb_a > 1L])
  }
  multi_note <- if (length(multi_hcc_pts) > 0L || length(multi_adj_pts) > 0L) {
    sprintf("Averaged multi-sample patients: %d multi-HCC (%s); %d multi-Adj (%s)",
            length(multi_hcc_pts), paste(multi_hcc_pts, collapse = ","),
            length(multi_adj_pts), paste(multi_adj_pts, collapse = ","))
  } else {
    "1:1 unique pairs"
  }

  s_clin <- s_mrna
  if (!is.null(p_df) && "PATIENT_ID" %in% colnames(p_df)) {
    p_dedup <- p_df[!duplicated(p_df$PATIENT_ID), , drop = FALSE]
    extra_cols <- setdiff(colnames(p_dedup), colnames(s_clin))
    if (length(extra_cols) > 0L) {
      m_idx <- match(s_clin$.PID, p_dedup$PATIENT_ID)
      for (ec in extra_cols) s_clin[[ec]] <- p_dedup[[ec]][m_idx]
    }
  }

  stage_field <- NA_character_; stage_raw <- NA_character_; has_stage <- FALSE
  if (cid == "HCCDB4" && "AJCC_STAGE" %in% colnames(s_clin)) {
    has_stage <- TRUE; stage_field <- "sample:AJCC_STAGE"; stage_raw <- fmt_table(s_clin$AJCC_STAGE[s_clin$TYPE == "HCC"])
  } else if (cid == "HCCDB13" && "BCLC_STAGE" %in% colnames(s_clin)) {
    has_stage <- TRUE; stage_field <- "patient:BCLC_STAGE"; stage_raw <- fmt_table(s_clin$BCLC_STAGE[s_clin$TYPE == "HCC"])
  } else if (cid == "HCCDB15" && "AJCC_STAGE1" %in% colnames(s_clin)) {
    has_stage <- TRUE; stage_field <- "patient:AJCC_STAGE1"; stage_raw <- fmt_table(s_clin$AJCC_STAGE1[s_clin$TYPE == "HCC"])
  } else if (cid == "HCCDB18" && "TNM_STAGE_T1" %in% colnames(s_clin)) {
    has_stage <- TRUE; stage_field <- "patient:TNM_STAGE_T1"; stage_raw <- fmt_table(s_clin$TNM_STAGE_T1[s_clin$TYPE == "HCC"])
  } else if (cid == "HCCDB24" && "AJCC_Stage" %in% colnames(s_clin)) {
    has_stage <- TRUE; stage_field <- "patient:AJCC_Stage"; stage_raw <- fmt_table(s_clin$AJCC_Stage[s_clin$TYPE == "HCC"])
  } else if (cid == "HCCDB25" && "AJCC_Stage1" %in% colnames(s_clin)) {
    has_stage <- TRUE; stage_field <- "patient:AJCC_Stage1"; stage_raw <- fmt_table(s_clin$AJCC_Stage1[s_clin$TYPE == "HCC"])
  } else if (cid == "HCCDB30" && "AJCC_STAGE" %in% colnames(s_clin)) {
    has_stage <- TRUE; stage_field <- "patient:AJCC_STAGE"; stage_raw <- fmt_table(s_clin$AJCC_STAGE[s_clin$TYPE == "HCC"])
  }

  has_etiol <- FALSE; etiol_field <- NA_character_; etiol_raw <- NA_character_
  masld_flag <- FALSE; masld_note <- "No MASLD/NAFLD/NASH annotation"
  if (cid == "HCCDB18" && "VIRUS" %in% colnames(s_clin)) {
    has_etiol <- TRUE; etiol_field <- "patient:VIRUS (+ALCOHOL)"; etiol_raw <- fmt_table(s_clin$VIRUS[s_clin$TYPE == "HCC"])
  } else if (cid == "HCCDB25") {
    has_etiol <- TRUE; etiol_field <- "cohort_design:CHCC-HBV (Gao 2019)"
    etiol_raw <- sprintf("HBV:%d (100%% HBV-related HCC cohort)", length(hcc_cols))
  } else if (cid == "HCCDB27" && "HBV_STATE" %in% colnames(s_clin)) {
    has_etiol <- TRUE; etiol_field <- "sample:HBV_STATE"; etiol_raw <- fmt_table(s_clin$HBV_STATE[s_clin$TYPE == "HCC"])
  } else if (cid == "HCCDB30" && "Risk factor" %in% colnames(s_clin)) {
    has_etiol <- TRUE; etiol_field <- "patient:Risk factor"; etiol_raw <- fmt_table(s_clin[["Risk factor"]][s_clin$TYPE == "HCC"])
    masld_flag <- TRUE
    masld_note <- "NASH:1 present in Risk factor (n=1; insufficient for standalone MASLD test; pooled into Non-viral)"
  }

  audit_list[[cid]] <- data.frame(
    cohort = cid, platform = plat, has_patient_file = !is.null(p_df),
    n_HCC_mrna = n_hcc, n_Adjacent_mrna = n_adj, patient_id_key = pid_col_s,
    n_matched_pairs = n_pairs, paired_analysis_eligible = (n_pairs >= 10L),
    multi_sample_patients_note = multi_note,
    has_stage = has_stage, stage_field = stage_field, stage_raw_counts_HCC = stage_raw,
    has_etiology = has_etiol, etiology_field = etiol_field, etiology_raw_counts_HCC = etiol_raw,
    MASLD_annotated = masld_flag, MASLD_note = masld_note,
    stringsAsFactors = FALSE
  )

  # Patient-paired Wilcoxon signed-rank test (n_pairs >= 10)
  if (n_pairs >= 10L) {
    hcc_map <- split(s_mrna$SAMPLE_ID[s_mrna$TYPE == "HCC" & s_mrna$.PID %in% matched_pids],
                     s_mrna$.PID[s_mrna$TYPE == "HCC" & s_mrna$.PID %in% matched_pids])
    adj_map <- split(s_mrna$SAMPLE_ID[s_mrna$TYPE == "Adjacent" & s_mrna$.PID %in% matched_pids],
                     s_mrna$.PID[s_mrna$TYPE == "Adjacent" & s_mrna$.PID %in% matched_pids])
    p_rows <- vector("list", nrow(master_genes))
    for (g in seq_len(nrow(master_genes))) {
      ridx <- gene_matches$row_idx[g]
      msym <- gene_matches$matched_symbol[g]
      if (is.na(ridx)) {
        p_rows[[g]] <- data.frame(
          cohort = cid, is_tcga = (cid == "HCCDB15"), platform = plat,
          display_label = master_genes$display_label[g], matched_symbol = NA_character_,
          ptm_category = master_genes$ptm_category[g], tier = master_genes$tier[g],
          n_pairs_cohort = n_pairs, n_pairs_valid = 0L,
          mean_HCC_paired = NA_real_, mean_Adj_paired = NA_real_,
          mean_paired_delta = NA_real_, median_paired_delta = NA_real_,
          p_value = NA_real_, padj = NA_real_, status = "not_detected",
          stringsAsFactors = FALSE
        )
        next
      }
      row_vec <- mat[ridx, ]
      xh <- vapply(matched_pids, function(pid) mean(as.numeric(row_vec[hcc_map[[pid]]]), na.rm = TRUE), numeric(1))
      xa <- vapply(matched_pids, function(pid) mean(as.numeric(row_vec[adj_map[[pid]]]), na.rm = TRUE), numeric(1))
      ok <- is.finite(xh) & is.finite(xa)
      xh_ok <- xh[ok]; xa_ok <- xa[ok]; diffs <- xh_ok - xa_ok
      if (length(diffs) < 10L) {
        st <- "insufficient_pairs"; pv <- NA_real_
      } else if (all(diffs == 0) || stats::var(diffs) == 0) {
        st <- "zero_variance"; pv <- NA_real_
      } else {
        wt <- suppressWarnings(stats::wilcox.test(xh_ok, xa_ok, paired = TRUE, alternative = "two.sided", exact = FALSE))
        st <- "tested"; pv <- as.numeric(wt$p.value)
      }
      p_rows[[g]] <- data.frame(
        cohort = cid, is_tcga = (cid == "HCCDB15"), platform = plat,
        display_label = master_genes$display_label[g], matched_symbol = msym,
        ptm_category = master_genes$ptm_category[g], tier = master_genes$tier[g],
        n_pairs_cohort = n_pairs, n_pairs_valid = length(diffs),
        mean_HCC_paired = if (length(xh_ok) > 0) mean(xh_ok) else NA_real_,
        mean_Adj_paired = if (length(xa_ok) > 0) mean(xa_ok) else NA_real_,
        mean_paired_delta = if (length(diffs) > 0) mean(diffs) else NA_real_,
        median_paired_delta = if (length(diffs) > 0) median(diffs) else NA_real_,
        p_value = pv, padj = NA_real_, status = st,
        stringsAsFactors = FALSE
      )
    }
    pdf_c <- bind_rows(p_rows)
    pdf_c$padj <- stats::p.adjust(pdf_c$p_value, method = "BH")
    pdf_c$direction <- ifelse(
      is.na(pdf_c$padj), NA_character_,
      ifelse(pdf_c$padj < 0.05 & pdf_c$mean_paired_delta > 0, "Up",
             ifelse(pdf_c$padj < 0.05 & pdf_c$mean_paired_delta < 0, "Down", "NS"))
    )
    paired_list[[cid]] <- pdf_c
  }

  # Stage stratification
  if (has_stage) {
    hcc_df <- s_clin[s_clin$TYPE == "HCC", , drop = FALSE]
    if (cid == "HCCDB13") {
      stg_grp <- map_bclc_stage(hcc_df$BCLC_STAGE)
      early_cols <- hcc_df$SAMPLE_ID[!is.na(stg_grp) & stg_grp == "Early_BCLC_0_A"]
      late_cols  <- hcc_df$SAMPLE_ID[!is.na(stg_grp) & stg_grp == "Late_BCLC_B_C"]
      mnote <- sprintf("BCLC 0/A -> Early (n=%d), BCLC B/C -> Late (n=%d); excluded NA/unmapped HCC=%d",
                       length(early_cols), length(late_cols), sum(is.na(stg_grp)))
      dom <- "Stage_BCLC"
    } else {
      raw_stg_col <- switch(
        cid,
        "HCCDB4"  = hcc_df$AJCC_STAGE,
        "HCCDB15" = hcc_df$AJCC_STAGE1,
        "HCCDB18" = hcc_df$TNM_STAGE_T1,
        "HCCDB24" = hcc_df$AJCC_Stage,
        "HCCDB25" = hcc_df$AJCC_Stage1,
        "HCCDB30" = hcc_df$AJCC_STAGE
      )
      stg_grp <- map_ajcc_stage(raw_stg_col)
      early_cols <- hcc_df$SAMPLE_ID[!is.na(stg_grp) & stg_grp == "Early_I_II"]
      late_cols  <- hcc_df$SAMPLE_ID[!is.na(stg_grp) & stg_grp == "Late_III_IV"]
      mnote <- sprintf("Stage I-II -> Early (n=%d), Stage III-IV -> Late (n=%d); excluded NA/unmapped HCC=%d",
                       length(early_cols), length(late_cols), sum(is.na(stg_grp)))
      dom <- "Stage_AJCC_TNM"
    }
    if (length(early_cols) >= 3L && length(adj_cols) >= 3L) {
      strat_list[[paste0(cid, "_Early_vs_Adj")]] <- run_unpaired_panel(
        mat, gene_matches, master_genes, early_cols, adj_cols, cid, dom,
        "Early", "Early_HCC_vs_Adjacent", "all_cohort_Adjacent", mnote
      )
    }
    if (length(late_cols) >= 3L && length(adj_cols) >= 3L) {
      strat_list[[paste0(cid, "_Late_vs_Adj")]] <- run_unpaired_panel(
        mat, gene_matches, master_genes, late_cols, adj_cols, cid, dom,
        "Late", "Late_HCC_vs_Adjacent", "all_cohort_Adjacent", mnote
      )
    }
    if (length(late_cols) >= 3L && length(early_cols) >= 3L) {
      strat_list[[paste0(cid, "_Late_vs_Early")]] <- run_unpaired_panel(
        mat, gene_matches, master_genes, late_cols, early_cols, cid, dom,
        "Late_vs_Early", "Late_HCC_vs_Early_HCC", "Early_HCC", mnote
      )
    }
  }

  # Etiology stratification
  if (has_etiol) {
    hcc_df <- s_clin[s_clin$TYPE == "HCC", , drop = FALSE]
    adj_df <- s_clin[s_clin$TYPE == "Adjacent", , drop = FALSE]
    etiol_subgroups <- list()
    if (cid == "HCCDB18") {
      etiol_subgroups[["HBV"]] <- list(hcc = hcc_df$SAMPLE_ID[hcc_df$VIRUS == "HBV" & !is.na(hcc_df$VIRUS)],
                                       adj = adj_df$SAMPLE_ID[adj_df$VIRUS == "HBV" & !is.na(adj_df$VIRUS)],
                                       note = "VIRUS=='HBV' (excluded 4 HBV,HCV co-infected)")
      etiol_subgroups[["HCV"]] <- list(hcc = hcc_df$SAMPLE_ID[hcc_df$VIRUS == "HCV" & !is.na(hcc_df$VIRUS)],
                                       adj = adj_df$SAMPLE_ID[adj_df$VIRUS == "HCV" & !is.na(adj_df$VIRUS)],
                                       note = "VIRUS=='HCV' (excluded 4 HBV,HCV co-infected)")
      etiol_subgroups[["Non_viral"]] <- list(hcc = hcc_df$SAMPLE_ID[hcc_df$VIRUS == "NBNC" & !is.na(hcc_df$VIRUS)],
                                             adj = adj_df$SAMPLE_ID[adj_df$VIRUS == "NBNC" & !is.na(adj_df$VIRUS)],
                                             note = "VIRUS=='NBNC' (Non-B Non-C)")
    } else if (cid == "HCCDB25") {
      etiol_subgroups[["HBV"]] <- list(hcc = hcc_cols, adj = adj_cols,
                                       note = "Known 100% HBV-related HCC cohort (Gao et al., 2019 CHCC-HBV)")
    } else if (cid == "HCCDB27") {
      etiol_subgroups[["HBV"]] <- list(hcc = hcc_df$SAMPLE_ID[hcc_df$HBV_STATE == "Positive" & !is.na(hcc_df$HBV_STATE)],
                                       adj = adj_df$SAMPLE_ID[adj_df$HBV_STATE == "Positive" & !is.na(adj_df$HBV_STATE)],
                                       note = "HBV_STATE=='Positive' (100% HBV cohort)")
    } else if (cid == "HCCDB30") {
      rf_h <- hcc_df[["Risk factor"]]; rf_a <- adj_df[["Risk factor"]]
      etiol_subgroups[["HBV"]] <- list(hcc = hcc_df$SAMPLE_ID[!is.na(rf_h) & rf_h == "HBV"],
                                       adj = adj_df$SAMPLE_ID[!is.na(rf_a) & rf_a == "HBV"],
                                       note = "Risk factor=='HBV'")
      etiol_subgroups[["HCV"]] <- list(hcc = hcc_df$SAMPLE_ID[!is.na(rf_h) & rf_h == "HCV"],
                                       adj = adj_df$SAMPLE_ID[!is.na(rf_a) & rf_a == "HCV"],
                                       note = "Risk factor=='HCV'")
      nv_vals <- c("nonBnonC", "Alcoholic", "NASH")
      etiol_subgroups[["Non_viral"]] <- list(hcc = hcc_df$SAMPLE_ID[!is.na(rf_h) & rf_h %in% nv_vals],
                                             adj = adj_df$SAMPLE_ID[!is.na(rf_a) & rf_a %in% nv_vals],
                                             note = "Risk factor in {nonBnonC, Alcoholic, NASH (n=1)}")
    }
    for (sg_name in names(etiol_subgroups)) {
      sg <- etiol_subgroups[[sg_name]]
      if (length(sg$hcc) < 3L) next
      if (length(sg$adj) >= 5L) {
        ref_cols <- sg$adj; ref_lbl <- "subgroup_paired_Adjacent"
      } else {
        ref_cols <- adj_cols; ref_lbl <- sprintf("all_cohort_Adjacent (subgroup_Adj n=%d < 5)", length(sg$adj))
      }
      strat_list[[paste0(cid, "_Etiol_", sg_name)]] <- run_unpaired_panel(
        mat, gene_matches, master_genes, sg$hcc, ref_cols, cid, "Etiology",
        sg_name, paste0(sg_name, "_HCC_vs_Adjacent"), ref_lbl, sg$note
      )
    }
  }
}

# ==============================================================================
# STEP 3: Export All 9 Statistical CSV Tables & sessionInfo.txt
# ==============================================================================
cat("[Step 3] Exporting all statistical CSV tables and sessionInfo.txt...\n")

cohort_check_df     <- bind_rows(cohort_check_list)
per_cohort_stats_df <- bind_rows(per_cohort_stats_list)
r41_per_cohort_df   <- bind_rows(r41_per_cohort_list)
audit_df            <- bind_rows(audit_list)
paired_per_cohort_df <- bind_rows(paired_list)
strat_df            <- bind_rows(strat_list)

n_included_all <- sum(cohort_check_df$included)
n_included_nontcga <- sum(cohort_check_df$included & cohort_check_df$cohort != "HCCDB15")

gene_recurrence_df <- bind_rows(lapply(seq_len(nrow(master_genes)), function(g) {
  disp <- master_genes$display_label[g]
  sub_all  <- per_cohort_stats_df[per_cohort_stats_df$display_label == disp, , drop = FALSE]
  sub_nt   <- sub_all[!sub_all$is_tcga, , drop = FALSE]
  sub_tcga <- sub_all[sub_all$is_tcga, , drop = FALSE]
  vd_nt <- sub_nt$delta_mean[!is.na(sub_nt$delta_mean)]
  data.frame(
    display_label = disp,
    hccdb_query_symbol = master_genes$hccdb_query_symbol[g],
    ptm_category = master_genes$ptm_category[g],
    tier = master_genes$tier[g],
    table2_protein_status = master_genes$table2_protein_status[g],
    n_detected_all = sum(sub_all$status != "not_detected"),
    n_tested_all   = sum(sub_all$status == "tested"),
    n_up_all       = sum(!is.na(sub_all$padj) & sub_all$padj < 0.05 & sub_all$delta_mean > 0),
    n_down_all     = sum(!is.na(sub_all$padj) & sub_all$padj < 0.05 & sub_all$delta_mean < 0),
    n_detected_nontcga = sum(sub_nt$status != "not_detected"),
    n_tested_nontcga   = sum(sub_nt$status == "tested"),
    n_up_nontcga       = sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$delta_mean > 0),
    n_down_nontcga     = sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$delta_mean < 0),
    n_up_nontcga_fc05   = sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$delta_mean >= 0.5),
    n_down_nontcga_fc05 = sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$delta_mean <= -0.5),
    median_delta_nontcga = if (length(vd_nt) > 0L) median(vd_nt) else NA_real_,
    tcga_delta_mean = if (nrow(sub_tcga) == 1L) sub_tcga$delta_mean[1] else NA_real_,
    tcga_padj       = if (nrow(sub_tcga) == 1L) sub_tcga$padj[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}))

r41_summary_df <- bind_rows(lapply(seq_len(nrow(master_genes)), function(g) {
  disp <- master_genes$display_label[g]
  sub_all  <- r41_per_cohort_df[r41_per_cohort_df$display_label == disp, , drop = FALSE]
  sub_nt   <- sub_all[!sub_all$is_tcga, , drop = FALSE]
  sub_tcga <- sub_all[sub_all$is_tcga, , drop = FALSE]
  data.frame(
    display_label = disp,
    hccdb_query_symbol = master_genes$hccdb_query_symbol[g],
    ptm_category = master_genes$ptm_category[g],
    tier = master_genes$tier[g],
    n_detected_all = sum(sub_all$status != "not_detected"),
    n_tested_all   = sum(sub_all$status == "tested"),
    n_pos_all      = sum(!is.na(sub_all$padj) & sub_all$padj < 0.05 & sub_all$rho > 0),
    n_neg_all      = sum(!is.na(sub_all$padj) & sub_all$padj < 0.05 & sub_all$rho < 0),
    median_rho_all = if (any(!is.na(sub_all$rho))) median(sub_all$rho, na.rm = TRUE) else NA_real_,
    n_detected_nontcga = sum(sub_nt$status != "not_detected"),
    n_tested_nontcga   = sum(sub_nt$status == "tested"),
    n_pos_nontcga      = sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$rho > 0),
    n_neg_nontcga      = sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$rho < 0),
    median_rho_nontcga = if (any(!is.na(sub_nt$rho))) median(sub_nt$rho, na.rm = TRUE) else NA_real_,
    tcga_rho  = if (nrow(sub_tcga) == 1L) sub_tcga$rho[1] else NA_real_,
    tcga_padj = if (nrow(sub_tcga) == 1L) sub_tcga$padj[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}))

paired_summary_df <- bind_rows(lapply(seq_len(nrow(master_genes)), function(g) {
  disp <- master_genes$display_label[g]
  sub_all <- paired_per_cohort_df[paired_per_cohort_df$display_label == disp, , drop = FALSE]
  sub_nt  <- sub_all[!sub_all$is_tcga, , drop = FALSE]
  data.frame(
    display_label = disp,
    hccdb_query_symbol = master_genes$hccdb_query_symbol[g],
    ptm_category = master_genes$ptm_category[g],
    tier = master_genes$tier[g],
    n_detected_paired_all = sum(sub_all$status != "not_detected"),
    n_tested_paired_all   = sum(sub_all$status == "tested"),
    n_up_paired_all       = sum(!is.na(sub_all$padj) & sub_all$padj < 0.05 & sub_all$mean_paired_delta > 0),
    n_down_paired_all     = sum(!is.na(sub_all$padj) & sub_all$padj < 0.05 & sub_all$mean_paired_delta < 0),
    median_paired_delta   = if (any(!is.na(sub_all$median_paired_delta))) median(sub_all$median_paired_delta, na.rm = TRUE) else NA_real_,
    median_mean_paired_delta = if (any(!is.na(sub_all$mean_paired_delta))) median(sub_all$mean_paired_delta, na.rm = TRUE) else NA_real_,
    n_detected_paired_nontcga = sum(sub_nt$status != "not_detected"),
    n_up_paired_nontcga       = sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$mean_paired_delta > 0),
    n_down_paired_nontcga     = sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$mean_paired_delta < 0),
    median_paired_delta_nontcga = if (any(!is.na(sub_nt$median_paired_delta))) median(sub_nt$median_paired_delta, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}))

write.csv(cohort_check_df,      file.path(res_dir, "cohort_check.csv"), row.names = FALSE)
write.csv(per_cohort_stats_df,  file.path(res_dir, "per_cohort_stats.csv"), row.names = FALSE)
write.csv(gene_recurrence_df,   file.path(res_dir, "gene_recurrence_summary.csv"), row.names = FALSE)
write.csv(r41_per_cohort_df,    file.path(res_dir, "supp_R4_1_per_cohort_are_corr.csv"), row.names = FALSE)
write.csv(r41_summary_df,       file.path(res_dir, "supp_R4_1_are_corr_summary.csv"), row.names = FALSE)
write.csv(audit_df,             file.path(res_dir, "supp_R4_5_clinical_metadata_audit.csv"), row.names = FALSE)
write.csv(paired_per_cohort_df, file.path(res_dir, "supp_R4_5_paired_per_cohort.csv"), row.names = FALSE)
write.csv(paired_summary_df,    file.path(res_dir, "supp_R4_5_paired_summary.csv"), row.names = FALSE)
write.csv(strat_df,             file.path(res_dir, "supp_R4_5_stratified_stats.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), con = file.path(res_dir, "sessionInfo.txt"))

cat("  Exported all 9 CSV tables + sessionInfo.txt to results/.\n\n")

# ==============================================================================
# STEP 4: Pure-R Composite Figure 5A-C (PDF + 300 dpi PNG)
# ==============================================================================
cat("[Step 4] Building composite Figure 5A-C (Figure5_complete.pdf & 300 dpi .png)...\n")

# --- PANEL A: 49-gene x 18-cohort ComplexHeatmap ---
included_cohorts_df <- cohort_check_df %>%
  filter(included) %>%
  mutate(cohort_num = get_cohort_num(cohort), is_rnaseq = cohort %in% rnaseq_cohorts) %>%
  arrange(desc(is_rnaseq), cohort_num)

ordered_cohorts <- included_cohorts_df$cohort
col_display_labels <- ifelse(
  ordered_cohorts == "HCCDB15", "HCCDB15 (TCGA-LIHC)",
  ifelse(ordered_cohorts == "HCCDB25", "HCCDB25 (Gao 2019)", ordered_cohorts)
)

gene_labels <- master_genes$display_label
delta_mat <- matrix(NA_real_, nrow = length(gene_labels), ncol = length(ordered_cohorts),
                    dimnames = list(gene_labels, ordered_cohorts))
padj_mat  <- matrix(NA_real_, nrow = length(gene_labels), ncol = length(ordered_cohorts),
                    dimnames = list(gene_labels, ordered_cohorts))
for (k in seq_len(nrow(per_cohort_stats_df))) {
  rg <- per_cohort_stats_df$display_label[k]
  cc <- per_cohort_stats_df$cohort[k]
  if (rg %in% gene_labels && cc %in% ordered_cohorts) {
    delta_mat[rg, cc] <- per_cohort_stats_df$delta_mean[k]
    padj_mat[rg, cc]  <- per_cohort_stats_df$padj[k]
  }
}
delta_mat_clamped <- pmax(pmin(delta_mat, 2), -2)
delta_mat_clamped[is.na(delta_mat)] <- NA_real_

ptm_levels <- unique(master_genes$ptm_category)
ptm_pretty_map <- c(
  "Target"                             = "Target",
  "Ubiquitination"                     = "Ubiquitination (E3)",
  "Ubiquitination_KEAP1_Sequestration" = "KEAP1 Sequestration",
  "Deubiquitination"                   = "Deubiquitination",
  "Phosphorylation"                    = "Phosphorylation",
  "Acetylation_Deacetylation"          = "Acetylation / Deacetylation",
  "SUMOylation"                        = "SUMOylation",
  "Glycosylation_Deglycation"          = "Glycosylation / Deglycation",
  "NRF2_Transcriptional_Coactivation"  = "Transcriptional Coactivation",
  "PARylation"                         = "PARylation",
  "Acetylation_Indirect"               = "Indirect (Acetylation)"
)
row_split_factor <- factor(ptm_pretty_map[master_genes$ptm_category],
                           levels = unname(ptm_pretty_map[ptm_levels]))

n_rnaseq <- sum(included_cohorts_df$is_rnaseq)
n_microarray <- sum(!included_cohorts_df$is_rnaseq)
rnaseq_grp_lbl <- sprintf("RNA-seq Cohorts (n=%d)", n_rnaseq)
microarray_grp_lbl <- sprintf("Microarray Cohorts (n=%d)", n_microarray)
col_split_factor <- factor(
  ifelse(included_cohorts_df$is_rnaseq, rnaseq_grp_lbl, microarray_grp_lbl),
  levels = c(rnaseq_grp_lbl, microarray_grp_lbl)
)

raw_prot <- trimws(master_genes$table2_protein_status)
prot_status_clean <- ifelse(
  grepl("^Up_Phospho", raw_prot), "Phospho Up",
  ifelse(grepl("^Up_", raw_prot), "Protein Up",
         ifelse(grepl("^Down_", raw_prot), "Protein Down", "NS / Not detected"))
)
prot_levels <- c("Protein Up", "Phospho Up", "Protein Down", "NS / Not detected")
prot_status_factor <- factor(prot_status_clean, levels = prot_levels)
prot_colors <- c(
  "Protein Up"        = "#D6604D",
  "Phospho Up"        = "#FDB863",
  "Protein Down"      = "#4393C3",
  "NS / Not detected" = "#EAEAEA"
)

bar_mat <- cbind(Up = gene_recurrence_df$n_up_all, Down = gene_recurrence_df$n_down_all)
rownames(bar_mat) <- gene_labels
bar_colors <- c("Up" = "#B2182B", "Down" = "#2166AC")

right_ha <- rowAnnotation(
  sig_cohorts = anno_barplot(
    bar_mat, beside = FALSE, gp = gpar(fill = bar_colors, col = NA),
    border = TRUE, bar_width = 0.72, ylim = c(0, n_included_all),
    axis_param = list(side = "bottom", at = c(0, 6, 12, 18), gp = gpar(fontsize = 7.8)),
    width = unit(27, "mm")
  ),
  published_proteomics = anno_simple(
    as.character(prot_status_factor), col = prot_colors,
    gp = gpar(col = "white", lwd = 0.5), width = unit(6.2, "mm")
  ),
  annotation_label = c(sprintf("Sig. Cohorts\n(n=%d)", n_included_all), "Published\nProteomics"),
  gap = unit(2.5, "mm"),
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 8.2, fontface = "bold"),
  show_legend = FALSE
)

na_fill_col <- "#D9D9D9"
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#FFFFFF", "#B2182B"))

ht_A <- Heatmap(
  delta_mat_clamped,
  name = "delta_mean",
  col = col_fun,
  na_col = na_fill_col,
  cluster_rows = FALSE, cluster_row_slices = FALSE,
  cluster_columns = FALSE, cluster_column_slices = FALSE,
  row_split = row_split_factor, column_split = col_split_factor,
  row_gap = unit(1.3, "mm"), column_gap = unit(3.0, "mm"),
  border = TRUE, rect_gp = gpar(col = "#EBEBEB", lwd = 0.45),
  row_labels = gene_labels, row_names_side = "left",
  row_names_gp = gpar(fontsize = 8.0, fontface = "italic"),
  row_title_side = "left", row_title_rot = 0,
  row_title_gp = gpar(fontsize = 8.2, fontface = "bold", col = "#222222"),
  column_labels = col_display_labels, column_names_side = "bottom",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8.0,
                         fontface = ifelse(ordered_cohorts %in% c("HCCDB15", "HCCDB25"), "bold", "plain")),
  column_title_gp = gpar(fontsize = 9.5, fontface = "bold"),
  right_annotation = right_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    qval <- padj_mat[i, j]; dval <- delta_mat_clamped[i, j]
    if (!is.na(qval) && qval < 0.05) {
      txt_col <- ifelse(!is.na(dval) && abs(dval) > 1.15, "white", "#111111")
      grid.text("*", x = x, y = y - unit(0.4, "mm"),
                gp = gpar(fontsize = 9.5, fontface = "bold", col = txt_col))
    }
  },
  width = ncol(delta_mat_clamped) * unit(6.9, "mm"),
  height = nrow(delta_mat_clamped) * unit(3.85, "mm"),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "mRNA delta_mean\n(HCC - Adjacent)",
    at = c(-2, -1, 0, 1, 2),
    labels = c("<= -2.0", "-1.0", "0", "+1.0", ">= +2.0"),
    title_gp = gpar(fontsize = 8.2, fontface = "bold"),
    labels_gp = gpar(fontsize = 7.8),
    legend_height = unit(32, "mm"),
    grid_width = unit(4.0, "mm"),
    border = "grey45"
  )
)

lgd_cell_key <- Legend(
  labels = c("padj < 0.05 (within-cohort BH)", "Not detected (NA)"),
  title = "Significance & Detection",
  title_gp = gpar(fontsize = 8.2, fontface = "bold"),
  labels_gp = gpar(fontsize = 7.8),
  grid_height = unit(3.8, "mm"), grid_width = unit(3.8, "mm"), row_gap = unit(0.9, "mm"),
  graphics = list(
    function(x, y, w, h) {
      grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey55"))
      grid.text("*", x, y - unit(0.35, "mm"), gp = gpar(fontsize = 10.5, fontface = "bold", col = "#111111"))
    },
    function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = na_fill_col, col = "grey55"))
  )
)

lgd_bar <- Legend(
  labels = c("Sig. Up (padj < 0.05, delta > 0)", "Sig. Down (padj < 0.05, delta < 0)"),
  title = sprintf("Sig. Cohorts (n=%d)", n_included_all),
  title_gp = gpar(fontsize = 8.2, fontface = "bold"),
  labels_gp = gpar(fontsize = 7.8),
  grid_height = unit(3.8, "mm"), grid_width = unit(3.8, "mm"), row_gap = unit(0.9, "mm"),
  legend_gp = gpar(fill = unname(bar_colors)), border = "grey55"
)

lgd_prot <- Legend(
  labels = prot_levels,
  title = "Published Proteomics\n/ Phosphoproteomics",
  title_gp = gpar(fontsize = 8.2, fontface = "bold"),
  labels_gp = gpar(fontsize = 7.8),
  grid_height = unit(3.8, "mm"), grid_width = unit(3.8, "mm"), row_gap = unit(0.9, "mm"),
  legend_gp = gpar(fill = unname(prot_colors[prot_levels])), border = "grey55"
)

draw_panel_A <- function() {
  ht_opt(HEATMAP_LEGEND_PADDING = unit(11, "mm"), ANNOTATION_LEGEND_PADDING = unit(11, "mm"))
  draw(
    ht_A,
    heatmap_legend_side = "right",
    annotation_legend_list = list(lgd_cell_key, lgd_bar, lgd_prot),
    annotation_legend_side = "right",
    merge_legend = TRUE,
    legend_gap = unit(4.8, "mm"),
    padding = unit(c(5, 8, 6, 6), "mm"),
    column_title = "Transcriptomic Dysregulation of 49 NRF2 PTM Regulators Across 18 HCCDB Cohorts (HCC vs. Adjacent)",
    column_title_gp = gpar(fontsize = 11.2, fontface = "bold")
  )
}

# --- PANEL B: 49-gene Tumor-Only Median Spearman rho vs 10-Gene ARE Score ---
df_B <- r41_summary_df %>%
  arrange(desc(median_rho_all)) %>%
  mutate(
    rank_idx = row_number(),
    net_sig_cohorts = n_pos_all - n_neg_all,
    sig_bucket = case_when(
      n_pos_all >= 10 ~ "Pos >= 10 cohorts",
      n_pos_all >= 4  ~ "Pos 4-9 cohorts",
      n_neg_all >= 4  ~ "Neg >= 4 cohorts",
      TRUE            ~ "Weak / Mixed (<4)"
    ),
    sig_bucket = factor(sig_bucket, levels = c(
      "Pos >= 10 cohorts", "Pos 4-9 cohorts", "Weak / Mixed (<4)", "Neg >= 4 cohorts"
    )),
    is_key = display_label %in% panelB_label_genes,
    lbl_text = ifelse(is_key, sprintf("%s (%+.2f)", display_label, median_rho_all), NA_character_)
  )

n_weak <- sum(abs(df_B$median_rho_all) < 0.20, na.rm = TRUE)

bucket_cols <- c(
  "Pos >= 10 cohorts" = "#B2182B",
  "Pos 4-9 cohorts"   = "#EF8A62",
  "Weak / Mixed (<4)" = "#8C8C8C",
  "Neg >= 4 cohorts"  = "#2166AC"
)

panel_B <- ggplot(df_B, aes(x = rank_idx, y = median_rho_all)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.20, ymax = 0.20,
           fill = "#E9EEF2", alpha = 0.85) +
  geom_hline(yintercept = 0, color = "#444444", linewidth = 0.5) +
  geom_hline(yintercept = c(-0.20, 0.20), linetype = "dashed", color = "#6C7A89", linewidth = 0.45) +
  annotate("text", x = 38.5, y = 0.24,
           label = sprintf("|rho| < 0.20 (%d/49 genes)", n_weak),
           size = 2.9, fontface = "bold", color = "#455A64", hjust = 0.5) +
  geom_segment(aes(xend = rank_idx, y = 0, yend = median_rho_all, color = sig_bucket),
               linewidth = 0.72, alpha = 0.9) +
  geom_point(aes(color = sig_bucket, size = is_key), alpha = 0.95) +
  scale_size_manual(values = c("FALSE" = 1.7, "TRUE" = 2.8), guide = "none") +
  scale_color_manual(values = bucket_cols, name = "Sig. Cohorts (padj < 0.05)") +
  geom_text_repel(
    data = df_B %>% filter(is_key, median_rho_all >= 0),
    aes(label = lbl_text, color = sig_bucket),
    size = 2.65, fontface = "bold.italic",
    seed = 2026,
    nudge_y = ifelse(df_B$median_rho_all[df_B$is_key & df_B$median_rho_all >= 0] > 0.5, 0.0, 0.14),
    nudge_x = ifelse(df_B$median_rho_all[df_B$is_key & df_B$median_rho_all >= 0] > 0.5, 5.0, 0.0),
    ylim = c(0.24, 0.74),
    box.padding = 0.35, point.padding = 0.18,
    min.segment.length = 0, segment.size = 0.32, segment.alpha = 0.7,
    max.overlaps = Inf, force = 6, force_pull = 0.4,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = df_B %>% filter(is_key, median_rho_all < 0),
    aes(label = lbl_text, color = sig_bucket),
    size = 2.65, fontface = "bold.italic",
    seed = 2026,
    nudge_y = -0.09,
    ylim = c(-0.31, -0.07),
    box.padding = 0.32, point.padding = 0.18,
    min.segment.length = 0, segment.size = 0.32, segment.alpha = 0.7,
    max.overlaps = Inf, force = 6, force_pull = 0.4,
    show.legend = FALSE
  ) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 49), expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(limits = c(-0.33, 0.76), breaks = seq(-0.2, 0.7, by = 0.2)) +
  labs(
    title = "Tumor-Only Correlation with 10-Gene NRF2 ARE Score",
    subtitle = "Median Spearman rho across 18 HCCDB cohorts (excluding SQSTM1 & NFE2L2 from score)",
    x = "49 NRF2 PTM Regulators (Ranked by Median Spearman rho)",
    y = "Median Spearman rho"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(size = 10.2, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7.8, color = "#444444", hjust = 0),
    axis.title = element_text(size = 8.6, face = "bold"),
    axis.text = element_text(size = 7.8, color = "#222222"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 7.8, face = "bold"),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(3.5, "mm"),
    legend.margin = margin(t = 0, b = 0),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 10)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.6)))

# --- PANEL C: 9-Gene Tripartite Sensitivity Heatmap (Paired | Stage | Etiology) ---
gene_group_map <- c(
  "SQSTM1" = "Recurrent\nUp",    "SAE1"   = "Recurrent\nUp",    "PARP1"  = "Recurrent\nUp",
  "KEAP1"  = "KEAP1–CUL3\n–RBX1","CUL3"   = "KEAP1–CUL3\n–RBX1","RBX1"   = "KEAP1–CUL3\n–RBX1",
  "NFE2L2" = "Target /\nDown",   "DCAF11" = "Target /\nDown",   "FYN"    = "Target /\nDown"
)
gene_levels_C <- rev(core9_genes)
group_levels_C <- c("Recurrent\nUp", "KEAP1–CUL3\n–RBX1", "Target /\nDown")

# 1) Left block: Paired (n=1,003 across 11 cohorts)
df_C_paired <- paired_summary_df %>%
  filter(display_label %in% core9_genes) %>%
  mutate(
    block = "Paired",
    col_id = "Paired (n=1,003)",
    val_raw = median_paired_delta,
    val_clamped = pmax(pmin(val_raw, 2), -2),
    dom_n = ifelse(median_paired_delta >= 0, n_up_paired_all, n_down_paired_all),
    cell_txt = sprintf("%d/%d", dom_n, n_detected_paired_all),
    txt_col = ifelse(!is.na(val_clamped) & abs(val_clamped) > 1.15, "white", "#111111"),
    is_paired = TRUE
  ) %>%
  select(display_label, block, col_id, val_clamped, cell_txt, txt_col, is_paired)

# 2) Middle block: Stage-Stratified (7 cohorts: RNA-seq DB15/18/25/30 left, Microarray DB4/13/24 right)
stage_cohort_order <- c("HCCDB15", "HCCDB18", "HCCDB25", "HCCDB30", "HCCDB4", "HCCDB13", "HCCDB24")
short_cid <- function(x) sub("^HCCDB", "DB", x)

stage_col_levels <- unlist(lapply(stage_cohort_order, function(c0) {
  sc <- ifelse(c0 == "HCCDB13", "DB13*", short_cid(c0))
  c(paste0(sc, "-E"), paste0(sc, "-L"))
}))

df_C_stage <- strat_df %>%
  filter(
    display_label %in% core9_genes,
    strat_domain %in% c("Stage_AJCC_TNM", "Stage_BCLC"),
    comparison %in% c("Early_HCC_vs_Adjacent", "Late_HCC_vs_Adjacent")
  ) %>%
  mutate(
    sc = ifelse(cohort == "HCCDB13", "DB13*", short_cid(cohort)),
    suf = ifelse(subgroup == "Early", "-E", "-L"),
    block = "Stage-Stratified (7 Cohorts: RNA-seq | Microarray)",
    col_id = paste0(sc, suf),
    val_clamped = pmax(pmin(delta_mean, 2), -2),
    cell_txt = ifelse(!is.na(padj) & padj < 0.05, "*", ""),
    txt_col = ifelse(!is.na(val_clamped) & abs(val_clamped) > 1.15, "white", "#111111"),
    is_paired = FALSE
  ) %>%
  select(display_label, block, col_id, val_clamped, cell_txt, txt_col, is_paired)

# 3) Right block: Etiology-Stratified (4 cohorts)
etiol_col_levels <- c(
  "DB18-HBV", "DB18-HCV", "DB18-NV",
  "DB25-HBV", "DB27-HBV",
  "DB30-HBV", "DB30-HCV", "DB30-NV"
)

df_C_etiol <- strat_df %>%
  filter(display_label %in% core9_genes, strat_domain == "Etiology") %>%
  mutate(
    sg_short = ifelse(subgroup == "Non_viral", "NV", subgroup),
    block = "Etiology-Stratified (4 Cohorts)",
    col_id = paste0(short_cid(cohort), "-", sg_short),
    val_clamped = pmax(pmin(delta_mean, 2), -2),
    cell_txt = ifelse(!is.na(padj) & padj < 0.05, "*", ""),
    txt_col = ifelse(!is.na(val_clamped) & abs(val_clamped) > 1.15, "white", "#111111"),
    is_paired = FALSE
  ) %>%
  select(display_label, block, col_id, val_clamped, cell_txt, txt_col, is_paired)

block_levels_C <- c(
  "Paired",
  "Stage-Stratified (7 Cohorts: RNA-seq | Microarray)",
  "Etiology-Stratified (4 Cohorts)"
)

df_C_all <- bind_rows(df_C_paired, df_C_stage, df_C_etiol) %>%
  mutate(
    display_label = factor(display_label, levels = gene_levels_C),
    gene_group = factor(gene_group_map[as.character(display_label)], levels = group_levels_C),
    block = factor(block, levels = block_levels_C),
    col_id = factor(col_id, levels = c("Paired (n=1,003)", stage_col_levels, etiol_col_levels))
  )

# Vertical black divider between RNA-seq (DB15/18/25/30 = 8 cols) and Microarray (DB4/13*/24) in Stage facet
vline_stage_df <- data.frame(
  block = factor("Stage-Stratified (7 Cohorts: RNA-seq | Microarray)", levels = block_levels_C),
  xint = 8.5
)

panel_C <- ggplot(df_C_all, aes(x = col_id, y = display_label, fill = val_clamped)) +
  geom_tile(color = "#E5E5E5", linewidth = 0.4) +
  geom_text(
    aes(label = cell_txt, color = txt_col, size = is_paired),
    fontface = "bold", vjust = 0.55, show.legend = FALSE
  ) +
  geom_vline(data = vline_stage_df, aes(xintercept = xint), color = "#111111", linewidth = 0.9) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
    midpoint = 0, limits = c(-2, 2), na.value = na_fill_col,
    guide = "none"
  ) +
  scale_color_identity() +
  scale_size_manual(values = c("TRUE" = 2.6, "FALSE" = 3.8)) +
  facet_grid(gene_group ~ block, scales = "free", space = "free", switch = "y") +
  labs(
    title = "Sensitivity Analyses of 9 Focal Genes: Patient-Paired, Stage & Etiology Stratification",
    subtitle = "Shared mRNA delta_mean color scale [-2, +2] with Panel A; Paired cells show sig./detected cohorts; * padj < 0.05 (BH)",
    x = NULL, y = NULL,
    caption = "E = Early (Stage I-II), L = Late (Stage III-IV); DB13* uses BCLC (0/A vs B/C); NV = Non-viral. MASLD/NASH unannotated in 17/18 cohorts (n = 1 in DB30)."
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(size = 10.2, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7.7, color = "#444444", hjust = 0),
    plot.caption = element_text(size = 7.3, color = "#333333", hjust = 0, face = "italic", margin = margin(t = 4)),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 7.4, color = "#222222"),
    axis.text.y = element_text(face = "italic", size = 8.2, color = "#111111"),
    strip.background = element_rect(fill = "#F3F5F7", color = "#B0BEC5", linewidth = 0.5),
    strip.text.x = element_text(size = 7.8, face = "bold", color = "#111111", margin = margin(t = 3, b = 3)),
    strip.text.y.left = element_text(angle = 0, size = 7.6, face = "bold", color = "#222222"),
    strip.placement = "outside",
    panel.spacing.x = unit(2.6, "mm"),
    panel.spacing.y = unit(1.6, "mm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 8, r = 10, b = 8, l = 6)
  )

# --- Assemble Composite Figure 5A-C ---
render_complete_figure5 <- function() {
  grob_A <- grid.grabExpr(draw_panel_A(), width = 14.4, height = 10.6)
  bottom_row <- plot_grid(
    panel_B, panel_C,
    ncol = 2, rel_widths = c(0.40, 0.60),
    labels = c("B", "C"),
    label_size = 15, label_fontface = "bold",
    label_x = 0.008, label_y = 0.995,
    align = "none"
  )
  full_fig <- plot_grid(
    grob_A, bottom_row,
    ncol = 1, rel_heights = c(0.65, 0.35),
    labels = c("A", ""),
    label_size = 15, label_fontface = "bold",
    label_x = 0.008, label_y = 0.995
  )
  print(full_fig)
}

fig_w <- 14.4
fig_h <- 16.4
pdf_out <- file.path(res_dir, "Figure5_complete.pdf")
png_out <- file.path(res_dir, "Figure5_complete.png")

grDevices::pdf(pdf_out, width = fig_w, height = fig_h, useDingbats = FALSE)
render_complete_figure5()
grDevices::dev.off()

grDevices::png(png_out, width = fig_w, height = fig_h, units = "in", res = 300)
render_complete_figure5()
grDevices::dev.off()

cat("  Saved:", pdf_out, "\n")
cat("  Saved:", png_out, "(300 dpi)\n")
cat("\nUnified HCCDB analysis & Figure 5A-C pipeline completed successfully!\n")
