#!/usr/bin/env Rscript
# ==============================================================================
# run_hccdb_revision.R
# Independent within-cohort transcriptomic analysis of 49 NRF2 PTM regulators
# across HCCDB cohorts (Tumor vs. Adjacent) and generation of Figure 5.
# ==============================================================================

# 0. Configure library path & check required packages (DO NOT auto-install)
custom_lib <- "D:/ruanjian/R-4.6.1/my_library"
if (dir.exists(custom_lib)) {
  .libPaths(c(custom_lib, .libPaths()))
}

required_pkgs <- c("data.table", "dplyr", "ComplexHeatmap", "circlize", "grid")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    sprintf(
      "Missing required R package(s) in %s: %s. Please confirm before installing.",
      custom_lib, paste(missing_pkgs, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

cat("======================================================================\n")
cat("NRF2 Review Revision - HCCDB Multi-Cohort Independent Re-analysis\n")
cat("R version:", R.version.string, "\n")
cat("Library paths:", paste(.libPaths(), collapse = " | "), "\n")
cat("======================================================================\n\n")

# ------------------------------------------------------------------------------
# Paths & Directories (HCCDB source directory is strictly READ-ONLY)
# ------------------------------------------------------------------------------
raw_hccdb_dir <- "D:/project_TACE-PGCC-cell/04_Bioinformatics/PUBLIC_DATA/HCCDB"
gene_list_file <- "master_gene_list.csv"
tmp_dir <- "tmp_hccdb"
res_dir <- "results"

if (!dir.exists(raw_hccdb_dir)) {
  stop("Read-only HCCDB directory not found: ", raw_hccdb_dir, call. = FALSE)
}
if (!file.exists(gene_list_file)) {
  stop("Gene list file not found in working directory: ", gene_list_file, call. = FALSE)
}

dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

# Extract numeric cohort index for natural ordering (e.g., HCCDB1, HCCDB3, ..., HCCDB30)
get_cohort_num <- function(cid) {
  as.integer(sub("^HCCDB", "", cid))
}

# Parse HCCDB sample metadata table (handles transposed row-attribute format with
# potentially ragged trailing metadata rows, as well as standard column format)
parse_sample_table <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  split_lines <- lapply(lines, function(x) trimws(strsplit(x, "\t", fixed = TRUE)[[1]]))
  first_tokens <- vapply(split_lines, function(x) if (length(x) > 0) x[1] else "", character(1))

  sid_row <- which(toupper(first_tokens) == "SAMPLE_ID")
  type_row <- which(toupper(first_tokens) == "TYPE")

  if (length(sid_row) >= 1 && length(type_row) >= 1) {
    sids <- split_lines[[sid_row[1]]][-1]
    types <- split_lines[[type_row[1]]][-1]
    n <- min(length(sids), length(types))
    df <- data.frame(
      SAMPLE_ID = sids[seq_len(n)],
      TYPE = types[seq_len(n)],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    df <- df[nzchar(df$SAMPLE_ID), , drop = FALSE]
    return(df)
  }

  # Fallback: standard column-wise sample table
  dt <- fread(filepath, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
  cnames_up <- toupper(trimws(colnames(dt)))
  sid_col <- which(cnames_up == "SAMPLE_ID")[1]
  type_col <- which(cnames_up == "TYPE")[1]
  if (is.na(sid_col) || is.na(type_col)) {
    stop("Could not locate SAMPLE_ID and TYPE in sample file: ", filepath, call. = FALSE)
  }
  df <- data.frame(
    SAMPLE_ID = trimws(as.character(dt[[sid_col]])),
    TYPE = trimws(as.character(dt[[type_col]])),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  df <- df[nzchar(df$SAMPLE_ID), , drop = FALSE]
  df
}

# Format raw sample type counts concisely
format_type_counts <- function(type_vec) {
  type_clean <- ifelse(is.na(type_vec) | !nzchar(type_vec), "<NA>", type_vec)
  tb <- table(type_clean)
  # Put HCC and Adjacent first if present, then alphabetical
  ord_names <- c(
    intersect(c("HCC", "Adjacent"), names(tb)),
    sort(setdiff(names(tb), c("HCC", "Adjacent")))
  )
  tb <- tb[ord_names]
  paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = "; ")
}

# Read mRNA_level3 table with check.names = FALSE and identify Symbol + sample columns
read_mrna_matrix <- function(filepath) {
  dt <- fread(
    filepath,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    data.table = FALSE,
    na.strings = c("NA", "NaN", "", "null", "NULL")
  )
  cnames <- colnames(dt)
  cnames_low <- tolower(trimws(cnames))

  # Identify leading annotation columns (commonly Entrez_ID + Symbol, or Symbol alone)
  sym_col <- which(cnames_low %in% c("symbol", "gene_symbol", "genesymbol", "gene"))[1]
  entrez_col <- which(grepl("^entrez", cnames_low))[1]

  if (is.na(sym_col)) {
    # Fallback heuristic: if col 1 is numeric Entrez ID, col 2 is Symbol; else col 1
    if (ncol(dt) >= 2 && is.numeric(dt[[1]]) && is.character(dt[[2]])) {
      entrez_col <- 1L
      sym_col <- 2L
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

  list(
    symbols = symbols,
    mat = mat,
    meta_col_names = cnames[meta_cols]
  )
}

# ------------------------------------------------------------------------------
# Step 1. Reuse existing tmp_hccdb/ if non-empty; otherwise unzip from read-only dir
# ------------------------------------------------------------------------------
existing_tmp_files <- if (dir.exists(tmp_dir)) {
  list.files(tmp_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
} else {
  character(0)
}

if (length(existing_tmp_files) > 0L) {
  cat(sprintf(
    "[Step 1] Reusing %d existing extracted .txt files in '%s/' (skipping unzip; read-only source untouched).\n",
    length(existing_tmp_files), tmp_dir
  ))
  sample_files_vec <- existing_tmp_files[grepl("\\.sample\\.txt$", basename(existing_tmp_files))]
  sample_tbl <- data.frame(
    cohort = sub("\\.sample\\.txt$", "", basename(sample_files_vec)),
    sample_path = sample_files_vec,
    stringsAsFactors = FALSE
  )

  # Identify mRNA_level3 files (all .txt files excluding .sample.txt and .patient.txt)
  mrna_candidates <- existing_tmp_files[
    !grepl("\\.(sample|patient)\\.txt$", basename(existing_tmp_files))
  ]
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
  })) %>%
    filter(!is.na(cohort))

  cohort_files <- full_join(sample_tbl, mrna_tbl, by = "cohort") %>%
    mutate(cohort_num = get_cohort_num(cohort)) %>%
    arrange(cohort_num)
  cat(sprintf("  Mapped %d cohorts directly from '%s/'.\n\n", nrow(cohort_files), tmp_dir))
} else {
  cat("[Step 1] Extracting all .zip archives from read-only HCCDB folder to tmp_hccdb/...\n")
  zip_files <- sort(list.files(raw_hccdb_dir, pattern = "\\.zip$", full.names = TRUE))
  if (length(zip_files) == 0) {
    stop("No .zip files found in: ", raw_hccdb_dir, call. = FALSE)
  }

  zip_manifest <- lapply(zip_files, function(zf) {
    toc <- unzip(zf, list = TRUE)
    unzip(zf, exdir = tmp_dir, overwrite = TRUE)
    inner_files <- toc$Name[!grepl("/$", toc$Name) & toc$Length > 0]
    data.frame(
      zip_name = basename(zf),
      extracted_path = file.path(tmp_dir, inner_files[1]),
      stringsAsFactors = FALSE
    )
  })
  zip_manifest_df <- bind_rows(zip_manifest)
  cat(sprintf("  Extracted %d zip archives into '%s/'.\n\n", nrow(zip_manifest_df), tmp_dir))

  sample_zips <- zip_manifest_df %>%
    filter(grepl("\\.sample\\.zip$", zip_name)) %>%
    mutate(cohort = sub("\\.sample\\.zip$", "", zip_name))

  mrna_zips <- zip_manifest_df %>%
    filter(grepl("_mRNA_level3\\.zip$", zip_name)) %>%
    mutate(cohort = sub("_mRNA_level3\\.zip$", "", zip_name))

  cohort_files <- full_join(
    sample_zips %>% select(cohort, sample_path = extracted_path),
    mrna_zips %>% select(cohort, mrna_path = extracted_path),
    by = "cohort"
  ) %>%
    mutate(cohort_num = get_cohort_num(cohort)) %>%
    arrange(cohort_num)
}

# Known RNA-seq cohorts in HCCDB
rnaseq_cohorts <- c("HCCDB15", "HCCDB18", "HCCDB25", "HCCDB30")

# Read master_gene_list.csv (49 genes)
master_genes <- read.csv(
  gene_list_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
master_genes$display_label <- trimws(master_genes$display_label)
master_genes$hccdb_query_symbol <- trimws(master_genes$hccdb_query_symbol)
cat(sprintf("[Step 2] Loaded %d target genes from '%s'.\n\n", nrow(master_genes), gene_list_file))

# Legacy alias fallback map (official symbol tried first; fallbacks tried in order)
legacy_alias_map <- list(
  "OGA"     = c("MGEA5"),
  "KAT8"    = c("MYST1", "MOF"),
  "NAA10"   = c("ARD1A", "ARD1"),
  "DCAF11"  = c("WDR23"),
  "SKP1"    = c("SKP1A"),
  "USP17L2" = c("DUB3"),
  "SYVN1"   = c("HRD1")
)

# ------------------------------------------------------------------------------
# Step 1b & Step 2. Screen each cohort independently & run within-cohort stats
# ------------------------------------------------------------------------------
cohort_check_list <- list()
per_cohort_stats_list <- list()

cat("[Step 2] Screening cohorts and running independent within-cohort statistics...\n")

for (i in seq_len(nrow(cohort_files))) {
  cid <- cohort_files$cohort[i]
  spath <- cohort_files$sample_path[i]
  mpath <- cohort_files$mrna_path[i]
  platform <- ifelse(cid %in% rnaseq_cohorts, "RNA-seq", "Microarray")

  if (is.na(spath) || !file.exists(spath) || is.na(mpath) || !file.exists(mpath)) {
    cohort_check_list[[cid]] <- data.frame(
      cohort = cid,
      platform = platform,
      sample_file = ifelse(is.na(spath), NA_character_, basename(spath)),
      mrna_file = ifelse(is.na(mpath), NA_character_, basename(mpath)),
      raw_type_counts = NA_character_,
      n_mrna_samples_total = NA_integer_,
      n_HCC = 0L,
      n_Adjacent = 0L,
      expr_min = NA_real_,
      expr_median = NA_real_,
      expr_max = NA_real_,
      log2_transformed = FALSE,
      included = FALSE,
      exclusion_reason = "Missing sample or mRNA_level3 file",
      stringsAsFactors = FALSE
    )
    next
  }

  s_df <- parse_sample_table(spath)
  raw_counts_str <- format_type_counts(s_df$TYPE)

  mrna_obj <- read_mrna_matrix(mpath)
  mat <- mrna_obj$mat
  syms <- mrna_obj$symbols
  mrna_sids <- colnames(mat)

  # STRICT group matching: ONLY exact "HCC" and exact "Adjacent" labels are allowed
  sids_hcc_raw <- s_df$SAMPLE_ID[s_df$TYPE == "HCC"]
  sids_adj_raw <- s_df$SAMPLE_ID[s_df$TYPE == "Adjacent"]

  hcc_cols <- mrna_sids[mrna_sids %in% sids_hcc_raw]
  adj_cols <- mrna_sids[mrna_sids %in% sids_adj_raw]

  n_hcc <- length(hcc_cols)
  n_adj <- length(adj_cols)

  # Expression distribution check on raw matrix
  expr_min_raw <- suppressWarnings(min(mat, na.rm = TRUE))
  expr_med_raw <- suppressWarnings(median(mat, na.rm = TRUE))
  expr_max_raw <- suppressWarnings(max(mat, na.rm = TRUE))

  # Expression scale transformation:
  # 1) Explicit log2 transformation for linear-scale HCCDB3 (GSE25097):
  #    log2(x) if strictly positive (min > 0), otherwise log2(x + 1)
  # 2) General rule for other cohorts: log2(x + 1) when max > 100 and min >= 0
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
    cat(sprintf("  [%s] %s (raw min=%.4f, max=%.3f -> post min=%.3f, max=%.3f).\n",
                cid, scale_note, expr_min_raw, expr_max_raw,
                min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)))
  } else if (is.finite(expr_max_raw) && is.finite(expr_min_raw) && expr_max_raw > 100 && expr_min_raw >= 0) {
    mat <- log2(mat + 1)
    did_log2 <- TRUE
    scale_note <- "Auto log2(x+1) applied (max > 100 & min >= 0)"
    cat(sprintf("  [%s] Applied log2(x+1) transformation (raw min=%.3f, max=%.3f).\n",
                cid, expr_min_raw, expr_max_raw))
  }

  # Inclusion criteria: both HCC and Adjacent present with >= 3 matched samples each
  included <- (n_hcc >= 3L && n_adj >= 3L)
  if (included) {
    other_types <- setdiff(unique(s_df$TYPE), c("HCC", "Adjacent"))
    notes_vec <- character(0)
    if (length(other_types) > 0) {
      notes_vec <- c(notes_vec, sprintf("ignored non-target sample types: %s", paste(other_types, collapse = ", ")))
    }
    if (nzchar(scale_note)) {
      notes_vec <- c(notes_vec, scale_note)
    }
    excl_reason <- if (length(notes_vec) > 0) {
      sprintf("Included (%s)", paste(notes_vec, collapse = "; "))
    } else {
      "Included"
    }
  } else {
    if (n_hcc < 3L && n_adj < 3L) {
      excl_reason <- sprintf("Insufficient HCC (%d) and Adjacent (%d) samples (<3)", n_hcc, n_adj)
    } else if (n_adj == 0L) {
      non_hcc <- setdiff(unique(s_df$TYPE), "HCC")
      if (length(non_hcc) > 0) {
        excl_reason <- sprintf("No Adjacent samples (has non-Adjacent control label(s): %s; not mixed)",
                               paste(non_hcc, collapse = ", "))
      } else {
        excl_reason <- "No Adjacent samples (HCC tumor-only cohort)"
      }
    } else if (n_adj < 3L) {
      excl_reason <- sprintf("Insufficient Adjacent samples (n_Adjacent=%d < 3)", n_adj)
    } else {
      excl_reason <- sprintf("Insufficient HCC samples (n_HCC=%d < 3)", n_hcc)
    }
  }

  cohort_check_list[[cid]] <- data.frame(
    cohort = cid,
    platform = platform,
    sample_file = basename(spath),
    mrna_file = basename(mpath),
    raw_type_counts = raw_counts_str,
    n_mrna_samples_total = ncol(mat),
    n_HCC = n_hcc,
    n_Adjacent = n_adj,
    expr_min = round(expr_min_raw, 6),
    expr_median = round(expr_med_raw, 6),
    expr_max = round(expr_max_raw, 6),
    log2_transformed = did_log2,
    included = included,
    exclusion_reason = excl_reason,
    stringsAsFactors = FALSE
  )

  cat(sprintf(
    "  %-7s | %-10s | raw=[%s] | matched HCC=%3d, Adj=%3d | min=%6.2f med=%5.2f max=%6.2f | log2=%s | %s\n",
    cid, platform, raw_counts_str, n_hcc, n_adj,
    expr_min_raw, expr_med_raw, expr_max_raw,
    ifelse(did_log2, "YES", "NO "),
    ifelse(included, "INCLUDED", paste0("EXCLUDED: ", excl_reason))
  ))

  if (!included) {
    next
  }

  # ----------------------------------------------------------------------------
  # Within-cohort deduplication & statistical testing for the 49 target genes
  # ----------------------------------------------------------------------------
  # Also log matrix-wide duplicate symbol summary if any exist
  valid_sym_mask <- !is.na(syms) & nzchar(syms)
  n_dup_matrix_syms <- sum(duplicated(syms[valid_sym_mask]))
  all_row_means <- rowMeans(mat, na.rm = TRUE)

  gene_rows <- vector("list", nrow(master_genes))
  dup_target_logs <- character(0)

  for (g in seq_len(nrow(master_genes))) {
    disp_sym <- master_genes$display_label[g]
    query_sym <- master_genes$hccdb_query_symbol[g]
    alias_syms <- if (disp_sym %in% names(legacy_alias_map)) legacy_alias_map[[disp_sym]] else character(0)
    # Priority order: official symbol (display_label) first, then hccdb_query_symbol & legacy aliases in order
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
      # Missing gene: strictly retain NA (never impute 0)
      gene_rows[[g]] <- data.frame(
        cohort = cid,
        is_tcga = (cid == "HCCDB15"),
        platform = platform,
        display_label = disp_sym,
        hccdb_query_symbol = query_sym,
        matched_symbol = NA_character_,
        ptm_category = master_genes$ptm_category[g],
        tier = master_genes$tier[g],
        n_HCC_cohort = n_hcc,
        n_Adjacent_cohort = n_adj,
        n_HCC_nonNA = 0L,
        n_Adjacent_nonNA = 0L,
        mean_HCC = NA_real_,
        mean_Adjacent = NA_real_,
        delta_mean = NA_real_,
        median_HCC = NA_real_,
        median_Adjacent = NA_real_,
        delta_median = NA_real_,
        p_value = NA_real_,
        padj = NA_real_,
        status = "not_detected",
        stringsAsFactors = FALSE
      )
      next
    }

    if (length(hit_idx) > 1L) {
      means_cand <- all_row_means[hit_idx]
      # Handle edge case where all candidate rows are all-NA
      if (all(is.na(means_cand))) {
        best_pos <- 1L
      } else {
        best_pos <- which.max(means_cand)
      }
      chosen_idx <- hit_idx[best_pos]
      dup_target_logs <- c(
        dup_target_logs,
        sprintf(
          "    [Dedup %s] %s (matched: %s, %d rows) -> kept row %d (Symbol='%s', all-sample mean=%.4f; dropped mean(s)=%s)",
          cid, disp_sym, paste(unique(syms[hit_idx]), collapse = "/"), length(hit_idx),
          chosen_idx, syms[chosen_idx], all_row_means[chosen_idx],
          paste(sprintf("%.4f", means_cand[-best_pos]), collapse = ", ")
        )
      )
    } else {
      chosen_idx <- hit_idx[1L]
    }

    matched_sym <- syms[chosen_idx]
    x_hcc <- as.numeric(mat[chosen_idx, hcc_cols])
    x_adj <- as.numeric(mat[chosen_idx, adj_cols])

    x_hcc_ok <- x_hcc[!is.na(x_hcc)]
    x_adj_ok <- x_adj[!is.na(x_adj)]

    n_hcc_ok <- length(x_hcc_ok)
    n_adj_ok <- length(x_adj_ok)

    mean_hcc <- if (n_hcc_ok > 0L) mean(x_hcc_ok) else NA_real_
    mean_adj <- if (n_adj_ok > 0L) mean(x_adj_ok) else NA_real_
    med_hcc  <- if (n_hcc_ok > 0L) median(x_hcc_ok) else NA_real_
    med_adj  <- if (n_adj_ok > 0L) median(x_adj_ok) else NA_real_

    d_mean <- if (n_hcc_ok > 0L && n_adj_ok > 0L) (mean_hcc - mean_adj) else NA_real_
    d_med  <- if (n_hcc_ok > 0L && n_adj_ok > 0L) (med_hcc - med_adj) else NA_real_

    if (n_hcc_ok < 3L || n_adj_ok < 3L) {
      st <- "insufficient_nonNA"
      pval <- NA_real_
    } else {
      comb_vals <- c(x_hcc_ok, x_adj_ok)
      if (length(unique(comb_vals)) <= 1L || stats::var(comb_vals) == 0) {
        st <- "zero_variance"
        pval <- NA_real_
      } else {
        wt <- suppressWarnings(
          stats::wilcox.test(x_hcc_ok, x_adj_ok, alternative = "two.sided", exact = FALSE)
        )
        pval <- as.numeric(wt$p.value)
        st <- "tested"
      }
    }

    gene_rows[[g]] <- data.frame(
      cohort = cid,
      is_tcga = (cid == "HCCDB15"),
      platform = platform,
      display_label = disp_sym,
      hccdb_query_symbol = query_sym,
      matched_symbol = matched_sym,
      ptm_category = master_genes$ptm_category[g],
      tier = master_genes$tier[g],
      n_HCC_cohort = n_hcc,
      n_Adjacent_cohort = n_adj,
      n_HCC_nonNA = n_hcc_ok,
      n_Adjacent_nonNA = n_adj_ok,
      mean_HCC = mean_hcc,
      mean_Adjacent = mean_adj,
      delta_mean = d_mean,
      median_HCC = med_hcc,
      median_Adjacent = med_adj,
      delta_median = d_med,
      p_value = pval,
      padj = NA_real_,
      status = st,
      stringsAsFactors = FALSE
    )
  }

  if (length(dup_target_logs) > 0L) {
    cat(paste(dup_target_logs, collapse = "\n"), "\n")
  } else {
    cat(sprintf(
      "    [Dedup %s] 0 duplicate rows among detected target genes (matrix-wide duplicate symbols: %d).\n",
      cid, n_dup_matrix_syms
    ))
  }

  cohort_df <- bind_rows(gene_rows)
  # Within-cohort Benjamini-Hochberg adjustment across tested target genes
  cohort_df$padj <- stats::p.adjust(cohort_df$p_value, method = "BH")
  cohort_df$direction <- ifelse(
    is.na(cohort_df$padj),
    NA_character_,
    ifelse(
      cohort_df$padj < 0.05 & cohort_df$delta_mean > 0, "Up",
      ifelse(cohort_df$padj < 0.05 & cohort_df$delta_mean < 0, "Down", "NS")
    )
  )

  per_cohort_stats_list[[cid]] <- cohort_df
}

# ------------------------------------------------------------------------------
# Step 3. Export summary tables & sessionInfo
# ------------------------------------------------------------------------------
cat("\n[Step 3] Exporting summary CSV tables and sessionInfo.txt...\n")

cohort_check_df <- bind_rows(cohort_check_list)
write.csv(cohort_check_df, file.path(res_dir, "cohort_check.csv"), row.names = FALSE)

per_cohort_stats_df <- bind_rows(per_cohort_stats_list)
write.csv(per_cohort_stats_df, file.path(res_dir, "per_cohort_stats.csv"), row.names = FALSE)

n_included_all <- sum(cohort_check_df$included)
n_included_nontcga <- sum(cohort_check_df$included & cohort_check_df$cohort != "HCCDB15")
cat(sprintf("  Included cohorts: %d total (%d independent non-TCGA + %d TCGA [HCCDB15]).\n",
            n_included_all, n_included_nontcga, n_included_all - n_included_nontcga))

# Summarize recurrence per gene (strictly preserving master_gene_list.csv order)
summary_rows <- lapply(seq_len(nrow(master_genes)), function(g) {
  disp_sym <- master_genes$display_label[g]
  sub_all <- per_cohort_stats_df[per_cohort_stats_df$display_label == disp_sym, , drop = FALSE]
  sub_nt  <- sub_all[!sub_all$is_tcga, , drop = FALSE]
  sub_tcga <- sub_all[sub_all$is_tcga, , drop = FALSE]

  # Detected = gene exists in the cohort matrix with >=1 non-NA observation (status != "not_detected")
  n_det_all  <- sum(sub_all$status != "not_detected")
  n_test_all <- sum(sub_all$status == "tested")
  n_up_all   <- sum(!is.na(sub_all$padj) & sub_all$padj < 0.05 & sub_all$delta_mean > 0)
  n_down_all <- sum(!is.na(sub_all$padj) & sub_all$padj < 0.05 & sub_all$delta_mean < 0)

  n_det_nt  <- sum(sub_nt$status != "not_detected")
  n_test_nt <- sum(sub_nt$status == "tested")
  n_up_nt   <- sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$delta_mean > 0)
  n_down_nt <- sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$delta_mean < 0)

  n_up_nt_fc05   <- sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$delta_mean >= 0.5)
  n_down_nt_fc05 <- sum(!is.na(sub_nt$padj) & sub_nt$padj < 0.05 & sub_nt$delta_mean <= -0.5)

  valid_deltas_nt <- sub_nt$delta_mean[!is.na(sub_nt$delta_mean)]
  med_delta_nt <- if (length(valid_deltas_nt) > 0L) median(valid_deltas_nt) else NA_real_

  tcga_d <- if (nrow(sub_tcga) == 1L) sub_tcga$delta_mean[1] else NA_real_
  tcga_q <- if (nrow(sub_tcga) == 1L) sub_tcga$padj[1] else NA_real_

  data.frame(
    display_label = disp_sym,
    hccdb_query_symbol = master_genes$hccdb_query_symbol[g],
    ptm_category = master_genes$ptm_category[g],
    tier = master_genes$tier[g],
    table2_protein_status = master_genes$table2_protein_status[g],
    n_detected_all = n_det_all,
    n_tested_all = n_test_all,
    n_up_all = n_up_all,
    n_down_all = n_down_all,
    n_detected_nontcga = n_det_nt,
    n_tested_nontcga = n_test_nt,
    n_up_nontcga = n_up_nt,
    n_down_nontcga = n_down_nt,
    n_up_nontcga_fc05 = n_up_nt_fc05,
    n_down_nontcga_fc05 = n_down_nt_fc05,
    median_delta_nontcga = med_delta_nt,
    tcga_delta_mean = tcga_d,
    tcga_padj = tcga_q,
    stringsAsFactors = FALSE
  )
})

gene_recurrence_df <- bind_rows(summary_rows)
write.csv(gene_recurrence_df, file.path(res_dir, "gene_recurrence_summary.csv"), row.names = FALSE)

# Save sessionInfo.txt
writeLines(capture.output(sessionInfo()), con = file.path(res_dir, "sessionInfo.txt"))
cat("  Saved: results/cohort_check.csv\n")
cat("  Saved: results/per_cohort_stats.csv\n")
cat("  Saved: results/gene_recurrence_summary.csv\n")
cat("  Saved: results/sessionInfo.txt\n\n")

# Print verification summary for the 7 alias-rescued target genes
cat("[Step 3b] Alias-fallback rescue summary for the 7 specified genes (across 18 included cohorts):\n")
for (ag in names(legacy_alias_map)) {
  sub_g <- per_cohort_stats_df[per_cohort_stats_df$display_label == ag, , drop = FALSE]
  det_g <- sub_g[sub_g$status != "not_detected", , drop = FALSE]
  n_det_all <- nrow(det_g)
  n_det_nt  <- sum(!det_g$is_tcga)
  if (n_det_all > 0L) {
    tb <- table(det_g$matched_symbol)
    hit_str <- paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = ", ")
    alias_only <- det_g[det_g$matched_symbol != ag, , drop = FALSE]
    alias_cohorts <- if (nrow(alias_only) > 0L) {
      paste(sprintf("%s(%s)", alias_only$cohort, alias_only$matched_symbol), collapse = ", ")
    } else {
      "none"
    }
  } else {
    hit_str <- "none"
    alias_cohorts <- "none"
  }
  cat(sprintf("  %-8s | detected=%2d/18 (non-TCGA=%2d/17) | matched=[%s] | alias hits: %s\n",
              ag, n_det_all, n_det_nt, hit_str, alias_cohorts))
}
cat("\n")

# ------------------------------------------------------------------------------
# Step 4. Plot Figure 5 (PDF + 300 dpi PNG)
# ------------------------------------------------------------------------------
cat("[Step 4] Generating Figure 5 (Figure5_HCCDB_heatmap.pdf & 300 dpi .png)...\n")

included_cohorts_df <- cohort_check_df %>%
  filter(included) %>%
  mutate(
    cohort_num = get_cohort_num(cohort),
    is_rnaseq = cohort %in% rnaseq_cohorts
  ) %>%
  arrange(desc(is_rnaseq), cohort_num)

ordered_cohorts <- included_cohorts_df$cohort

# Explicit column labels highlighting HCCDB15 (TCGA-LIHC) and HCCDB25 (Gao 2019)
col_display_labels <- ifelse(
  ordered_cohorts == "HCCDB15", "HCCDB15 (TCGA-LIHC)",
  ifelse(ordered_cohorts == "HCCDB25", "HCCDB25 (Gao 2019)", ordered_cohorts)
)

# Build 49 x N_cohorts matrices for delta_mean and padj
gene_labels <- master_genes$display_label
delta_mat <- matrix(NA_real_, nrow = length(gene_labels), ncol = length(ordered_cohorts),
                    dimnames = list(gene_labels, ordered_cohorts))
padj_mat  <- matrix(NA_real_, nrow = length(gene_labels), ncol = length(ordered_cohorts),
                    dimnames = list(gene_labels, ordered_cohorts))

for (k in seq_len(nrow(per_cohort_stats_df))) {
  r_gene <- per_cohort_stats_df$display_label[k]
  c_coh  <- per_cohort_stats_df$cohort[k]
  if (r_gene %in% gene_labels && c_coh %in% ordered_cohorts) {
    delta_mat[r_gene, c_coh] <- per_cohort_stats_df$delta_mean[k]
    padj_mat[r_gene, c_coh]  <- per_cohort_stats_df$padj[k]
  }
}

# Truncate delta_mean to [-2, 2] for heatmap color mapping while preserving NA
delta_mat_clamped <- pmax(pmin(delta_mat, 2), -2)
delta_mat_clamped[is.na(delta_mat)] <- NA_real_

# Row split by ptm_category (preserving master_gene_list.csv block order; PCDH20 isolated at bottom)
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
row_split_factor <- factor(
  ptm_pretty_map[master_genes$ptm_category],
  levels = unname(ptm_pretty_map[ptm_levels])
)

# Column split: RNA-seq cohorts on the left, Microarray cohorts on the right
n_rnaseq <- sum(included_cohorts_df$is_rnaseq)
n_microarray <- sum(!included_cohorts_df$is_rnaseq)
rnaseq_grp_lbl <- sprintf("RNA-seq Cohorts (n=%d)", n_rnaseq)
microarray_grp_lbl <- sprintf("Microarray Cohorts (n=%d)", n_microarray)

col_split_factor <- factor(
  ifelse(included_cohorts_df$is_rnaseq, rnaseq_grp_lbl, microarray_grp_lbl),
  levels = c(rnaseq_grp_lbl, microarray_grp_lbl)
)

# Derive clean proteomics status directly from master_gene_list.csv$table2_protein_status
raw_prot <- trimws(master_genes$table2_protein_status)
prot_status_clean <- ifelse(
  grepl("^Up_Phospho", raw_prot), "Phospho Up",
  ifelse(
    grepl("^Up_", raw_prot), "Protein Up",
    ifelse(
      grepl("^Down_", raw_prot), "Protein Down",
      "NS / Not detected"
    )
  )
)
prot_levels <- c("Protein Up", "Phospho Up", "Protein Down", "NS / Not detected")
prot_status_factor <- factor(prot_status_clean, levels = prot_levels)

prot_colors <- c(
  "Protein Up"        = "#D6604D",
  "Phospho Up"        = "#FDB863",
  "Protein Down"      = "#4393C3",
  "NS / Not detected" = "#EAEAEA"
)

# Right annotation 1: All 18 cohorts significant Up / Down counts (padj < 0.05)
bar_mat <- cbind(
  Up   = gene_recurrence_df$n_up_all,
  Down = gene_recurrence_df$n_down_all
)
rownames(bar_mat) <- gene_labels
bar_colors <- c("Up" = "#B2182B", "Down" = "#2166AC")

bar_track_title <- sprintf("Sig. Cohorts\n(n=%d)", n_included_all)
right_ha <- rowAnnotation(
  sig_cohorts = anno_barplot(
    bar_mat,
    beside = FALSE,
    gp = gpar(fill = bar_colors, col = NA),
    border = TRUE,
    bar_width = 0.72,
    ylim = c(0, n_included_all),
    axis_param = list(
      side = "bottom",
      at = c(0, 6, 12, 18),
      gp = gpar(fontsize = 8)
    ),
    width = unit(28, "mm")
  ),
  published_proteomics = anno_simple(
    as.character(prot_status_factor),
    col = prot_colors,
    gp = gpar(col = "white", lwd = 0.6),
    width = unit(6.5, "mm")
  ),
  annotation_label = c(bar_track_title, "Published\nProteomics"),
  gap = unit(2.5, "mm"),
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 8.5, fontface = "bold"),
  show_legend = FALSE
)

# Heatmap color function: smooth Blue - White - Red gradient clamped at [-2, 2], NA in light grey
na_fill_col <- "#D9D9D9"
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#FFFFFF", "#B2182B"))

ht <- Heatmap(
  delta_mat_clamped,
  name = "delta_mean",
  col = col_fun,
  na_col = na_fill_col,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  row_split = row_split_factor,
  column_split = col_split_factor,
  row_gap = unit(1.6, "mm"),
  column_gap = unit(3.0, "mm"),
  border = TRUE,
  rect_gp = gpar(col = "#EBEBEB", lwd = 0.5),

  # Row labels & slice titles
  row_labels = gene_labels,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8.5, fontface = "italic"),
  row_title_side = "left",
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 8.5, fontface = "bold", col = "#222222"),

  # Column labels & slice titles
  column_labels = col_display_labels,
  column_names_side = "bottom",
  column_names_rot = 45,
  column_names_gp = gpar(
    fontsize = 8.5,
    fontface = ifelse(ordered_cohorts %in% c("HCCDB15", "HCCDB25"), "bold", "plain")
  ),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),

  right_annotation = right_ha,

  # Overlay asterisk for BH-adjusted P < 0.05 within cohort
  cell_fun = function(j, i, x, y, width, height, fill) {
    qval <- padj_mat[i, j]
    dval <- delta_mat_clamped[i, j]
    if (!is.na(qval) && qval < 0.05) {
      txt_col <- ifelse(!is.na(dval) && abs(dval) > 1.15, "white", "#111111")
      grid.text("*", x = x, y = y - unit(0.45, "mm"),
                gp = gpar(fontsize = 10, fontface = "bold", col = txt_col))
    }
  },

  width = ncol(delta_mat_clamped) * unit(6.8, "mm"),
  height = nrow(delta_mat_clamped) * unit(4.4, "mm"),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "mRNA delta_mean\n(HCC - Adjacent)",
    at = c(-2, -1, 0, 1, 2),
    labels = c("<= -2.0", "-1.0", "0", "+1.0", ">= +2.0"),
    title_gp = gpar(fontsize = 8.5, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(34, "mm"),
    grid_width = unit(4.2, "mm"),
    border = "grey45"
  )
)

# Compact, crisp auxiliary legends (passed flat to draw() to avoid viewport stretching)
lgd_cell_key <- Legend(
  labels = c("padj < 0.05 (within-cohort BH)", "Not detected (NA)"),
  title = "Significance & Detection",
  title_gp = gpar(fontsize = 8.5, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  grid_height = unit(4.0, "mm"),
  grid_width = unit(4.0, "mm"),
  row_gap = unit(1.0, "mm"),
  graphics = list(
    function(x, y, w, h) {
      grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey55"))
      grid.text("*", x, y - unit(0.4, "mm"), gp = gpar(fontsize = 11, fontface = "bold", col = "#111111"))
    },
    function(x, y, w, h) {
      grid.rect(x, y, w, h, gp = gpar(fill = na_fill_col, col = "grey55"))
    }
  )
)

lgd_bar <- Legend(
  labels = c("Sig. Up (padj < 0.05, delta > 0)", "Sig. Down (padj < 0.05, delta < 0)"),
  title = sprintf("Sig. Cohorts (n=%d)", n_included_all),
  title_gp = gpar(fontsize = 8.5, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  grid_height = unit(4.0, "mm"),
  grid_width = unit(4.0, "mm"),
  row_gap = unit(1.0, "mm"),
  legend_gp = gpar(fill = unname(bar_colors)),
  border = "grey55"
)

lgd_prot <- Legend(
  labels = prot_levels,
  title = "Published Proteomics\n/ Phosphoproteomics",
  title_gp = gpar(fontsize = 8.5, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  grid_height = unit(4.0, "mm"),
  grid_width = unit(4.0, "mm"),
  row_gap = unit(1.0, "mm"),
  legend_gp = gpar(fill = unname(prot_colors[prot_levels])),
  border = "grey55"
)

render_figure5 <- function() {
  ht_opt(
    HEATMAP_LEGEND_PADDING = unit(11, "mm"),
    ANNOTATION_LEGEND_PADDING = unit(11, "mm")
  )
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_list = list(lgd_cell_key, lgd_bar, lgd_prot),
    annotation_legend_side = "right",
    merge_legend = TRUE,
    legend_gap = unit(5.5, "mm"),
    padding = unit(c(8, 8, 8, 8), "mm"),
    column_title = sprintf(
      "Transcriptomic Dysregulation of 49 NRF2 PTM Regulators Across %d HCCDB Cohorts (HCC vs. Adjacent)",
      n_included_all
    ),
    column_title_gp = gpar(fontsize = 11.5, fontface = "bold")
  )
}

fig_width_in  <- 13.2
fig_height_in <- 12.6

pdf_path <- file.path(res_dir, "Figure5_HCCDB_heatmap.pdf")
png_path <- file.path(res_dir, "Figure5_HCCDB_heatmap.png")

grDevices::pdf(pdf_path, width = fig_width_in, height = fig_height_in, useDingbats = FALSE)
render_figure5()
grDevices::dev.off()

grDevices::png(png_path, width = fig_width_in, height = fig_height_in, units = "in", res = 300)
render_figure5()
grDevices::dev.off()

cat("  Saved:", pdf_path, "\n")
cat("  Saved:", png_path, "(300 dpi)\n")
cat("\nAll tasks completed successfully!\n")
