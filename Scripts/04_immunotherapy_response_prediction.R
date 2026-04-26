#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE, scipen = 999)
})

# ============================================================
# Rana Salihoglu
#
#   1) Robust CLI validation and explicit reproducibility metadata
#   2) Output root now follows the documented project-root convention
#   3) More conservative phenotype parsing with negative-label priority
#   4) Cohort-wise standardization to limit cross-cohort scale artifacts
#   5) Explicit pathway catalog construction with overlap evidence
#   6) Per-cohort effect estimation + meta-level consistency statistics
#   7) Weighted Stouffer signed meta-analysis for direction-aware synthesis
#   8) Leave-one-cohort-out robustness summaries
#   9) Integrated response-propensity score built from cohort-wise z scores
#  10) Stronger QC, logging, failure messages, and manuscript-facing outputs
#  11) Observed-response benchmarking where labels exist
#  12) Random-effects heterogeneity summaries
#  13) Optional covariate-adjusted association analyses using TME surrogates

# ============================================================

# -----------------------------
# CLI parsing
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  idx <- hit[length(hit)]
  if (idx >= length(args)) return(default)
  args[idx + 1]
}

script_flag <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- tryCatch({
  if (length(script_flag) > 0) normalizePath(sub("^--file=", "", script_flag[1]), winslash = "/", mustWork = FALSE) else NA_character_
}, error = function(e) NA_character_)
script_basename <- if (!is.na(script_path) && nzchar(script_path)) {
  tools::file_path_sans_ext(basename(script_path))
} else {
  "04_immunotherapy_response_prediction_q1_from_uploaded_genes_v6_q1_final"
}

project_root <- path.expand(arg_value("--project_root", "~/Desktop/PROJECT_16/test"))
gene_file <- path.expand(arg_value("--gene_file", "possible_prognostic_genes_FULLTRAIN.csv"))
outdir_name <- arg_value("--outdir_name", script_basename)
outdir_name <- gsub("[^A-Za-z0-9._-]+", "_", outdir_name)
if (!nzchar(outdir_name)) outdir_name <- script_basename
species_name <- arg_value("--species", "Homo sapiens")
path_manifest <- path.expand(arg_value("--path_manifest", "full_paths(2).txt"))
install_missing <- tolower(arg_value("--install_missing", "false")) %in% c("true", "1", "yes", "y")
min_gene_set <- suppressWarnings(as.integer(arg_value("--min_gene_set", "10")))
max_pathways <- suppressWarnings(as.integer(arg_value("--max_pathways", "36")))
min_samples_per_group <- suppressWarnings(as.integer(arg_value("--min_samples_per_group", "4")))
min_labeled_samples <- suppressWarnings(as.integer(arg_value("--min_labeled_samples", "12")))
top_gene_n <- suppressWarnings(as.integer(arg_value("--top_gene_n", "50")))
min_present_genes_per_score <- suppressWarnings(as.integer(arg_value("--min_present_genes_per_score", "2")))
seed <- suppressWarnings(as.integer(arg_value("--seed", "1234")))
set.seed(seed)

if (!is.finite(min_gene_set) || min_gene_set < 5) min_gene_set <- 10
if (!is.finite(max_pathways) || max_pathways < 5) max_pathways <- 36
if (!is.finite(min_samples_per_group) || min_samples_per_group < 2) min_samples_per_group <- 4
if (!is.finite(min_labeled_samples) || min_labeled_samples < 6) min_labeled_samples <- 12
if (!is.finite(top_gene_n) || top_gene_n < 5) top_gene_n <- 50
if (!is.finite(min_present_genes_per_score) || min_present_genes_per_score < 1) min_present_genes_per_score <- 2
if (!is.finite(seed)) seed <- 1234

normalize_maybe <- function(x) {
  tryCatch(normalizePath(path.expand(x), winslash = "/", mustWork = FALSE), error = function(e) path.expand(x))
}

resolve_existing_file <- function(target, candidates = character()) {
  cand <- unique(c(target, candidates))
  cand <- cand[!is.na(cand) & nzchar(cand)]
  for (x in cand) {
    if (file.exists(x) && !dir.exists(x)) return(normalize_maybe(x))
  }
  base <- basename(target)
  if (nzchar(base) && base != target) {
    for (x in cand) {
      parent <- dirname(x)
      alt <- file.path(parent, base)
      if (file.exists(alt) && !dir.exists(alt)) return(normalize_maybe(alt))
    }
  }
  NA_character_
}

read_path_manifest <- function(path_manifest, gene_file, project_root) {
  manifest_candidates <- unique(c(
    path_manifest,
    file.path(getwd(), basename(path_manifest)),
    file.path("/mnt/data", basename(path_manifest)),
    file.path(dirname(normalize_maybe(gene_file)), basename(path_manifest)),
    file.path(dirname(normalize_maybe(project_root)), basename(path_manifest))
  ))
  mf <- resolve_existing_file(path_manifest, manifest_candidates)
  if (is.na(mf)) return(character())
  lines <- tryCatch(readLines(mf, warn = FALSE), error = function(e) character())
  lines <- trimws(lines)
  unique(lines[nzchar(lines)])
}

manifest_paths <- read_path_manifest(path_manifest, gene_file, project_root)

infer_project_root <- function(project_root, manifest_paths, gene_file) {
  direct_candidates <- unique(c(
    project_root,
    dirname(gene_file),
    dirname(normalize_maybe(gene_file)),
    dirname(file.path(getwd(), gene_file)),
    dirname(file.path("/mnt/data", basename(gene_file)))
  ))
  for (cand in direct_candidates) {
    if (is.na(cand) || !nzchar(cand)) next
    cand_n <- normalize_maybe(cand)
    looks_like_root <- dir.exists(file.path(cand_n, "OUT_v44")) ||
      dir.exists(file.path(cand_n, "GEO_PREP")) ||
      file.exists(file.path(cand_n, "ferrdb_driver.csv"))
    if (looks_like_root) return(cand_n)
  }
  if (length(manifest_paths) > 0) {
    root_hits <- unique(dirname(dirname(manifest_paths[grepl("/(OUT_v44|GEO_PREP|bioinf_scripts|scripts)/", manifest_paths)])))
    root_hits <- root_hits[nzchar(root_hits)]
    if (length(root_hits) > 0) return(normalize_maybe(root_hits[[1]]))
  }
  normalize_maybe(project_root)
}

project_root <- infer_project_root(project_root, manifest_paths, gene_file)

gene_candidates <- unique(c(
  gene_file,
  file.path(getwd(), gene_file),
  file.path(getwd(), basename(gene_file)),
  file.path("/mnt/data", gene_file),
  file.path("/mnt/data", basename(gene_file)),
  file.path(project_root, basename(gene_file)),
  file.path(project_root, "OUT_v44", "genes", basename(gene_file)),
  file.path(project_root, "GEO_PREP", "GSE81089", "WGCNA_GSE81089_out_v12_signature_fix", "tables", basename(gene_file)),
  manifest_paths[basename(manifest_paths) == basename(gene_file)],
  manifest_paths[grepl(paste0("/", basename(gene_file), "$"), manifest_paths)]
))
resolved_gene_file <- resolve_existing_file(gene_file, gene_candidates)
if (!is.na(resolved_gene_file)) gene_file <- resolved_gene_file

cat("[INFO] project_root :", project_root, "\n")
cat("[INFO] gene_file    :", gene_file, "\n")
cat("[INFO] path_manifest:", ifelse(length(manifest_paths) > 0, "loaded", "not_found_or_empty"), "\n")
cat("[INFO] outdir_name  :", outdir_name, "\n")
cat("[INFO] script_name  :", script_basename, "\n")
cat("[INFO] seed         :", seed, "\n")

# -----------------------------
# Package management
# -----------------------------
cran_pkgs <- c(
  "data.table", "dplyr", "tibble", "tidyr", "stringr", "purrr", "readr",
  "ggplot2", "forcats", "scales", "pheatmap", "jsonlite", "glue"
)
bioc_pkgs <- c("GSVA", "msigdbr")
optional_pkgs <- c("openxlsx", "ggrepel", "pROC", "IOBR")

ensure_pkg <- function(pkg, bioc = FALSE, install_missing = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  if (!install_missing) return(FALSE)
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (bioc || pkg %in% c("GSVA", "msigdbr", "IOBR")) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  requireNamespace(pkg, quietly = TRUE)
}

missing_core <- c(
  cran_pkgs[!vapply(cran_pkgs, ensure_pkg, logical(1), bioc = FALSE, install_missing = install_missing)],
  bioc_pkgs[!vapply(bioc_pkgs, ensure_pkg, logical(1), bioc = TRUE, install_missing = install_missing)]
)
if (length(missing_core) > 0) {
  stop("Required packages missing: ", paste(missing_core, collapse = ", "),
       ". Re-run with --install_missing true or install them manually.")
}
invisible(vapply(optional_pkgs, ensure_pkg, logical(1), bioc = FALSE, install_missing = install_missing))

library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(purrr)
library(readr)
library(ggplot2)
library(forcats)
library(scales)
library(pheatmap)
library(jsonlite)
library(glue)
library(GSVA)
library(msigdbr)

has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_pROC <- requireNamespace("pROC", quietly = TRUE)
has_iobr <- requireNamespace("IOBR", quietly = TRUE)

# -----------------------------
# Paths
# -----------------------------
root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
out_dir <- file.path(root, outdir_name)
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(fig_dir) && !dir.exists(fig_dir)) stop("fig_dir exists but is not a directory: ", fig_dir)
if (file.exists(tab_dir) && !dir.exists(tab_dir)) stop("tab_dir exists but is not a directory: ", tab_dir)
if (file.exists(log_dir) && !dir.exists(log_dir)) stop("log_dir exists but is not a directory: ", log_dir)

paths <- list(
  geo_root = file.path(root, "GEO_PREP"),
  out_v44 = file.path(root, "OUT_v44"),
  ferr_driver = file.path(root, "ferrdb_driver.csv"),
  ferr_suppressor = file.path(root, "ferrdb_suppressor.csv"),
  ferr_marker = file.path(root, "ferrdb_marker.csv"),
  wgcna_tables = file.path(root, "GEO_PREP", "GSE81089", "WGCNA_GSE81089_out_v12_signature_fix", "tables"),
  enrich_hallmark = file.path(root, "pathway_enrichment_Q1_v2", "tables", "table_gsea_meta_hallmark.tsv"),
  enrich_integrated = file.path(root,  "pathway_enrichment_Q1_v2", "tables", "table_integrated_pathway_summary.tsv")
)

# -----------------------------
# Helpers
# -----------------------------
msg <- function(...) cat(glue(..., .envir = parent.frame()), "\n")

safe_dir_create <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) stop("Invalid directory path")
  if (file.exists(path) && !dir.exists(path)) stop("Path exists but is not a directory: ", path)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Could not create directory: ", path)
  invisible(path)
}

read_flex <- function(path, ...) {
  if (is.null(path) || is.na(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    fread(path, sep = "\t", data.table = FALSE, ...)
  } else {
    fread(path, data.table = FALSE, ...)
  }
}

safe_write <- function(df, path) {
  if (is.null(df)) return(invisible(NULL))
  safe_dir_create(dirname(path))
  readr::write_tsv(as_tibble(df), path)
}

sanitize_filename <- function(x, max_n = 180) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_+", "", x)
  x <- sub("_+$", "", x)
  if (!nzchar(x)) x <- "figure"
  if (nchar(x) > max_n) x <- substr(x, 1, max_n)
  x
}

plot_save <- function(p, file_stub, width = 10, height = 7) {
  safe_dir_create(dirname(file_stub))
  file_stub <- file.path(dirname(file_stub), sanitize_filename(basename(file_stub)))
  ggplot2::ggsave(paste0(file_stub, ".png"), plot = p, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(paste0(file_stub, ".pdf"), plot = p, width = width, height = height, bg = "white")
}

safe_pheatmap_save <- function(mat, file_stub, width = 2200, height = 1600, res = 220, ...) {
  safe_dir_create(dirname(file_stub))
  file_stub <- file.path(dirname(file_stub), sanitize_filename(basename(file_stub)))
  out_png <- paste0(file_stub, ".png")
  out_pdf <- paste0(file_stub, ".pdf")
  tryCatch({
    pheatmap::pheatmap(mat, filename = out_png, width = width / res, height = height / res, ...)
  }, error = function(e) {
    grDevices::png(out_png, width = width, height = height, res = res)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    pheatmap::pheatmap(mat, ...)
  })
  tryCatch({
    grDevices::pdf(out_pdf, width = width / res, height = height / res, onefile = TRUE)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    pheatmap::pheatmap(mat, ...)
  }, error = function(e) {
    message("[WARN] PDF heatmap write failed: ", out_pdf, " :: ", conditionMessage(e))
  })
  invisible(c(out_png, out_pdf))
}

first_existing <- function(vec) {
  hit <- vec[file.exists(vec)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

clean_names <- function(x) make.names(x, unique = TRUE)

cohort_z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (sum(is.finite(x)) < 3) return(rep(NA_real_, length(x)))
  z <- as.numeric(scale(x))
  z[!is.finite(z)] <- NA_real_
  z
}

cohort_rank01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) < 2) return(out)
  r <- rank(x[ok], ties.method = "average")
  out[ok] <- (r - 1) / max(1, sum(ok) - 1)
  out
}

safe_min_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_max_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_neglog10 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  out[ok] <- -log10(pmax(x[ok], 1e-300))
  out
}

safe_rescale01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(0, length(x))
  ok <- is.finite(x)
  if (sum(ok) < 2) return(out)
  rx <- range(x[ok], na.rm = TRUE)
  if (!all(is.finite(rx)) || diff(rx) == 0) {
    out[ok] <- 0.5
    return(out)
  }
  out[ok] <- (x[ok] - rx[1]) / diff(rx)
  out
}

top_n_rows <- function(df, n = 20, order_col = NULL, decreasing = TRUE) {
  if (is.null(df) || nrow(df) == 0) return(df)
  n_take <- min(n, nrow(df))
  if (is.null(order_col) || !order_col %in% names(df)) return(df[seq_len(n_take), , drop = FALSE])
  ord <- order(df[[order_col]], decreasing = decreasing, na.last = TRUE)
  df[ord[seq_len(n_take)], , drop = FALSE]
}

safe_factor_reorder <- function(labels, values, desc = FALSE) {
  x <- tibble(label = as.character(labels), value = suppressWarnings(as.numeric(values))) %>%
    filter(!is.na(label), label != "") %>%
    group_by(label) %>%
    summarise(order_value = mean(value[is.finite(value)], na.rm = TRUE), .groups = "drop")
  if (nrow(x) == 0) return(factor(labels))
  x$order_value[!is.finite(x$order_value)] <- 0
  ord <- x$label[order(x$order_value, decreasing = desc)]
  factor(as.character(labels), levels = ord)
}

clean_feature_label <- function(x, max_len = 72) {
  x <- gsub("_", " ", as.character(x), fixed = TRUE)
  x <- gsub("\\s+", " ", x)
  stringr::str_trim(stringr::str_trunc(x, width = max_len))
}

compact_pathway_label <- function(gs_name, gs_description = NA_character_, collection = NA_character_) {
  raw <- paste(gs_name, gs_description)
  raw_low <- tolower(raw)
  out <- dplyr::case_when(
    str_detect(raw_low, "alpha interferon|response to alpha interferon|interferon alpha") ~ "IFN-alpha response",
    str_detect(raw_low, "ifng|gamma interferon|interferon gamma") ~ "IFN-gamma response",
    str_detect(raw_low, "tgfb|transforming growth factor beta") ~ "TGF-beta response",
    str_detect(raw_low, "angiogenesis|blood vessels") ~ "Angiogenesis",
    str_detect(raw_low, "complement") ~ "Complement",
    str_detect(raw_low, "inflammatory response|inflam") ~ "Inflammatory response",
    str_detect(raw_low, "g2/m checkpoint|g2m") ~ "G2/M checkpoint",
    str_detect(raw_low, "plasmacytoid dendritic") & str_detect(raw_low, "il3") ~ "pDC IL3 program",
    str_detect(raw_low, "plasmacytoid dendritic") ~ "Plasmacytoid dendritic cell program",
    str_detect(raw_low, "memory cd8") & str_detect(raw_low, "effector") ~ "Memory vs effector CD8 program",
    str_detect(raw_low, "effector cd8") & str_detect(raw_low, "memory") ~ "Effector vs memory CD8 program",
    str_detect(raw_low, "naive cd8") & str_detect(raw_low, "memory") ~ "Naive vs memory CD8 program",
    str_detect(raw_low, "naive cd8") & str_detect(raw_low, "effector") ~ "Naive vs effector CD8 program",
    str_detect(raw_low, "krlg1") & str_detect(raw_low, "effector cd8") ~ "KLRG1-int effector CD8 program",
    str_detect(raw_low, "marginal zone b cells") ~ "Marginal zone B-cell program",
    str_detect(raw_low, "macrophages") & str_detect(raw_low, "m-csf") ~ "Macrophage M-CSF program",
    str_detect(raw_low, "pbmc|peripheral blood mononuclear") & str_detect(raw_low, "14d|14 d|day 14") ~ "PBMC day14 response",
    str_detect(raw_low, "pbmc|peripheral blood mononuclear") & str_detect(raw_low, "10d|10 d|day 10") ~ "PBMC day10 response",
    str_detect(raw_low, "pbmc|peripheral blood mononuclear") & str_detect(raw_low, "7d|7 d|day 7") ~ "PBMC day7 response",
    str_detect(raw_low, "pbmc|peripheral blood mononuclear") & str_detect(raw_low, "3d|3 d|day 3") ~ "PBMC day3 response",
    str_detect(raw_low, "healthy cd4") ~ "Healthy CD4 T-cell program",
    str_detect(raw_low, "t cell") ~ clean_feature_label(gs_name, max_len = 48),
    TRUE ~ clean_feature_label(ifelse(is.na(gs_description) | gs_description == "", gs_name, gs_description), max_len = 48)
  )
  out <- gsub("\\s+", " ", out)
  out <- stringr::str_trim(out)
  if (!nzchar(out)) out <- clean_feature_label(gs_name, max_len = 48)
  out
}

pathway_family_label <- function(gs_name, gs_description = NA_character_) {
  raw_low <- tolower(paste(gs_name, gs_description))
  dplyr::case_when(
    str_detect(raw_low, "interferon|ifng") ~ "Interferon",
    str_detect(raw_low, "tgfb") ~ "TGF-beta",
    str_detect(raw_low, "angiogenesis|blood vessels|strom|fibroblast") ~ "Stromal/angiogenic",
    str_detect(raw_low, "complement|inflam|chemokine") ~ "Inflammation/complement",
    str_detect(raw_low, "dendritic|macrophage|myeloid|nk|b cells") ~ "Innate/myeloid/B-cell",
    str_detect(raw_low, "cd8|t cell|antigen|mhc|hla|cytotoxic|checkpoint|exhaust") ~ "T-cell adaptive",
    TRUE ~ "Other immune"
  )
}

jaccard_index <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  inter <- length(intersect(a, b))
  union <- length(unique(c(a, b)))
  if (union == 0) return(0)
  inter / union
}

prune_redundant_gene_sets <- function(tbl, max_n = 30, max_jaccard = 0.75) {
  if (is.null(tbl) || nrow(tbl) == 0) return(tbl)
  keep <- logical(nrow(tbl))
  selected <- integer(0)
  for (i in seq_len(nrow(tbl))) {
    if (length(selected) == 0) {
      keep[i] <- TRUE
      selected <- c(selected, i)
      next
    }
    this_genes <- tbl$genes[[i]]
    max_ol <- max(vapply(selected, function(j) jaccard_index(this_genes, tbl$genes[[j]]), numeric(1)), na.rm = TRUE)
    if (!is.finite(max_ol) || max_ol < max_jaccard) {
      keep[i] <- TRUE
      selected <- c(selected, i)
    }
    if (sum(keep) >= max_n) break
  }
  tbl[keep, , drop = FALSE]
}

find_gene_col <- function(df) {
  cand <- names(df)[tolower(names(df)) %in% c("gene", "symbol", "gene_symbol", "genesymbol", "official_symbol", "hub_gene", "x")]
  if (length(cand) > 0) cand[[1]] else names(df)[[1]]
}

looks_like_gene_symbol <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(FALSE)
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  if (length(x) == 0) return(FALSE)
  mean(str_detect(head(x, 200), "^[A-Za-z0-9._-]{2,25}$")) >= 0.7
}

looks_like_sample_id <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(FALSE)
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  if (length(x) == 0) return(FALSE)
  mean(str_detect(head(x, 200), "^(GSM|TCGA|SRR|ERR|SAM|[A-Za-z]{1,8}[-_][A-Za-z0-9._-]+|[A-Za-z0-9._-]{4,})$")) >= 0.5
}

detect_sample_col <- function(df) {
  if (is.null(df) || ncol(df) == 0) return(names(df)[1])
  nm <- names(df)
  low <- tolower(nm)
  idx <- which(str_detect(low, "sample|sample_id|sampleid|patient|patient_id|patientid|specimen|barcode|submitter|gsm|id$"))
  if (length(idx) > 0) return(nm[idx[1]])
  nm[1]
}

collapse_duplicate_genes <- function(mat, gene_names) {
  gene_names <- clean_names(gene_names)
  ok <- !is.na(gene_names) & gene_names != ""
  mat <- mat[ok, , drop = FALSE]
  gene_names <- gene_names[ok]
  if (nrow(mat) == 0) return(mat)
  split_idx <- split(seq_along(gene_names), gene_names)
  out <- vapply(split_idx, function(idx) {
    if (length(idx) == 1) mat[idx, ] else colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
  }, FUN.VALUE = numeric(ncol(mat)))
  out <- t(out)
  colnames(out) <- colnames(mat)
  rownames(out) <- names(split_idx)
  mode(out) <- "numeric"
  out
}

parse_expression_table <- function(expr_df) {
  if (is.null(expr_df) || nrow(expr_df) < 5 || ncol(expr_df) < 5) return(NULL)
  expr_df <- as.data.frame(expr_df, stringsAsFactors = FALSE)
  nm <- names(expr_df)
  first_col <- nm[1]
  first_vals <- expr_df[[first_col]]

  if (looks_like_sample_id(first_vals)) {
    sample_col <- detect_sample_col(expr_df)
    gene_cols <- setdiff(names(expr_df), sample_col)
    dat <- expr_df[, gene_cols, drop = FALSE]
    for (k in names(dat)) dat[[k]] <- suppressWarnings(as.numeric(dat[[k]]))
    keep_cols <- names(dat)[vapply(dat, function(x) any(is.finite(x)), logical(1))]
    if (length(keep_cols) >= 50) {
      mat <- as.matrix(dat[, keep_cols, drop = FALSE])
      rownames(mat) <- as.character(expr_df[[sample_col]])
      mode(mat) <- "numeric"
      colnames(mat) <- clean_names(colnames(mat))
      mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
      return(mat)
    }
  }

  if (looks_like_gene_symbol(first_vals)) {
    gene_col <- first_col
    sample_cols <- setdiff(names(expr_df), gene_col)
    dat <- expr_df[, sample_cols, drop = FALSE]
    for (k in names(dat)) dat[[k]] <- suppressWarnings(as.numeric(dat[[k]]))
    keep_cols <- names(dat)[vapply(dat, function(x) any(is.finite(x)), logical(1))]
    if (length(keep_cols) >= 5) {
      mat <- as.matrix(dat[, keep_cols, drop = FALSE])
      mode(mat) <- "numeric"
      rownames(mat) <- as.character(expr_df[[gene_col]])
      mat <- collapse_duplicate_genes(mat, rownames(mat))
      out <- t(mat)
      out <- out[, colSums(is.finite(out)) > 0, drop = FALSE]
      colnames(out) <- clean_names(colnames(out))
      rownames(out) <- keep_cols
      out <- out[!duplicated(rownames(out)), , drop = FALSE]
      return(out)
    }
  }

  sample_col <- detect_sample_col(expr_df)
  gene_cols <- setdiff(names(expr_df), sample_col)
  dat <- expr_df[, gene_cols, drop = FALSE]
  for (k in names(dat)) dat[[k]] <- suppressWarnings(as.numeric(dat[[k]]))
  keep_cols <- names(dat)[vapply(dat, function(x) any(is.finite(x)), logical(1))]
  if (length(keep_cols) >= 50) {
    mat <- as.matrix(dat[, keep_cols, drop = FALSE])
    rownames(mat) <- as.character(expr_df[[sample_col]])
    mode(mat) <- "numeric"
    colnames(mat) <- clean_names(colnames(mat))
    mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
    return(mat)
  }
  NULL
}

cohort_row_z <- function(mat) {
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(mat)
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

simple_ssgsea <- function(expr_mat, gene_sets, min_gene_set = 10, alpha = 0.25) {
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= min_gene_set]
  if (length(gene_sets) == 0) return(matrix(numeric(), nrow = 0, ncol = ncol(expr_mat), dimnames = list(character(), colnames(expr_mat))))
  ranks <- apply(expr_mat, 2, rank, ties.method = "average")
  rownames(ranks) <- rownames(expr_mat)
  out <- matrix(NA_real_, nrow = length(gene_sets), ncol = ncol(expr_mat), dimnames = list(names(gene_sets), colnames(expr_mat)))
  all_genes <- rownames(expr_mat)
  for (i in seq_along(gene_sets)) {
    gs <- intersect(unique(gene_sets[[i]]), all_genes)
    if (length(gs) < min_gene_set) next
    idx <- all_genes %in% gs
    Nh <- sum(idx)
    Nm <- length(all_genes) - Nh
    if (Nh < 1 || Nm < 1) next
    for (j in seq_len(ncol(expr_mat))) {
      r <- ranks[, j]
      ord <- order(r, decreasing = TRUE)
      hits <- idx[ord]
      ranked_r <- abs(r[ord]) ^ alpha
      denom <- sum(ranked_r[hits])
      if (!is.finite(denom) || denom <= 0) next
      Phit <- cumsum(ifelse(hits, ranked_r / denom, 0))
      Pmiss <- cumsum(ifelse(!hits, 1 / Nm, 0))
      out[i, j] <- sum(Phit - Pmiss)
    }
  }
  out
}

safe_gsva_ssgsea <- function(expr_mat, gene_sets, min_gene_set = 10) {
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= min_gene_set]
  if (length(gene_sets) == 0) {
    return(matrix(numeric(), nrow = 0, ncol = ncol(expr_mat), dimnames = list(character(), colnames(expr_mat))))
  }
  expr_mat <- as.matrix(expr_mat)
  mode(expr_mat) <- "numeric"
  expr_mat[!is.finite(expr_mat)] <- NA_real_

  res <- tryCatch({
    if (exists("ssgseaParam", where = asNamespace("GSVA"), mode = "function")) {
      par_obj <- GSVA::ssgseaParam(exprData = expr_mat, geneSets = gene_sets, minSize = min_gene_set, alpha = 0.25, normalize = TRUE)
      suppressWarnings(GSVA::gsva(par_obj))
    } else {
      suppressWarnings(GSVA::gsva(expr_mat, gene_sets, method = "ssgsea", ssgsea.norm = TRUE, kcdf = "Gaussian", abs.ranking = TRUE))
    }
  }, error = function(e) NULL)

  if (is.null(res)) {
    msg("[WARN] GSVA failed; falling back to internal ranking scorer")
    return(simple_ssgsea(expr_mat, gene_sets, min_gene_set = min_gene_set))
  }
  res
}

compute_auc_safe <- function(labels01, scores) {
  if (!has_pROC) return(NA_real_)
  ok <- is.finite(scores) & !is.na(labels01)
  labels01 <- labels01[ok]
  scores <- scores[ok]
  if (length(unique(labels01)) < 2) return(NA_real_)
  tryCatch(as.numeric(pROC::auc(labels01, scores, quiet = TRUE)), error = function(e) NA_real_)
}

clean_group_values <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- str_replace_all(x, "[_./-]", " ")
  x <- str_squish(x)
  x
}

response_patterns <- list(
  pos = c(
    "complete response", "partial response", "durable clinical benefit", "clinical benefit",
    "objective response", "response", "responder", "benefit", "sensitive", "dcb",
    "^cr$", "^pr$", "^yes$", "^1$"
  ),
  neg = c(
    "progressive disease", "stable disease", "no response", "non response", "nonresponder",
    "non responder", "no benefit", "resistant", "ncb", "^sd$", "^pd$", "^nr$", "^no$", "^0$"
  )
)

infer_response_from_pheno <- function(pheno_df) {
  if (is.null(pheno_df) || nrow(pheno_df) == 0) return(NULL)
  nms <- names(pheno_df)
  nms_low <- tolower(nms)

  response_col_idx <- unique(c(
    grep("response|responder|recist|benefit|clinical.*benefit|bor|best.*overall.*response|best.response|therapy.*response|ici.*response|immun.*response|anti.pd|anti pd|anti.pdl1|anti pdl1|anti.ctla|anti ctla|outcome|efficacy|dcb|ncb|status|label", nms_low)
  ))
  if (length(response_col_idx) == 0) response_col_idx <- seq_along(nms)

  sample_col_idx <- unique(c(
    grep("sample|sample_id|gsm|geo_accession|geo accession|submitter|patient|patient_id|case|id$", nms_low)
  ))
  sample_col <- if (length(sample_col_idx) > 0) nms[[sample_col_idx[1]]] else nms[[1]]

  best <- NULL
  best_hits <- -1
  for (col in nms[response_col_idx]) {
    vals_raw <- as.character(pheno_df[[col]])
    vals <- clean_group_values(vals_raw)
    response01 <- rep(NA_integer_, length(vals))

    neg_hit <- rep(FALSE, length(vals))
    pos_hit <- rep(FALSE, length(vals))
    for (pat in response_patterns$neg) neg_hit <- neg_hit | str_detect(vals, regex(pat))
    for (pat in response_patterns$pos) pos_hit <- pos_hit | str_detect(vals, regex(pat))

    response01[neg_hit] <- 0L
    response01[pos_hit & !neg_hit] <- 1L

    numv <- suppressWarnings(as.numeric(vals_raw))
    uniq_num <- sort(unique(numv[is.finite(numv)]))
    if (length(uniq_num) == 2 && all(uniq_num %in% c(0, 1))) {
      response01[is.finite(numv)] <- as.integer(numv[is.finite(numv)])
    }

    n_hits <- sum(!is.na(response01))
    if (n_hits > best_hits && length(unique(stats::na.omit(response01))) >= 2) {
      best_hits <- n_hits
      best <- tibble(
        sample_id = as.character(pheno_df[[sample_col]]),
        response01 = response01,
        response_label = ifelse(response01 == 1L, "Responder", ifelse(response01 == 0L, "NonResponder", NA_character_)),
        response_source_col = col,
        response_mode = "observed"
      )
    }
  }
  best
}

compute_binary_effects <- function(df, value_col = "score", response_col = "response01", min_samples_per_group = 4) {
  x <- as_tibble(df)
  vals <- suppressWarnings(as.numeric(x[[value_col]]))
  y <- suppressWarnings(as.integer(x[[response_col]]))
  ok <- is.finite(vals) & !is.na(y)
  vals <- vals[ok]
  y <- y[ok]
  if (length(unique(y)) < 2) {
    return(c(p_value = NA_real_, delta = NA_real_, auc = NA_real_, odds_ratio = NA_real_, beta = NA_real_, mean_resp = NA_real_, mean_nonresp = NA_real_, rank_biserial = NA_real_, n_total = length(vals), n_resp = sum(y == 1), n_nonresp = sum(y == 0)))
  }
  resp <- vals[y == 1]
  nonresp <- vals[y == 0]
  if (length(resp) < min_samples_per_group || length(nonresp) < min_samples_per_group) {
    return(c(p_value = NA_real_, delta = mean(resp, na.rm = TRUE) - mean(nonresp, na.rm = TRUE), auc = NA_real_, odds_ratio = NA_real_, beta = NA_real_, mean_resp = mean(resp, na.rm = TRUE), mean_nonresp = mean(nonresp, na.rm = TRUE), rank_biserial = NA_real_, n_total = length(vals), n_resp = length(resp), n_nonresp = length(nonresp)))
  }
  wt <- tryCatch(suppressWarnings(wilcox.test(resp, nonresp, exact = FALSE)), error = function(e) NULL)
  auc <- compute_auc_safe(y, vals)
  rank_biserial <- if (is.finite(auc)) 2 * auc - 1 else NA_real_
  fit <- tryCatch(glm(y ~ vals, family = binomial()), error = function(e) NULL)
  beta <- if (!is.null(fit) && length(coef(fit)) >= 2) suppressWarnings(as.numeric(coef(fit)[2])) else NA_real_
  or <- if (is.finite(beta)) exp(beta) else NA_real_
  c(
    p_value = if (!is.null(wt)) wt$p.value else NA_real_,
    delta = mean(resp, na.rm = TRUE) - mean(nonresp, na.rm = TRUE),
    auc = auc,
    odds_ratio = or,
    beta = beta,
    mean_resp = mean(resp, na.rm = TRUE),
    mean_nonresp = mean(nonresp, na.rm = TRUE),
    rank_biserial = rank_biserial,
    n_total = length(vals),
    n_resp = length(resp),
    n_nonresp = length(nonresp)
  )
}

compute_two_group_effects <- function(values, group, min_samples_per_group = 4, group1 = NULL, group0 = NULL) {
  vals <- suppressWarnings(as.numeric(values))
  grp <- as.character(group)
  ok <- is.finite(vals) & !is.na(grp) & grp != ""
  vals <- vals[ok]
  grp <- grp[ok]
  if (length(unique(grp)) < 2) {
    return(c(p_value = NA_real_, delta = NA_real_, rank_biserial = NA_real_, n_total = length(vals), n_group1 = NA_real_, n_group0 = NA_real_))
  }
  lv <- sort(unique(grp))
  if (is.null(group0) || !group0 %in% lv) group0 <- lv[1]
  if (is.null(group1) || !group1 %in% lv) group1 <- setdiff(lv, group0)[1]
  if (!is.finite(match(group1, lv)) || !is.finite(match(group0, lv)) || identical(group1, group0)) {
    group0 <- lv[1]
    group1 <- lv[2]
  }
  x0 <- vals[grp == group0]
  x1 <- vals[grp == group1]
  if (length(x0) < min_samples_per_group || length(x1) < min_samples_per_group) {
    return(c(p_value = NA_real_, delta = mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE), rank_biserial = NA_real_, n_total = length(vals), n_group1 = length(x1), n_group0 = length(x0)))
  }
  wt <- tryCatch(suppressWarnings(wilcox.test(x1, x0, exact = FALSE)), error = function(e) NULL)
  U <- if (!is.null(wt) && !is.null(wt$statistic)) as.numeric(wt$statistic) - length(x1) * (length(x1) + 1) / 2 else NA_real_
  rbc <- if (is.finite(U) && length(x1) > 0 && length(x0) > 0) 2 * U / (length(x1) * length(x0)) - 1 else NA_real_
  c(
    p_value = if (!is.null(wt)) wt$p.value else NA_real_,
    delta = mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE),
    rank_biserial = rbc,
    n_total = length(vals),
    n_group1 = length(x1),
    n_group0 = length(x0)
  )
}

meta_signed_stouffer <- function(pvals, effects, weights = NULL) {
  pvals <- suppressWarnings(as.numeric(pvals))
  effects <- suppressWarnings(as.numeric(effects))
  if (is.null(weights)) weights <- rep(1, length(pvals))
  weights <- suppressWarnings(as.numeric(weights))
  ok <- is.finite(pvals) & pvals > 0 & pvals <= 1 & is.finite(effects) & is.finite(weights) & weights > 0
  if (sum(ok) == 0) return(c(meta_p = NA_real_, meta_z = NA_real_, direction = NA_real_))
  z <- qnorm(pmax(1e-300, pvals[ok]) / 2, lower.tail = FALSE)
  s <- sign(effects[ok])
  s[s == 0] <- 1
  wz <- sum(weights[ok] * s * z) / sqrt(sum(weights[ok]^2))
  c(meta_p = 2 * pnorm(-abs(wz)), meta_z = wz, direction = sign(wz))
}

loo_meta_rank <- function(df, feature_col, effect_col, p_col, weight_col) {
  feats <- unique(df[[feature_col]])
  cohorts <- unique(df$cohort)
  if (length(cohorts) < 2) return(NULL)
  out <- list()
  idx <- 1L
  for (coh in cohorts) {
    sub <- df[df$cohort != coh, , drop = FALSE]
    if (nrow(sub) == 0) next
    tmp <- sub %>%
      group_by(.data[[feature_col]]) %>%
      summarise(
        meta_p = meta_signed_stouffer(.data[[p_col]], .data[[effect_col]], .data[[weight_col]])[["meta_p"]],
        mean_effect = mean(.data[[effect_col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(rank_score = safe_neglog10(meta_p) + abs(mean_effect))
    tmp$left_out_cohort <- coh
    names(tmp)[1] <- feature_col
    out[[idx]] <- tmp
    idx <- idx + 1L
  }
  bind_rows(out)
}



approx_meta_heterogeneity <- function(effect, n_case, n_control, fallback_n = NULL) {
  yi <- suppressWarnings(as.numeric(effect))
  n1 <- suppressWarnings(as.numeric(n_case))
  n0 <- suppressWarnings(as.numeric(n_control))
  nt <- suppressWarnings(as.numeric(fallback_n))
  vi <- rep(NA_real_, length(yi))
  ok2 <- is.finite(n1) & n1 > 1 & is.finite(n0) & n0 > 1
  vi[ok2] <- 1 / n1[ok2] + 1 / n0[ok2]
  ok1 <- !ok2 & is.finite(nt) & nt > 2
  vi[ok1] <- 4 / nt[ok1]
  ok <- is.finite(yi) & is.finite(vi) & vi > 0
  if (sum(ok) == 0) {
    return(c(k = 0, fixed_beta = NA_real_, fixed_se = NA_real_, fixed_p = NA_real_, random_beta = NA_real_, random_se = NA_real_, random_p = NA_real_, Q = NA_real_, Q_p = NA_real_, I2 = NA_real_, tau2 = NA_real_))
  }
  yi <- yi[ok]
  vi <- vi[ok]
  wi <- 1 / vi
  fixed_beta <- sum(wi * yi) / sum(wi)
  fixed_se <- sqrt(1 / sum(wi))
  Q <- sum(wi * (yi - fixed_beta)^2)
  dfq <- length(yi) - 1
  C <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (dfq > 0 && is.finite(C) && C > 0) max(0, (Q - dfq) / C) else 0
  wi_re <- 1 / (vi + tau2)
  random_beta <- sum(wi_re * yi) / sum(wi_re)
  random_se <- sqrt(1 / sum(wi_re))
  I2 <- if (dfq > 0 && is.finite(Q) && Q > 0) max(0, (Q - dfq) / Q) else 0
  c(
    k = length(yi),
    fixed_beta = fixed_beta,
    fixed_se = fixed_se,
    fixed_p = 2 * pnorm(-abs(fixed_beta / fixed_se)),
    random_beta = random_beta,
    random_se = random_se,
    random_p = 2 * pnorm(-abs(random_beta / random_se)),
    Q = Q,
    Q_p = if (dfq > 0) pchisq(Q, df = dfq, lower.tail = FALSE) else NA_real_,
    I2 = I2,
    tau2 = tau2
  )
}

select_adjustment_covariates <- function(deconv_long, max_covars = 4) {
  if (is.null(deconv_long) || nrow(deconv_long) == 0) return(character())
  patt <- c('cd8', 'cytotoxic', 't cell', 'immune', 'estimate', 'stromal', 'macro', 'fibro', 'myeloid', 'nk')
  cand <- deconv_long %>%
    mutate(feature_low = tolower(feature)) %>%
    filter(str_detect(feature_low, paste(patt, collapse = '|'))) %>%
    group_by(feature) %>%
    summarise(
      n_cohorts = n_distinct(cohort),
      n_samples = sum(is.finite(score_raw)),
      sd_score = sd(score_raw, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_cohorts), desc(n_samples), desc(sd_score))
  head(cand$feature, max_covars)
}

compute_adjusted_feature_associations <- function(long_df, deconv_long, analysis_mode, feature_col = 'feature', value_col = 'score_z', feature_subset = NULL, max_covars = 4, min_n = 20) {
  if (is.null(long_df) || nrow(long_df) == 0) return(NULL)
  covars <- select_adjustment_covariates(deconv_long, max_covars = max_covars)
  if (length(covars) == 0) return(NULL)
  deconv_cov <- deconv_long %>%
    filter(feature %in% covars) %>%
    select(sample_id, cohort, feature, score_z) %>%
    distinct() %>%
    pivot_wider(names_from = feature, values_from = score_z)
  dat <- long_df
  if (!is.null(feature_subset)) dat <- dat %>% filter(.data[[feature_col]] %in% feature_subset)
  if (nrow(dat) == 0) return(NULL)
  pieces <- split(dat, list(dat$cohort, dat[[feature_col]]), drop = TRUE)
  out <- bind_rows(lapply(pieces, function(df) {
    if (nrow(df) < min_n) return(NULL)
    feat_name <- unique(df[[feature_col]])[1]
    tmp <- df %>%
      select(sample_id, cohort, value = all_of(value_col), response01, integrated_benefit_score) %>%
      left_join(deconv_cov, by = c('sample_id', 'cohort'))
    keep_cov <- covars[covars %in% names(tmp)]
    cc <- complete.cases(tmp[, unique(c('value', keep_cov, if (analysis_mode == 'observed_response') 'response01' else 'integrated_benefit_score')), drop = FALSE])
    tmp <- tmp[cc, , drop = FALSE]
    if (nrow(tmp) < min_n) return(NULL)
    if (analysis_mode == 'observed_response') {
      if (length(unique(tmp$response01)) < 2) return(NULL)
      rhs <- paste(c('value', keep_cov), collapse = ' + ')
      fit <- tryCatch(glm(stats::as.formula(paste('response01 ~', rhs)), data = tmp, family = binomial()), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      co <- summary(fit)$coefficients
      if (!'value' %in% rownames(co)) return(NULL)
      tibble(
        cohort = unique(df$cohort)[1],
        feature = feat_name,
        model_type = 'adjusted_logistic',
        beta_adj = unname(co['value', 'Estimate']),
        se_adj = unname(co['value', 'Std. Error']),
        p_adj = unname(co['value', 'Pr(>|z|)']),
        n_total = nrow(tmp),
        n_case = sum(tmp$response01 == 1),
        n_control = sum(tmp$response01 == 0),
        covariates = paste(keep_cov, collapse = ';')
      )
    } else {
      rhs <- paste(c('value', keep_cov), collapse = ' + ')
      fit <- tryCatch(lm(stats::as.formula(paste('integrated_benefit_score ~', rhs)), data = tmp), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      co <- summary(fit)$coefficients
      if (!'value' %in% rownames(co)) return(NULL)
      tibble(
        cohort = unique(df$cohort)[1],
        feature = feat_name,
        model_type = 'adjusted_linear',
        beta_adj = unname(co['value', 'Estimate']),
        se_adj = unname(co['value', 'Std. Error']),
        p_adj = unname(co['value', 'Pr(>|t|)']),
        n_total = nrow(tmp),
        n_case = NA_real_,
        n_control = NA_real_,
        covariates = paste(keep_cov, collapse = ';')
      )
    }
  }))
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out %>%
    group_by(feature, model_type, covariates) %>%
    mutate(fdr_adj = p.adjust(p_adj, method = 'BH')) %>%
    ungroup()
}

validate_benefit_score_against_observed <- function(benefit_tbl, min_samples_per_group = 4) {
  if (is.null(benefit_tbl) || nrow(benefit_tbl) == 0) return(list(per_cohort = NULL, meta = NULL))
  val <- benefit_tbl %>% filter(!is.na(response01), is.finite(integrated_benefit_score))
  if (nrow(val) == 0) return(list(per_cohort = NULL, meta = NULL))
  per_cohort <- bind_rows(lapply(split(val, val$cohort), function(df) {
    if (nrow(df) < 2 * min_samples_per_group || length(unique(df$response01)) < 2) return(NULL)
    eff <- compute_binary_effects(df, value_col = 'integrated_benefit_score', response_col = 'response01', min_samples_per_group = min_samples_per_group)
    tibble(
      cohort = unique(df$cohort)[1],
      auc = as.numeric(eff[['auc']]),
      delta = as.numeric(eff[['delta']]),
      beta = as.numeric(eff[['beta']]),
      odds_ratio = as.numeric(eff[['odds_ratio']]),
      p_value = as.numeric(eff[['p_value']]),
      n_total = as.numeric(eff[['n_total']]),
      n_case = as.numeric(eff[['n_resp']]),
      n_control = as.numeric(eff[['n_nonresp']])
    )
  }))
  if (is.null(per_cohort) || nrow(per_cohort) == 0) return(list(per_cohort = NULL, meta = NULL))
  meta_eff <- approx_meta_heterogeneity(per_cohort$delta, per_cohort$n_case, per_cohort$n_control, per_cohort$n_total)
  meta <- tibble(
    metric = c('mean_auc', 'median_auc', 'weighted_mean_auc', 'meta_delta_fixed', 'meta_delta_random', 'meta_delta_random_p', 'heterogeneity_Q', 'heterogeneity_Q_p', 'heterogeneity_I2', 'tau2', 'n_validation_cohorts', 'n_validation_samples'),
    value = c(
      mean(per_cohort$auc, na.rm = TRUE),
      median(per_cohort$auc, na.rm = TRUE),
      weighted.mean(per_cohort$auc, w = pmax(1, per_cohort$n_total), na.rm = TRUE),
      meta_eff[['fixed_beta']],
      meta_eff[['random_beta']],
      meta_eff[['random_p']],
      meta_eff[['Q']],
      meta_eff[['Q_p']],
      meta_eff[['I2']],
      meta_eff[['tau2']],
      nrow(per_cohort),
      sum(per_cohort$n_total, na.rm = TRUE)
    )
  )
  list(per_cohort = per_cohort, meta = meta)
}

find_project_supported_terms <- function(paths) {
  project_terms <- character()
  for (pth in c(paths$enrich_hallmark, paths$enrich_integrated)) {
    dat <- read_flex(pth)
    if (!is.null(dat) && nrow(dat) > 0) {
      possible_cols <- names(dat)[tolower(names(dat)) %in% c("pathway", "term", "description", "gs_name", "name")]
      if (length(possible_cols) > 0) project_terms <- c(project_terms, as.character(dat[[possible_cols[1]]]))
    }
  }
  unique(clean_group_values(project_terms))
}

read_gene_table <- function(path) {
  df <- as_tibble(read_flex(path))
  if (is.null(df) || nrow(df) == 0) {
    stop("Uploaded gene file is missing or empty: ", path,
         " | checked basename: ", basename(path),
         " | working directory: ", getwd())
  }
  gcol <- find_gene_col(df)
  if (!"gene" %in% names(df)) df <- df %>% rename(gene = !!gcol)
  df <- df %>%
    mutate(gene = trimws(as.character(gene))) %>%
    filter(!is.na(gene), gene != "")
  if ("z_meta" %in% names(df) && !("p_meta" %in% names(df) && any(is.finite(df$p_meta)))) {
    df <- df %>% mutate(p_meta = 2 * pnorm(-abs(suppressWarnings(as.numeric(z_meta)))))
  }
  if (!"bag_frac" %in% names(df) && "bag_freq" %in% names(df) && "n_cohorts_used" %in% names(df)) {
    df <- df %>% mutate(bag_frac = suppressWarnings(as.numeric(bag_freq)) / pmax(1, suppressWarnings(as.numeric(n_cohorts_used))))
  }
  if (!"frac_cohorts_used" %in% names(df) && "n_cohorts_used" %in% names(df)) {
    mx <- max(suppressWarnings(as.numeric(df$n_cohorts_used)), na.rm = TRUE)
    if (is.finite(mx) && mx > 0) df <- df %>% mutate(frac_cohorts_used = suppressWarnings(as.numeric(n_cohorts_used)) / mx)
  }
  if (!"gene_priority_score" %in% names(df)) {
    z <- if ("z_meta" %in% names(df)) abs(suppressWarnings(as.numeric(df$z_meta))) else rep(0, nrow(df))
    bf <- if ("bag_frac" %in% names(df)) suppressWarnings(as.numeric(df$bag_frac)) else rep(0, nrow(df))
    fc <- if ("frac_cohorts_used" %in% names(df)) suppressWarnings(as.numeric(df$frac_cohorts_used)) else rep(0, nrow(df))
    p <- if ("p_meta" %in% names(df)) safe_neglog10(df$p_meta) else rep(0, nrow(df))
    df <- df %>% mutate(
      gene_priority_score =
        0.45 * safe_rescale01(bf) +
        0.20 * safe_rescale01(fc) +
        0.20 * safe_rescale01(z) +
        0.15 * safe_rescale01(replace(p, !is.finite(p), 0))
    )
  }
  df %>% distinct(gene, .keep_all = TRUE) %>% arrange(desc(gene_priority_score), desc(bag_frac), desc(frac_cohorts_used), gene)
}

read_support_gene_tbl <- function(path, source_label) {
  df <- read_flex(path)
  if (is.null(df) || nrow(df) == 0) return(tibble(gene = character(), source = character()))
  gcol <- find_gene_col(df)
  tibble(gene = trimws(as.character(df[[gcol]])), source = source_label) %>%
    filter(!is.na(gene), gene != "") %>% distinct()
}

# -----------------------------
# Read uploaded genes + contextual support
# -----------------------------
user_gene_tbl <- read_gene_table(gene_file)
if (nrow(user_gene_tbl) > top_gene_n) user_gene_tbl <- user_gene_tbl %>% slice_head(n = top_gene_n)
safe_write(user_gene_tbl, file.path(tab_dir, "table_uploaded_gene_list_ranked.tsv"))

ferr_driver <- read_support_gene_tbl(paths$ferr_driver, "FerrDb_driver")
ferr_suppressor <- read_support_gene_tbl(paths$ferr_suppressor, "FerrDb_suppressor")
ferr_marker <- read_support_gene_tbl(paths$ferr_marker, "FerrDb_marker")
ferr_all <- bind_rows(ferr_driver, ferr_suppressor, ferr_marker) %>% distinct()

locked_tbl <- read_support_gene_tbl(file.path(paths$out_v44, "genes", "locked_signature_rna_genes.csv"), "Locked_signature")
selected_tbl <- read_support_gene_tbl(file.path(paths$out_v44, "genes", "selected_rna_genes_fulltrain.csv"), "Selected_RNA")
wgcna_tbl <- bind_rows(
  read_support_gene_tbl(file.path(paths$wgcna_tables, "wgcna_genes_for_ML.csv"), "WGCNA_ML"),
  read_support_gene_tbl(file.path(paths$wgcna_tables, "conservative_ferroptosis_genes_GS_MM.tsv"), "WGCNA_conservative"),
  read_support_gene_tbl(file.path(paths$wgcna_tables, "module_hub_genes_MM_GS.tsv"), "WGCNA_hub")
) %>% distinct()

immune_like_pattern <- paste(
  c("^CD", "^CXCL", "^CCR", "^CTLA", "^PDCD", "^HLA", "^IL", "^IRF", "^JAK",
    "^TLR", "^AIM", "^GBP", "^BTK", "^BCL11B", "^PRKCB", "^NCKAP1L", "^DOCK8",
    "^PTPRC", "^FAS", "^FASLG", "^LCK", "^ZAP70", "^GZMB", "^PRF1", "^STAT1",
    "^IFI", "^IFIT", "^CXCR", "^CCL"),
  collapse = "|"
)

gene_context_tbl <- user_gene_tbl %>%
  mutate(
    gene = clean_names(gene),
    in_ferrdb = gene %in% ferr_all$gene,
    in_locked = gene %in% locked_tbl$gene,
    in_selected = gene %in% selected_tbl$gene,
    in_wgcna = gene %in% wgcna_tbl$gene,
    immune_context_flag = str_detect(gene, immune_like_pattern),
    support_count = as.integer(in_ferrdb) + as.integer(in_locked) + as.integer(in_selected) + as.integer(in_wgcna),
    context_layer = case_when(
      in_ferrdb & immune_context_flag ~ "Ferroptosis-immune interface",
      in_ferrdb & !immune_context_flag ~ "Ferroptosis-core supported",
      !in_ferrdb & immune_context_flag ~ "Immune-context dominant",
      TRUE ~ "Unresolved / exploratory"
    )
  ) %>%
  arrange(desc(gene_priority_score), desc(support_count), context_layer, gene)

safe_write(gene_context_tbl, file.path(tab_dir, "table_uploaded_genes_context_annotation.tsv"))

seed_genes <- unique(gene_context_tbl$gene)
ferro_core_genes <- gene_context_tbl %>% filter(context_layer %in% c("Ferroptosis-core supported", "Ferroptosis-immune interface")) %>% pull(gene) %>% unique()
immune_context_genes <- gene_context_tbl %>% filter(context_layer %in% c("Immune-context dominant", "Ferroptosis-immune interface")) %>% pull(gene) %>% unique()

# -----------------------------
# Cohort discovery and QC
# -----------------------------
if (!dir.exists(paths$geo_root)) stop("GEO_PREP directory not found under project root: ", paths$geo_root)
cohort_dirs <- list.dirs(paths$geo_root, recursive = FALSE, full.names = TRUE)
cohort_dirs <- cohort_dirs[basename(cohort_dirs) != "supp"]
if (length(cohort_dirs) == 0) stop("No cohort directories detected under GEO_PREP.")

cohort_data <- list()
qc_list <- list()
for (cohort_dir in cohort_dirs) {
  cohort_name <- basename(cohort_dir)
  expr_path <- first_existing(c(
    file.path(cohort_dir, "X_rna_symbol.csv"),
    file.path(cohort_dir, paste0(cohort_name, "_TPM.tsv.gz")),
    file.path(cohort_dir, paste0(cohort_name, "_TPM.tsv")),
    file.path(cohort_dir, paste0(cohort_name, "_expr.tsv.gz")),
    file.path(cohort_dir, paste0(cohort_name, "_expr.csv"))
  ))
  pheno_path <- first_existing(c(
    file.path(cohort_dir, paste0(cohort_name, "_pheno_expanded.csv")),
    file.path(cohort_dir, paste0(cohort_name, "_pheno_raw.csv")),
    file.path(cohort_dir, "pheno_expanded.csv"),
    file.path(cohort_dir, "pheno_raw.csv"),
    file.path(cohort_dir, "sample_map.csv")
  ))

  expr_df <- read_flex(expr_path)
  expr_mat <- parse_expression_table(expr_df)
  pheno_df <- read_flex(pheno_path)

  usable_expr <- !is.null(expr_mat) && nrow(expr_mat) >= 8 && ncol(expr_mat) >= 50
  usable <- FALSE
  reason <- "missing_or_invalid_expression"
  n_labeled <- 0L
  response_map <- NULL

  if (usable_expr) {
    sample_ids <- rownames(expr_mat)
    response_map <- infer_response_from_pheno(pheno_df)
    if (!is.null(response_map)) {
      common_ids <- intersect(sample_ids, response_map$sample_id)
      n_labeled <- sum(!is.na(response_map$response01[match(common_ids, response_map$sample_id)]))
      if (length(common_ids) >= min_labeled_samples &&
          length(unique(na.omit(response_map$response01[match(common_ids, response_map$sample_id)]))) == 2) {
        usable <- TRUE
        reason <- "observed_response"
      } else {
        reason <- "response_detected_but_underpowered"
      }
    } else {
      reason <- "response_not_detected"
    }

    cohort_data[[cohort_name]] <- list(
      cohort = cohort_name,
      expr = expr_mat,
      pheno = pheno_df,
      response_map = response_map,
      response_detected = usable,
      expr_path = expr_path,
      pheno_path = pheno_path
    )
  }

  qc_list[[cohort_name]] <- tibble(
    cohort = cohort_name,
    expr_path = ifelse(is.na(expr_path), "", expr_path),
    pheno_path = ifelse(is.na(pheno_path), "", pheno_path),
    n_samples_expr = if (is.null(expr_mat)) NA_integer_ else nrow(expr_mat),
    n_genes_expr = if (is.null(expr_mat)) NA_integer_ else ncol(expr_mat),
    response_detected = usable,
    n_labeled_response = n_labeled,
    reason = reason
  )
}

cohort_qc <- bind_rows(qc_list) %>% arrange(desc(response_detected), cohort)
safe_write(cohort_qc, file.path(tab_dir, "table_discovered_cohorts_qc.tsv"))
if (length(cohort_data) == 0) stop("No analyzable expression cohorts found under GEO_PREP.")

analysis_mode <- if (sum(cohort_qc$response_detected, na.rm = TRUE) > 0) "observed_response" else "response_propensity"
msg("[INFO] analysis mode: {analysis_mode}")

# -----------------------------
# Sample-level annotation
# -----------------------------
anno_all <- bind_rows(lapply(cohort_data, function(obj) {
  expr_mat <- obj$expr
  core_present <- intersect(ferro_core_genes, colnames(expr_mat))
  immune_present <- intersect(immune_context_genes, colnames(expr_mat))
  seed_present <- intersect(seed_genes, colnames(expr_mat))

  ferro_core_score_raw <- if (length(core_present) >= min_present_genes_per_score) rowMeans(expr_mat[, core_present, drop = FALSE], na.rm = TRUE) else rep(NA_real_, nrow(expr_mat))
  immune_context_score_raw <- if (length(immune_present) >= min_present_genes_per_score) rowMeans(expr_mat[, immune_present, drop = FALSE], na.rm = TRUE) else rep(NA_real_, nrow(expr_mat))
  seed_signature_score_raw <- if (length(seed_present) >= min_present_genes_per_score) rowMeans(expr_mat[, seed_present, drop = FALSE], na.rm = TRUE) else rep(NA_real_, nrow(expr_mat))

  seed_signature_score_z <- cohort_z(seed_signature_score_raw)
  risk_group <- ifelse(seed_signature_score_z >= median(seed_signature_score_z, na.rm = TRUE), "High", "Low")

  out <- tibble(
    sample_id = rownames(expr_mat),
    cohort = obj$cohort,
    n_seed_genes_present = length(seed_present),
    n_ferro_core_genes_present = length(core_present),
    n_immune_context_genes_present = length(immune_present),
    ferro_core_score_raw = ferro_core_score_raw,
    immune_context_score_raw = immune_context_score_raw,
    seed_signature_score_raw = seed_signature_score_raw,
    ferro_core_score_z = cohort_z(ferro_core_score_raw),
    immune_context_score_z = cohort_z(immune_context_score_raw),
    seed_signature_score_z = seed_signature_score_z,
    risk_group = risk_group,
    response01 = NA_integer_,
    response_label = NA_character_,
    response_source_col = NA_character_,
    response_mode = if (isTRUE(obj$response_detected)) "observed" else "latent"
  )
  if (!is.null(obj$response_map)) {
    hit <- match(out$sample_id, obj$response_map$sample_id)
    out$response01 <- obj$response_map$response01[hit]
    out$response_label <- obj$response_map$response_label[hit]
    out$response_source_col <- obj$response_map$response_source_col[hit]
    out$response_mode <- ifelse(!is.na(out$response01), "observed", out$response_mode)
  }
  out
}))
safe_write(anno_all, file.path(tab_dir, "table_sample_annotation.tsv"))

# -----------------------------
# MSigDB pathway universe and catalog selection
# -----------------------------
msig_h <- msigdbr::msigdbr(species = species_name, collection = "H")
msig_c7 <- msigdbr::msigdbr(species = species_name, collection = "C7")
msig_c8 <- msigdbr::msigdbr(species = species_name, collection = "C8")

msig_all <- bind_rows(
  msig_h %>% transmute(gs_name, gene_symbol, gs_description, collection = "Hallmark"),
  msig_c7 %>% transmute(gs_name, gene_symbol, gs_description, collection = "C7"),
  msig_c8 %>% transmute(gs_name, gene_symbol, gs_description, collection = "C8")
) %>%
  mutate(gene_symbol = clean_names(gene_symbol)) %>%
  distinct()

keyword_pattern <- paste(
  c("interferon", "ifn", "immune", "t cell", "cd8", "cytotoxic", "antigen",
    "checkpoint", "pd1", "pdl1", "ctla", "nk", "macrophage", "dendritic",
    "b cell", "myeloid", "strom", "fibroblast", "angiogenesis", "tgfb",
    "inflam", "chemokine", "exhaust", "hla", "mhc", "tumor microenvironment"),
  collapse = "|"
)

project_terms <- find_project_supported_terms(paths)

candidate_program_tbl <- msig_all %>%
  mutate(desc_low = tolower(paste(gs_name, gs_description))) %>%
  filter(str_detect(desc_low, keyword_pattern)) %>%
  group_by(gs_name, collection) %>%
  summarise(
    genes = list(unique(gene_symbol)),
    label = first(gs_description),
    n_genes = n_distinct(gene_symbol),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    overlap_n_seed = length(intersect(genes, seed_genes)),
    overlap_n_ferro = length(intersect(genes, ferro_core_genes)),
    overlap_n_immune = length(intersect(genes, immune_context_genes)),
    overlap_ratio_seed = overlap_n_seed / max(1, length(seed_genes)),
    project_supported = clean_group_values(gs_name) %in% project_terms | clean_group_values(label) %in% project_terms,
    pathway_family = pathway_family_label(gs_name, label),
    collection_priority = dplyr::case_when(collection == "Hallmark" ~ 3, collection == "C8" ~ 2, TRUE ~ 1),
    interpretability_penalty = ifelse(collection == "C7" & nchar(coalesce(label, gs_name)) > 120, 0.5, 0),
    pathway_score = overlap_n_seed + 0.5 * overlap_n_immune + 1.5 * overlap_n_ferro + 1.25 * collection_priority + ifelse(project_supported, 2, 0) - interpretability_penalty
  ) %>%
  ungroup() %>%
  arrange(desc(pathway_score), desc(project_supported), desc(overlap_n_seed), desc(collection_priority), pathway_family, gs_name) %>%
  prune_redundant_gene_sets(max_n = max_pathways, max_jaccard = 0.75) %>%
  mutate(
    feature = vapply(seq_len(n()), function(i) compact_pathway_label(gs_name[i], label[i], collection[i]), character(1)),
    feature_long = clean_feature_label(ifelse(is.na(label) | label == "", gs_name, label), max_len = 96)
  )

safe_write(candidate_program_tbl %>% select(gs_name, feature, feature_long, pathway_family, label, collection, n_genes, overlap_n_seed, overlap_n_ferro, overlap_n_immune, overlap_ratio_seed, project_supported, pathway_score),
           file.path(tab_dir, "table_selected_immunotherapy_pathway_catalog.tsv"))

program_gene_sets <- candidate_program_tbl$genes
names(program_gene_sets) <- candidate_program_tbl$gs_name
program_labels <- setNames(candidate_program_tbl$feature, candidate_program_tbl$gs_name)
if (length(program_gene_sets) == 0) stop("No package-derived immunotherapy pathways were selected.")

# -----------------------------
# Pathway scoring
# -----------------------------
program_long <- bind_rows(lapply(cohort_data, function(obj) {
  expr_use <- t(obj$expr)
  expr_use <- cohort_row_z(expr_use)
  available_sets <- lapply(program_gene_sets, function(gs) intersect(gs, rownames(expr_use)))
  available_sets <- available_sets[vapply(available_sets, length, integer(1)) >= min_gene_set]
  if (length(available_sets) < 3) return(NULL)
  scores <- safe_gsva_ssgsea(expr_use, available_sets, min_gene_set = min_gene_set)
  if (is.null(scores) || nrow(scores) == 0) return(NULL)
  as_tibble(as.data.frame(t(scores)), rownames = "sample_id") %>%
    pivot_longer(cols = -sample_id, names_to = "gs_name", values_to = "score_raw") %>%
    mutate(
      cohort = obj$cohort,
      feature = ifelse(gs_name %in% names(program_labels), unname(program_labels[gs_name]), gs_name),
      feature_source = "MSigDB_ssGSEA"
    ) %>%
    group_by(cohort, gs_name) %>%
    mutate(score_z = cohort_z(score_raw), score_rank01 = cohort_rank01(score_raw)) %>%
    ungroup() %>%
    left_join(anno_all, by = c("sample_id", "cohort"))
}))
if (is.null(program_long) || nrow(program_long) == 0) stop("No package-derived pathways could be scored.")
safe_write(program_long, file.path(tab_dir, "table_pathway_scores_by_sample_long.tsv"))

# -----------------------------
# Optional IOBR deconvolution
# -----------------------------
deconv_long <- NULL
if (has_iobr) {
  msg("[INFO] IOBR detected; attempting deconvolution")
  deconv_long <- bind_rows(lapply(cohort_data, function(obj) {
    expr_use <- t(obj$expr)
    expr_use <- as.matrix(expr_use)
    mode(expr_use) <- "numeric"
    expr_use <- expr_use[rowSums(is.finite(expr_use)) > 0, , drop = FALSE]
    if (nrow(expr_use) < 100 || ncol(expr_use) < 8) return(NULL)
    methods_try <- c("estimate", "mcpcounter", "xcell")
    all_long <- list()
    idx <- 1L
    for (m in methods_try) {
      tmp <- tryCatch({
        out <- IOBR::deconvo_tme(eset = expr_use, method = m, arrays = TRUE)
        out <- as.data.frame(out)
        if (!"ID" %in% names(out)) out$ID <- rownames(out)
        as_tibble(out) %>%
          rename(sample_id = ID) %>%
          pivot_longer(cols = -sample_id, names_to = "feature", values_to = "score_raw") %>%
          mutate(cohort = obj$cohort, feature_source = paste0("IOBR_", m)) %>%
          group_by(cohort, feature) %>%
          mutate(score_z = cohort_z(score_raw), score_rank01 = cohort_rank01(score_raw)) %>%
          ungroup() %>%
          left_join(anno_all, by = c("sample_id", "cohort"))
      }, error = function(e) NULL)
      all_long[[idx]] <- tmp
      idx <- idx + 1L
    }
    bind_rows(all_long)
  }))
  if (!is.null(deconv_long) && nrow(deconv_long) > 0) {
    safe_write(deconv_long, file.path(tab_dir, "table_deconvolution_scores_by_sample_long.tsv"))
  }
}

# -----------------------------
# Integrated benefit propensity
# -----------------------------
benefit_tbl <- program_long %>%
  group_by(sample_id, cohort) %>%
  summarise(io_pathway_score = mean(score_z[is.finite(score_z)], na.rm = TRUE), .groups = "drop") %>%
  left_join(anno_all, by = c("sample_id", "cohort")) %>%
  mutate(
    io_pathway_score_z = cohort_z(io_pathway_score),
    immune_context_score_z = cohort_z(immune_context_score_raw),
    seed_signature_score_z = cohort_z(seed_signature_score_raw),
    integrated_benefit_score = 0.60 * io_pathway_score_z + 0.25 * immune_context_score_z - 0.15 * seed_signature_score_z,
    integrated_benefit_score = cohort_z(integrated_benefit_score),
    integrated_benefit_group = ifelse(integrated_benefit_score >= median(integrated_benefit_score, na.rm = TRUE), "Benefit-High", "Benefit-Low")
  )
benefit_tbl$integrated_benefit_score[!is.finite(benefit_tbl$integrated_benefit_score)] <- NA_real_
safe_write(benefit_tbl, file.path(tab_dir, "table_integrated_benefit_score.tsv"))

anno_all <- anno_all %>% left_join(benefit_tbl %>% select(sample_id, cohort, io_pathway_score, io_pathway_score_z, integrated_benefit_score, integrated_benefit_group), by = c("sample_id", "cohort"))
program_long <- program_long %>% left_join(benefit_tbl %>% select(sample_id, cohort, integrated_benefit_score, integrated_benefit_group), by = c("sample_id", "cohort"))
if (!is.null(deconv_long) && nrow(deconv_long) > 0) {
  deconv_long <- deconv_long %>% left_join(benefit_tbl %>% select(sample_id, cohort, integrated_benefit_score, integrated_benefit_group), by = c("sample_id", "cohort"))
}

analysis_group_col <- if (analysis_mode == "observed_response") "response01" else "integrated_benefit_group"

# -----------------------------
# Grouped statistics
# -----------------------------
compute_grouped_effects <- function(long_df, group_col, source_label, value_col = "score_z", feature_col = "feature") {
  if (is.null(long_df) || nrow(long_df) == 0) return(NULL)
  pieces <- split(long_df, list(long_df$cohort, long_df[[feature_col]]), drop = TRUE)
  bind_rows(lapply(pieces, function(df) {
    if (nrow(df) == 0) return(NULL)
    if (group_col == "response01") {
      eff <- compute_binary_effects(df, value_col = value_col, response_col = "response01", min_samples_per_group = min_samples_per_group)
      tibble(
        cohort = unique(df$cohort)[1],
        feature = unique(df[[feature_col]])[1],
        feature_source = source_label,
        p_value = as.numeric(eff[["p_value"]]),
        delta = as.numeric(eff[["delta"]]),
        auc = as.numeric(eff[["auc"]]),
        odds_ratio = as.numeric(eff[["odds_ratio"]]),
        beta = as.numeric(eff[["beta"]]),
        rank_biserial = as.numeric(eff[["rank_biserial"]]),
        n_total = as.numeric(eff[["n_total"]]),
        n_case = as.numeric(eff[["n_resp"]]),
        n_control = as.numeric(eff[["n_nonresp"]])
      )
    } else {
      grp_vals <- unique(as.character(df[[group_col]]))
      preferred_group1 <- if ("Benefit-High" %in% grp_vals) "Benefit-High" else sort(grp_vals)[2]
      preferred_group0 <- if ("Benefit-Low" %in% grp_vals) "Benefit-Low" else sort(setdiff(grp_vals, preferred_group1))[1]
      eff <- compute_two_group_effects(
        df[[value_col]], df[[group_col]],
        min_samples_per_group = min_samples_per_group,
        group1 = preferred_group1,
        group0 = preferred_group0
      )
      tibble(
        cohort = unique(df$cohort)[1],
        feature = unique(df[[feature_col]])[1],
        feature_source = source_label,
        comparison = paste0(preferred_group1, "_minus_", preferred_group0),
        p_value = as.numeric(eff[["p_value"]]),
        delta = as.numeric(eff[["delta"]]),
        auc = NA_real_,
        odds_ratio = NA_real_,
        beta = NA_real_,
        rank_biserial = as.numeric(eff[["rank_biserial"]]),
        n_total = as.numeric(eff[["n_total"]]),
        n_case = as.numeric(eff[["n_group1"]]),
        n_control = as.numeric(eff[["n_group0"]])
      )
    }
  }))
}

program_effects <- compute_grouped_effects(program_long, analysis_group_col, "MSigDB_ssGSEA", value_col = "score_z", feature_col = "feature")
program_effects <- program_effects %>%
  mutate(
    fdr_global = p.adjust(p_value, method = "BH"),
    effect_weight = sqrt(pmax(1, n_total)),
    effect_sign = sign(delta)
  )
safe_write(program_effects, file.path(tab_dir, "table_pathway_effects_by_cohort.tsv"))

program_meta <- program_effects %>%
  group_by(feature, feature_source) %>%
  summarise(
    n_cohorts = sum(is.finite(delta) | is.finite(auc)),
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    sd_delta = sd(delta, na.rm = TRUE),
    mean_auc = mean(auc, na.rm = TRUE),
    best_auc = safe_max_num(auc),
    min_fdr_global = safe_min_num(fdr_global),
    consistency_fraction = mean(sign(delta[is.finite(delta)]) == sign(mean(delta, na.rm = TRUE)), na.rm = TRUE),
    meta_p = meta_signed_stouffer(p_value, delta, effect_weight)[["meta_p"]],
    meta_z = meta_signed_stouffer(p_value, delta, effect_weight)[["meta_z"]],
    .groups = "drop"
  ) %>%
  mutate(
    meta_fdr = p.adjust(meta_p, method = "BH"),
    meta_rank = safe_neglog10(meta_fdr) + abs(mean_delta) + ifelse(is.finite(mean_auc), pmax(0, mean_auc - 0.5), 0) + 0.5 * coalesce(consistency_fraction, 0)
  ) %>%
  arrange(meta_fdr, desc(meta_rank), desc(abs(mean_delta)))
safe_write(program_meta, file.path(tab_dir, "table_pathway_meta_summary.tsv"))

program_loo <- loo_meta_rank(program_effects, feature_col = "feature", effect_col = "delta", p_col = "p_value", weight_col = "effect_weight")
if (!is.null(program_loo) && nrow(program_loo) > 0) {
  safe_write(program_loo, file.path(tab_dir, "table_pathway_leave_one_cohort_out.tsv"))
}

gene_level_long <- bind_rows(lapply(cohort_data, function(obj) {
  expr <- as_tibble(obj$expr, rownames = "sample_id")
  present <- intersect(seed_genes, names(expr))
  if (length(present) < 3) return(NULL)
  expr %>%
    select(sample_id, all_of(present)) %>%
    pivot_longer(cols = -sample_id, names_to = "gene", values_to = "score_raw") %>%
    mutate(cohort = obj$cohort) %>%
    group_by(cohort, gene) %>%
    mutate(score_z = cohort_z(score_raw), score_rank01 = cohort_rank01(score_raw)) %>%
    ungroup() %>%
    left_join(anno_all, by = c("sample_id", "cohort")) %>%
    left_join(gene_context_tbl %>% select(gene, context_layer, gene_priority_score, support_count), by = "gene")
}))
safe_write(gene_level_long, file.path(tab_dir, "table_gene_scores_by_sample_long.tsv"))

gene_effects <- compute_grouped_effects(gene_level_long, analysis_group_col, "Uploaded_gene", value_col = "score_z", feature_col = "gene") %>%
  left_join(gene_context_tbl %>% select(gene, context_layer, gene_priority_score, support_count), by = c("feature" = "gene")) %>%
  rename(gene = feature) %>%
  mutate(
    fdr_global = p.adjust(p_value, method = "BH"),
    effect_weight = sqrt(pmax(1, n_total))
  )
safe_write(gene_effects, file.path(tab_dir, "table_gene_effects_by_cohort.tsv"))

gene_meta <- gene_effects %>%
  group_by(gene, context_layer, gene_priority_score, support_count) %>%
  summarise(
    n_cohorts = sum(is.finite(delta) | is.finite(auc)),
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    sd_delta = sd(delta, na.rm = TRUE),
    mean_auc = mean(auc, na.rm = TRUE),
    min_fdr_global = safe_min_num(fdr_global),
    consistency_fraction = mean(sign(delta[is.finite(delta)]) == sign(mean(delta, na.rm = TRUE)), na.rm = TRUE),
    meta_p = meta_signed_stouffer(p_value, delta, effect_weight)[["meta_p"]],
    meta_z = meta_signed_stouffer(p_value, delta, effect_weight)[["meta_z"]],
    .groups = "drop"
  ) %>%
  mutate(
    meta_fdr = p.adjust(meta_p, method = "BH"),
    meta_rank = safe_neglog10(meta_fdr) + 0.25 * safe_rescale01(gene_priority_score) + 0.15 * support_count + abs(mean_delta) + ifelse(is.finite(mean_auc), pmax(0, mean_auc - 0.5), 0) + 0.5 * coalesce(consistency_fraction, 0)
  ) %>%
  arrange(meta_fdr, desc(meta_rank), desc(abs(mean_delta)))
safe_write(gene_meta, file.path(tab_dir, "table_gene_meta_summary.tsv"))

gene_loo <- loo_meta_rank(gene_effects %>% rename(feature = gene), feature_col = "feature", effect_col = "delta", p_col = "p_value", weight_col = "effect_weight")
if (!is.null(gene_loo) && nrow(gene_loo) > 0) {
  safe_write(gene_loo, file.path(tab_dir, "table_gene_leave_one_cohort_out.tsv"))
}

deconv_effects <- NULL
deconv_meta <- NULL
if (!is.null(deconv_long) && nrow(deconv_long) > 0) {
  deconv_effects <- compute_grouped_effects(deconv_long, analysis_group_col, "IOBR_deconvolution", value_col = "score_z", feature_col = "feature") %>%
    mutate(
      fdr_global = p.adjust(p_value, method = "BH"),
      effect_weight = sqrt(pmax(1, n_total))
    )
  safe_write(deconv_effects, file.path(tab_dir, "table_deconvolution_effects_by_cohort.tsv"))
  deconv_meta <- deconv_effects %>%
    group_by(feature, feature_source) %>%
    summarise(
      n_cohorts = sum(is.finite(delta) | is.finite(auc)),
      mean_delta = mean(delta, na.rm = TRUE),
      mean_auc = mean(auc, na.rm = TRUE),
      min_fdr_global = safe_min_num(fdr_global),
      meta_p = meta_signed_stouffer(p_value, delta, effect_weight)[["meta_p"]],
      consistency_fraction = mean(sign(delta[is.finite(delta)]) == sign(mean(delta, na.rm = TRUE)), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      meta_fdr = p.adjust(meta_p, method = "BH"),
      meta_rank = safe_neglog10(meta_fdr) + abs(mean_delta) + 0.5 * coalesce(consistency_fraction, 0)
    ) %>%
    arrange(meta_fdr, desc(meta_rank), desc(abs(mean_delta)))
  safe_write(deconv_meta, file.path(tab_dir, "table_deconvolution_meta_summary.tsv"))
}



# -----------------------------
# Q1-strengthening validation / heterogeneity / confounder adjustment
# -----------------------------
validation_res <- validate_benefit_score_against_observed(benefit_tbl, min_samples_per_group = min_samples_per_group)
validation_tbl <- validation_res$per_cohort
validation_meta <- validation_res$meta
if (!is.null(validation_tbl) && nrow(validation_tbl) > 0) {
  safe_write(validation_tbl, file.path(tab_dir, "table_observed_response_validation_by_cohort.tsv"))
}
if (!is.null(validation_meta) && nrow(validation_meta) > 0) {
  safe_write(validation_meta, file.path(tab_dir, "table_observed_response_validation_meta.tsv"))
}

program_re_meta <- program_effects %>%
  group_by(feature, feature_source) %>%
  summarise(
    k = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['k']],
    fixed_beta = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['fixed_beta']],
    fixed_se = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['fixed_se']],
    fixed_p = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['fixed_p']],
    random_beta = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['random_beta']],
    random_se = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['random_se']],
    random_p = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['random_p']],
    Q = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['Q']],
    Q_p = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['Q_p']],
    I2 = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['I2']],
    tau2 = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['tau2']],
    .groups = 'drop'
  ) %>%
  mutate(random_fdr = p.adjust(random_p, method = 'BH')) %>%
  arrange(random_fdr, desc(abs(random_beta)))
safe_write(program_re_meta, file.path(tab_dir, 'table_pathway_random_effects_meta.tsv'))

gene_re_meta <- gene_effects %>%
  group_by(gene, context_layer, gene_priority_score, support_count) %>%
  summarise(
    k = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['k']],
    fixed_beta = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['fixed_beta']],
    fixed_se = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['fixed_se']],
    fixed_p = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['fixed_p']],
    random_beta = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['random_beta']],
    random_se = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['random_se']],
    random_p = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['random_p']],
    Q = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['Q']],
    Q_p = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['Q_p']],
    I2 = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['I2']],
    tau2 = approx_meta_heterogeneity(delta, n_case, n_control, n_total)[['tau2']],
    .groups = 'drop'
  ) %>%
  mutate(random_fdr = p.adjust(random_p, method = 'BH')) %>%
  arrange(random_fdr, desc(abs(random_beta)))
safe_write(gene_re_meta, file.path(tab_dir, 'table_gene_random_effects_meta.tsv'))

program_top_for_adjustment <- head(program_meta$feature, 15)
gene_top_for_adjustment <- head(gene_meta$gene, 20)
program_adjusted <- compute_adjusted_feature_associations(program_long, deconv_long, analysis_mode = analysis_mode, feature_col = 'feature', value_col = 'score_z', feature_subset = program_top_for_adjustment, max_covars = 4, min_n = 20)
gene_adjusted <- compute_adjusted_feature_associations(gene_level_long %>% rename(feature = gene), deconv_long, analysis_mode = analysis_mode, feature_col = 'feature', value_col = 'score_z', feature_subset = gene_top_for_adjustment, max_covars = 4, min_n = 20)
if (!is.null(program_adjusted) && nrow(program_adjusted) > 0) {
  safe_write(program_adjusted, file.path(tab_dir, 'table_pathway_adjusted_associations.tsv'))
  program_adjusted_meta <- program_adjusted %>%
    group_by(feature, model_type, covariates) %>%
    summarise(
      n_cohorts = n(),
      mean_beta_adj = mean(beta_adj, na.rm = TRUE),
      consistency_fraction = mean(sign(beta_adj[is.finite(beta_adj)]) == sign(mean(beta_adj, na.rm = TRUE)), na.rm = TRUE),
      meta_p = meta_signed_stouffer(p_adj, beta_adj, sqrt(pmax(1, n_total)))[['meta_p']],
      .groups = 'drop'
    ) %>%
    mutate(meta_fdr = p.adjust(meta_p, method = 'BH')) %>%
    arrange(meta_fdr, desc(abs(mean_beta_adj)))
  safe_write(program_adjusted_meta, file.path(tab_dir, 'table_pathway_adjusted_meta.tsv'))
} else {
  program_adjusted_meta <- NULL
}
if (!is.null(gene_adjusted) && nrow(gene_adjusted) > 0) {
  safe_write(gene_adjusted, file.path(tab_dir, 'table_gene_adjusted_associations.tsv'))
  gene_adjusted_meta <- gene_adjusted %>%
    group_by(feature, model_type, covariates) %>%
    summarise(
      n_cohorts = n(),
      mean_beta_adj = mean(beta_adj, na.rm = TRUE),
      consistency_fraction = mean(sign(beta_adj[is.finite(beta_adj)]) == sign(mean(beta_adj, na.rm = TRUE)), na.rm = TRUE),
      meta_p = meta_signed_stouffer(p_adj, beta_adj, sqrt(pmax(1, n_total)))[['meta_p']],
      .groups = 'drop'
    ) %>%
    mutate(meta_fdr = p.adjust(meta_p, method = 'BH')) %>%
    arrange(meta_fdr, desc(abs(mean_beta_adj)))
  safe_write(gene_adjusted_meta, file.path(tab_dir, 'table_gene_adjusted_meta.tsv'))
} else {
  gene_adjusted_meta <- NULL
}


# -----------------------------
# Figures
# -----------------------------
plot_df <- benefit_tbl %>%
  mutate(display_group = if (analysis_mode == "observed_response") response_label else integrated_benefit_group) %>%
  filter(is.finite(integrated_benefit_score))

if (nrow(plot_df) > 0) {
  if (analysis_mode == "observed_response") {
    p <- ggplot(plot_df %>% filter(!is.na(display_group)), aes(x = cohort, y = integrated_benefit_score, fill = display_group)) +
      geom_violin(scale = "width", trim = FALSE, alpha = 0.75) +
      geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9) +
      coord_flip() +
      labs(title = paste0("Integrated immunotherapy benefit score across cohorts (", analysis_mode, ")"),
           x = NULL, y = "Integrated benefit score (cohort-z)") +
      theme_bw(base_size = 12)
  } else {
    cohort_medians <- plot_df %>% group_by(cohort) %>% summarise(median_score = median(integrated_benefit_score, na.rm = TRUE), .groups = "drop")
    p <- ggplot(plot_df, aes(x = cohort, y = integrated_benefit_score)) +
      geom_violin(fill = "grey80", color = "grey25", scale = "width", trim = FALSE, alpha = 0.9) +
      geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "grey20") +
      geom_point(data = cohort_medians, aes(x = cohort, y = median_score), inherit.aes = FALSE, size = 2.5, color = "firebrick") +
      coord_flip() +
      labs(title = "Integrated immunotherapy benefit score across cohorts (latent propensity only)",
           subtitle = "Benefit-High vs Benefit-Low split is not plotted here to avoid circular internal validation.",
           x = NULL, y = "Integrated benefit score (cohort-z)") +
      theme_bw(base_size = 12)
  }
  plot_save(p, file.path(fig_dir, "FIG01_integrated_benefit_score_by_group"), 11, 7)
}

heat_df <- program_long %>%
  mutate(display_group = if (analysis_mode == "observed_response") response_label else integrated_benefit_group) %>%
  filter(!is.na(display_group), is.finite(score_z)) %>%
  group_by(cohort, display_group, feature) %>%
  summarise(mean_score = mean(score_z, na.rm = TRUE), .groups = "drop")

if (nrow(heat_df) > 0) {
  mat <- heat_df %>%
    unite("cohort_group", cohort, display_group, sep = "__") %>%
    pivot_wider(names_from = cohort_group, values_from = mean_score) %>%
    as.data.frame()
  rownames(mat) <- mat$feature
  mat$feature <- NULL
  if (nrow(mat) > 1 && ncol(mat) > 1) {
    safe_pheatmap_save(
      as.matrix(mat),
      file.path(fig_dir, "FIG02_immunotherapy_pathway_heatmap"),
      width = 2400, height = 1700, res = 220,
      scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
      fontsize_row = 9, fontsize_col = 10,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      main = paste0("Redundancy-pruned immunotherapy pathways (", analysis_mode, ")")
    )
  }
}

meta_plot_tbl <- top_n_rows(program_meta, n = 20, order_col = if (analysis_mode == "observed_response") "mean_auc" else "mean_delta", decreasing = TRUE) %>%
  mutate(feature = safe_factor_reorder(feature, if (analysis_mode == "observed_response") mean_auc else mean_delta))
if (nrow(meta_plot_tbl) > 0) {
  xvar <- if (analysis_mode == "observed_response") "mean_auc" else "mean_delta"
  xlab <- if (analysis_mode == "observed_response") "Mean AUC" else "Mean delta (cohort-z)"
  xint <- if (analysis_mode == "observed_response") 0.5 else 0
  p <- ggplot(meta_plot_tbl, aes(x = .data[[xvar]], y = feature, color = safe_neglog10(meta_fdr), size = n_cohorts)) +
    geom_vline(xintercept = xint, linetype = 2, color = "grey60") +
    geom_point() +
    scale_color_viridis_c(option = "C", na.value = "grey70") +
    labs(title = "Top immunotherapy pathways", x = xlab, y = NULL, color = "-log10 meta-FDR", size = "Cohorts") +
    theme_bw(base_size = 12)
  plot_save(p, file.path(fig_dir, "FIG03_top_pathway_meta_summary"), 10, 7)
}

gene_plot_tbl <- top_n_rows(gene_meta, n = 20, order_col = if (analysis_mode == "observed_response") "mean_auc" else "mean_delta", decreasing = TRUE) %>%
  mutate(gene = safe_factor_reorder(gene, if (analysis_mode == "observed_response") mean_auc else mean_delta))
if (nrow(gene_plot_tbl) > 0) {
  xvar <- if (analysis_mode == "observed_response") "mean_auc" else "mean_delta"
  xlab <- if (analysis_mode == "observed_response") "Mean AUC" else "Mean delta (cohort-z)"
  xint <- if (analysis_mode == "observed_response") 0.5 else 0
  p <- ggplot(gene_plot_tbl, aes(x = .data[[xvar]], y = gene, color = context_layer, size = support_count)) +
    geom_vline(xintercept = xint, linetype = 2, color = "grey60") +
    geom_point(alpha = 0.9) +
    labs(title = "Top uploaded genes in immunotherapy association analysis", x = xlab, y = NULL, color = "Context", size = "Support count") +
    theme_bw(base_size = 12)
  plot_save(p, file.path(fig_dir, "FIG04_top_gene_meta_summary"), 11, 7)
}

context_bar <- gene_context_tbl %>% count(context_layer, sort = TRUE)
if (nrow(context_bar) > 0) {
  p <- ggplot(context_bar, aes(x = reorder(context_layer, n), y = n, fill = context_layer)) +
    geom_col() + coord_flip() +
    labs(title = "Uploaded gene-list context composition", x = NULL, y = "Number of genes") +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  plot_save(p, file.path(fig_dir, "FIG05_gene_context_composition"), 8, 5)
}

consistency_tbl <- top_n_rows(program_meta, n = 25, order_col = "consistency_fraction", decreasing = TRUE) %>%
  mutate(feature = safe_factor_reorder(feature, consistency_fraction))
if (nrow(consistency_tbl) > 0) {
  p <- ggplot(consistency_tbl, aes(x = consistency_fraction, y = feature, color = safe_neglog10(meta_fdr), size = n_cohorts)) +
    geom_vline(xintercept = 0.5, linetype = 2, color = "grey60") +
    geom_point() +
    scale_color_viridis_c(option = "C", na.value = "grey70") +
    labs(title = "Cross-cohort directional consistency of pathways", x = "Consistency fraction", y = NULL, color = "-log10 meta-FDR", size = "Cohorts") +
    theme_bw(base_size = 12)
  plot_save(p, file.path(fig_dir, "FIG06_pathway_directional_consistency"), 10, 7)
}



if (!is.null(validation_tbl) && nrow(validation_tbl) > 0) {
  p <- ggplot(validation_tbl, aes(x = reorder(cohort, auc), y = auc)) +
    geom_col(fill = '#4C78A8') +
    geom_hline(yintercept = 0.5, linetype = 2, color = 'grey50') +
    coord_flip() +
    labs(title = 'Observed-response benchmarking of integrated benefit score', x = NULL, y = 'AUC vs observed response') +
    theme_bw(base_size = 12)
  plot_save(p, file.path(fig_dir, 'FIG07_observed_response_validation_auc'), 9, 5.5)
}

if (nrow(program_re_meta) > 0) {
  het_plot <- program_re_meta %>%
    slice_head(n = 15) %>%
    mutate(feature = safe_factor_reorder(feature, random_beta, desc = TRUE))
  p <- ggplot(het_plot, aes(x = random_beta, y = feature, size = k, color = I2)) +
    geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
    geom_point() +
    scale_color_viridis_c(option = 'C', na.value = 'grey70') +
    labs(title = 'Random-effects pathway meta-analysis', x = 'Random-effects beta', y = NULL, color = 'I2', size = 'Cohorts') +
    theme_bw(base_size = 12)
  plot_save(p, file.path(fig_dir, 'FIG08_pathway_random_effects_meta'), 10, 7)
}

if (exists('program_adjusted_meta') && !is.null(program_adjusted_meta) && nrow(program_adjusted_meta) > 0) {
  adj_plot <- program_adjusted_meta %>%
    slice_head(n = 12) %>%
    mutate(feature = safe_factor_reorder(feature, mean_beta_adj, desc = TRUE))
  p <- ggplot(adj_plot, aes(x = mean_beta_adj, y = feature, size = n_cohorts, color = safe_neglog10(meta_fdr))) +
    geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
    geom_point() +
    scale_color_viridis_c(option = 'C', na.value = 'grey70') +
    labs(title = 'Covariate-adjusted pathway associations', x = 'Adjusted mean beta', y = NULL, color = '-log10 adj meta-FDR', size = 'Cohorts') +
    theme_bw(base_size = 12)
  plot_save(p, file.path(fig_dir, 'FIG09_adjusted_pathway_associations'), 10, 7)
}


# -----------------------------
# Reporting
# -----------------------------
result_highlights <- tibble(
  metric = c(
    "n_uploaded_genes",
    "n_ferroptosis_core_supported",
    "n_immune_context_dominant",
    "n_ferroptosis_immune_interface",
    "n_discovered_expression_cohorts",
    "n_cohorts_with_observed_response_labels",
    "analysis_mode",
    "n_selected_pathways",
    "n_pathway_features_scored",
    "n_deconvolution_features_scored",
    "median_samples_per_cohort",
    "iobr_available",
    "top_pathway_feature",
    "top_gene_feature"
  ),
  value = c(
    nrow(user_gene_tbl),
    sum(gene_context_tbl$context_layer == "Ferroptosis-core supported", na.rm = TRUE),
    sum(gene_context_tbl$context_layer == "Immune-context dominant", na.rm = TRUE),
    sum(gene_context_tbl$context_layer == "Ferroptosis-immune interface", na.rm = TRUE),
    length(cohort_data),
    sum(cohort_qc$response_detected, na.rm = TRUE),
    analysis_mode,
    length(program_gene_sets),
    dplyr::n_distinct(program_long$feature),
    ifelse(is.null(deconv_long), 0, dplyr::n_distinct(deconv_long$feature)),
    median(vapply(cohort_data, function(obj) nrow(obj$expr), numeric(1)), na.rm = TRUE),
    has_iobr,
    ifelse(nrow(program_meta) > 0, program_meta$feature[1], NA_character_),
    ifelse(nrow(gene_meta) > 0, gene_meta$gene[1], NA_character_)
  )
)
safe_write(result_highlights, file.path(tab_dir, "table_result_highlights.tsv"))

manuscript_summary <- tibble(
  section = c(
    "Design",
    "Primary endpoint mode",
    "Gene-context framing",
    "Pathway strategy",
    "Meta-analysis",
    "Sensitivity analysis",
    "Optional microenvironment deconvolution"
  ),
  detail = c(
    "Cross-cohort expression analysis with automatic cohort discovery under GEO_PREP.",
    ifelse(analysis_mode == "observed_response", "Observed responder/non-responder labels were available in at least one cohort and used as the primary endpoint where detected.", "Observed response labels were insufficient; analysis therefore used an integrated immunotherapy benefit propensity score as the latent endpoint."),
    "Uploaded genes were partitioned into ferroptosis-core, immune-context dominant, and ferroptosis-immune interface layers to mitigate reviewer concerns about biologic confounding.",
    "Immune-relevant pathway programs were selected from MSigDB Hallmark/C7/C8 using overlap with the uploaded seed genes plus project-supported enrichment terms, then scored by ssGSEA.",
    "Feature-level synthesis uses per-cohort effect estimation followed by direction-aware weighted Stouffer meta-analysis and BH correction.",
    "Leave-one-cohort-out summaries quantify robustness of pathway and gene rankings to cohort exclusion.",
    ifelse(has_iobr, "IOBR-based deconvolution was attempted using estimate, MCPcounter, and xCell backends when data were suitable.", "IOBR was not installed; deconvolution outputs were therefore omitted without affecting the main analysis." )
  )
)
safe_write(manuscript_summary, file.path(tab_dir, "table_manuscript_summary.tsv"))

run_metadata <- tibble(
  key = c(
    "project_root", "gene_file", "script_name", "analysis_mode", "generated_at", "species_name",
    "install_missing", "min_gene_set", "max_pathways", "min_samples_per_group",
    "min_labeled_samples", "top_gene_n", "min_present_genes_per_score", "seed", "has_iobr", "has_pROC"
  ),
  value = c(
    root, normalizePath(gene_file, winslash = "/", mustWork = FALSE), script_basename, analysis_mode,
    as.character(Sys.time()), species_name, as.character(install_missing), as.character(min_gene_set),
    as.character(max_pathways), as.character(min_samples_per_group), as.character(min_labeled_samples),
    as.character(top_gene_n), as.character(min_present_genes_per_score), as.character(seed),
    as.character(has_iobr), as.character(has_pROC)
  )
)
safe_write(run_metadata, file.path(log_dir, "run_metadata.tsv"))

manifest <- tibble(path = list.files(out_dir, recursive = TRUE, full.names = FALSE))
safe_write(manifest, file.path(log_dir, "output_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

summary_lines <- c(
  "Immunotherapy Response Prediction / Propensity Analysis Summary",
  paste0("Project root: ", root),
  paste0("Gene file: ", normalizePath(gene_file, winslash = "/", mustWork = FALSE)),
  paste0("Generated at: ", Sys.time()),
  paste0("Analysis mode: ", analysis_mode),
  "",
  "Q1-oriented strengthening moves included in this script:",
  "- Missing p_meta is reconstructed from z_meta when required.",
  "- Uploaded genes are stratified into ferroptosis-core, immune-context, and interface layers.",
  "- Pathway scoring uses package-derived MSigDB H/C7/C8 sets with project-term support.",
  "- Cohort-wise z-standardization is applied before integrated synthesis to reduce scale heterogeneity.",
  "- In latent propensity mode, internal Benefit-High/Benefit-Low splits are not used for headline score-distribution plots to avoid circular visual validation.",
  "- Redundancy-pruned pathway catalogs and compact pathway labels improve reviewer-facing interpretability.",
  "- Primary results include per-cohort effects, signed weighted meta-analysis, BH control, and leave-one-cohort-out robustness.",
  "- Optional IOBR deconvolution is integrated when installed and data are suitable.
  - Observed-response benchmarking is exported whenever real responder labels are available.
  - Random-effects heterogeneity summaries (Q, I2, tau2) are reported for pathway and gene meta-analysis.
  - Covariate-adjusted association tables are generated when deconvolution-derived TME surrogates are available.",
  "- All outputs are written into a script-named project-root folder for manuscript traceability."
)
writeLines(summary_lines, file.path(log_dir, "analysis_summary.txt"))

if (has_openxlsx) {
  wb <- openxlsx::createWorkbook()
  add_sheet <- function(name, df) {
    if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, as.data.frame(df))
  }
  add_sheet("result_highlights", result_highlights)
  add_sheet("uploaded_gene_list", user_gene_tbl)
  add_sheet("gene_context", gene_context_tbl)
  add_sheet("cohort_qc", cohort_qc)
  add_sheet("sample_annotation", anno_all)
  add_sheet("pathway_catalog", candidate_program_tbl %>% select(gs_name, feature, label, collection, n_genes, overlap_n_seed, overlap_n_ferro, overlap_n_immune, pathway_score))
  add_sheet("program_effects", program_effects)
  add_sheet("program_meta", program_meta)
  if (!is.null(program_loo) && nrow(program_loo) > 0) add_sheet("program_loo", program_loo)
  add_sheet("gene_effects", gene_effects)
  add_sheet("gene_meta", gene_meta)
  if (!is.null(gene_loo) && nrow(gene_loo) > 0) add_sheet("gene_loo", gene_loo)
  if (!is.null(deconv_effects) && nrow(deconv_effects) > 0) add_sheet("deconv_effects", deconv_effects)
  if (!is.null(deconv_meta) && nrow(deconv_meta) > 0) add_sheet("deconv_meta", deconv_meta)
  add_sheet("manuscript_summary", manuscript_summary)
  openxlsx::saveWorkbook(wb, file.path(out_dir, "immunotherapy_response_prediction_summary.xlsx"), overwrite = TRUE)
}

msg("[DONE] Outputs written to {out_dir}")
msg("[DONE] Figures : {fig_dir}")
msg("[DONE] Tables  : {tab_dir}")
