#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
  options(width = 140)
})

# ============================================================
# Rana Salihoglu
# 
#   1) GSVA API compatibility layer for legacy and new releases.
#   2) Explicit warning/error logging for reproducibility.
#   3) More robust statistics (BH-FDR, effect sizes, Stouffer meta-Z).
#   4) Higher-grade figure defaults and deterministic execution.
#   5) Stronger input validation and safer fallbacks.
# ============================================================

# -----------------------------
# CLI
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit[length(hit)] == length(args)) return(default)
  args[hit[length(hit)] + 1]
}

project_root    <- path.expand(arg_value("--project_root", "~/Desktop/PROJECT_16/test"))
install_missing <- tolower(arg_value("--install_missing", "false")) %in% c("true", "1", "yes", "y")
min_gene_set    <- suppressWarnings(as.integer(arg_value("--min_gene_set", "10")))
max_programs    <- suppressWarnings(as.integer(arg_value("--max_programs", "24")))
organism        <- arg_value("--organism", "Homo sapiens")
seed_method     <- suppressWarnings(as.integer(arg_value("--seed", "1234")))

script_arg <- commandArgs(trailingOnly = FALSE)
script_path_raw <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_path <- if (length(script_path_raw) > 0) normalizePath(script_path_raw[[1]], winslash = "/", mustWork = FALSE) else NA_character_
script_name <- if (is.na(script_path)) "03_tme_analysis_papergrade_q1_iobr_driven_updated_full_v9.R" else basename(script_path)
script_stem <- tools::file_path_sans_ext(script_name)
outdir_name <- arg_value("--outdir_name", script_stem)

set.seed(seed_method)

cat("[INFO] project_root:", project_root, "\n")
cat("[INFO] script_name :", script_name, "\n")
cat("[INFO] outdir_name :", outdir_name, "\n")
cat("[INFO] seed        :", seed_method, "\n")

# -----------------------------
# Package management
# -----------------------------
cran_pkgs <- c(
  "data.table", "dplyr", "tibble", "tidyr", "stringr", "purrr", "readr",
  "ggplot2", "forcats", "scales", "pheatmap", "jsonlite", "glue", "msigdbr",
  "viridisLite"
)
bioc_pkgs <- c("GSVA", "clusterProfiler", "IOBR")
optional_cran_pkgs <- c("openxlsx", "ggrepel")
optional_bioc_pkgs <- c("GSEABase")

ensure_pkg <- function(pkg, bioc = FALSE, install_missing = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  if (!install_missing) return(FALSE)
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
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
  stop(
    "Required packages missing: ", paste(missing_core, collapse = ", "),
    ". Re-run with --install_missing true or install them manually."
  )
}
invisible(vapply(optional_cran_pkgs, ensure_pkg, logical(1), bioc = FALSE, install_missing = install_missing))
invisible(vapply(optional_bioc_pkgs, ensure_pkg, logical(1), bioc = TRUE, install_missing = install_missing))

has_iobr     <- requireNamespace("IOBR", quietly = TRUE)
has_gseabase <- requireNamespace("GSEABase", quietly = TRUE)

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
library(msigdbr)
library(GSVA)
library(clusterProfiler)
library(survival)
if (has_iobr) library(IOBR)

has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
has_ggrepel  <- requireNamespace("ggrepel", quietly = TRUE)

# -----------------------------
# Paths
# -----------------------------
root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
out_dir <- file.path("04_TME", outdir_name)
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  out_v44         = file.path(root, "OUT_v44"),
  geo_root        = file.path(root, "GEO_PREP"),
  wgcna_dir       = file.path(root, "GEO_PREP", "GSE81089", "WGCNA_GSE81089_out_v12_signature_fix"),
  ferr_driver     = file.path(root, "ferrdb_driver.csv"),
  ferr_suppressor = file.path(root, "ferrdb_suppressor.csv"),
  ferr_marker     = file.path(root, "ferrdb_marker.csv")
)

# -----------------------------
# Logging / conditions
# -----------------------------
.warn_log <- character()
.err_log  <- character()

append_log <- function(kind = c("WARN", "ERROR", "INFO"), text) {
  kind <- match.arg(kind)
  line <- paste0("[", kind, "] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", text)
  if (kind == "WARN") .warn_log <<- c(.warn_log, line)
  if (kind == "ERROR") .err_log <<- c(.err_log, line)
  cat(line, "\n")
}

msg <- function(...) append_log("INFO", glue(..., .envir = parent.frame()))
warn_msg <- function(...) append_log("WARN", glue(..., .envir = parent.frame()))
err_msg <- function(...) append_log("ERROR", glue(..., .envir = parent.frame()))

withCallingHandlers(
  expr = {
    # no-op; global warning collector activated below as helper wrapper is used locally
  },
  warning = function(w) {
    .warn_log <<- c(.warn_log, paste0("[WARN] ", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }
)

# -----------------------------
# Helpers
# -----------------------------
read_flex <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) fread(path, sep = "\t", data.table = FALSE, ...)
  else fread(path, data.table = FALSE, ...)
}

first_existing <- function(vec) {
  hit <- vec[file.exists(vec)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

safe_write <- function(df, path) {
  if (is.null(df)) return(invisible(NULL))
  readr::write_tsv(tibble::as_tibble(df), path)
}

top_n_rows <- function(df, n) {
  df <- tibble::as_tibble(df)
  if (nrow(df) == 0) return(df)
  dplyr::slice(df, seq_len(min(as.integer(n), nrow(df))))
}

tail_n_rows <- function(df, n) {
  df <- tibble::as_tibble(df)
  if (nrow(df) == 0) return(df)
  n_take <- min(as.integer(n), nrow(df))
  start_idx <- max(1L, nrow(df) - n_take + 1L)
  dplyr::slice(df, start_idx:nrow(df))
}

safe_min_num <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  x_num <- x_num[is.finite(x_num)]
  if (length(x_num) == 0) return(NA_real_)
  min(x_num)
}

safe_neglog10 <- function(x, floor = 1e-300) {
  x_num <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x_num))
  ok <- is.finite(x_num) & !is.na(x_num)
  out[ok] <- -log10(pmax(x_num[ok], floor))
  out
}

safe_factor_reorder <- function(x, by, desc = FALSE, fun = stats::median) {
  x_chr <- as.character(x)
  by_num <- suppressWarnings(as.numeric(by))
  keep <- !is.na(x_chr) & x_chr != "" & is.finite(by_num)
  if (!any(keep)) return(factor(x_chr, levels = unique(x_chr)))
  ord_tbl <- tibble::tibble(x = x_chr[keep], by = by_num[keep]) %>%
    dplyr::group_by(.data$x) %>%
    dplyr::summarise(by = fun(.data$by, na.rm = TRUE), .groups = "drop")
  if (desc) {
    ord_tbl <- ord_tbl %>% dplyr::arrange(dplyr::desc(.data$by), .data$x)
  } else {
    ord_tbl <- ord_tbl %>% dplyr::arrange(.data$by, .data$x)
  }
  factor(x_chr, levels = ord_tbl$x)
}

clean_label <- function(x, width = 42) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  stringr::str_trunc(x, width = width)
}

compact_program_label <- function(gs_name, description = NA_character_, collection = NA_character_, width = 42) {
  gs_name <- as.character(gs_name)
  description <- as.character(description)
  collection <- as.character(collection)

  n <- max(length(gs_name), length(description), length(collection))
  gs_name <- rep_len(gs_name, n)
  description <- rep_len(description, n)
  collection <- rep_len(collection, n)

  out <- gs_name
  is_h  <- !is.na(collection) & collection == "H"
  is_c8 <- !is.na(collection) & collection == "C8"
  use_desc <- !is.na(description) & description != "" & nchar(description) <= 70 & !is_h & !is_c8

  out[is_h] <- stringr::str_remove(gs_name[is_h], "^HALLMARK_")
  out[is_c8] <- stringr::str_remove(gs_name[is_c8], "^[A-Z0-9]+_")
  out[use_desc] <- description[use_desc]

  clean_label(out, width = width)
}

robust_z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(madv) || madv == 0) {
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(x)))
    return((x - mean(x, na.rm = TRUE)) / s)
  }
  (x - med) / madv
}

rank01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  n_ok <- sum(is.finite(x))
  if (n_ok <= 1) return(rep(0.5, length(x)))
  (rank(x, na.last = "keep", ties.method = "average") - 1) / (n_ok - 1)
}

meta_stouffer <- function(p, effect_sign = NULL) {
  p <- suppressWarnings(as.numeric(p))
  keep <- is.finite(p) & !is.na(p) & p > 0 & p <= 1
  if (!any(keep)) return(c(meta_z = NA_real_, meta_p = NA_real_))
  p <- p[keep]
  z <- stats::qnorm(p / 2, lower.tail = FALSE)
  if (!is.null(effect_sign)) {
    sgn <- suppressWarnings(as.numeric(effect_sign))[keep]
    sgn[!is.finite(sgn) | sgn == 0] <- 1
    z <- z * sign(sgn)
  }
  meta_z <- sum(z) / sqrt(length(z))
  meta_p <- 2 * stats::pnorm(abs(meta_z), lower.tail = FALSE)
  c(meta_z = meta_z, meta_p = meta_p)
}

cap_neglog10 <- function(x, cap = 20, floor = 1e-300) {
  y <- safe_neglog10(x, floor = floor)
  y[is.finite(y)] <- pmin(y[is.finite(y)], cap)
  y
}

plot_save <- function(p, file, width = 10, height = 7, dpi = 320) {
  ggplot2::ggsave(filename = file, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  pdf_file <- sub("\\.png$", ".pdf", file, ignore.case = TRUE)
  ggplot2::ggsave(filename = pdf_file, plot = p, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
}

save_pheatmap <- function(mat, file, width = 9, height = 7, ...) {
  dots <- list(...)
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"

  if (!is.null(dots$annotation_col)) {
    ann_col <- as.data.frame(dots$annotation_col, stringsAsFactors = FALSE)
    rownames(ann_col) <- rownames(as.data.frame(dots$annotation_col))
    if ("risk_group" %in% names(ann_col)) {
      ann_col$risk_group <- factor(as.character(ann_col$risk_group), levels = c("Low", "High"))
    }
    if ("cohort" %in% names(ann_col)) {
      ann_col$cohort <- factor(as.character(ann_col$cohort), levels = unique(as.character(ann_col$cohort)))
    }
    dots$annotation_col <- ann_col
  }

  pfun <- function() do.call(pheatmap::pheatmap, c(list(mat = mat), dots))

  png(file, width = width * 300, height = height * 300, res = 300)
  pfun()
  dev.off()

  pdf_file <- sub("\\.png$", ".pdf", file, ignore.case = TRUE)
  grDevices::cairo_pdf(pdf_file, width = width, height = height, onefile = TRUE)
  pfun()
  dev.off()
}

write_output_manifest <- function(base_dir) {
  files <- list.files(base_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  tibble::tibble(
    relative_path = sub(paste0("^", normalizePath(base_dir, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(files, winslash = "/", mustWork = FALSE)),
    bytes = file.info(files)$size
  ) %>%
    dplyr::arrange(.data$relative_path)
}

clean_expr_colnames <- function(x) make.names(x, unique = TRUE)

find_gene_col <- function(df) {
  cand <- names(df)[tolower(names(df)) %in% c("gene", "symbol", "gene_symbol", "genesymbol", "official_symbol", "hub_gene")]
  if (length(cand) > 0) cand[[1]] else names(df)[[1]]
}

extract_gene_tbl <- function(df, source_label) {
  if (is.null(df) || nrow(df) == 0) return(tibble(gene = character(), source = character()))
  gcol <- find_gene_col(df)
  tibble(gene = unique(trimws(as.character(df[[gcol]]))), source = source_label) %>%
    filter(!is.na(.data$gene), .data$gene != "")
}

read_ferrdb <- function(path, label) {
  df <- read_flex(path)
  if (is.null(df) || nrow(df) == 0) return(tibble(gene = character(), source = character()))
  gcol <- find_gene_col(df)
  tibble(gene = unique(trimws(as.character(df[[gcol]]))), source = label) %>%
    filter(!is.na(.data$gene), .data$gene != "")
}

cohort_row_z <- function(mat) {
  if (nrow(mat) == 0 || ncol(mat) == 0) return(mat)
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

simple_ssgsea <- function(expr_mat, gene_sets, min_gene_set = 10, alpha = 0.25) {
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= min_gene_set]
  if (length(gene_sets) == 0) return(matrix(numeric(), nrow = 0, ncol = ncol(expr_mat)))
  ranks <- apply(expr_mat, 2, rank, ties.method = "average")
  rownames(ranks) <- rownames(expr_mat)
  out <- matrix(NA_real_, nrow = length(gene_sets), ncol = ncol(expr_mat),
                dimnames = list(names(gene_sets), colnames(expr_mat)))
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
      Phit <- cumsum(ifelse(hits, ranked_r / sum(ranked_r[hits]), 0))
      Pmiss <- cumsum(ifelse(!hits, 1 / Nm, 0))
      es <- Phit - Pmiss
      out[i, j] <- sum(es)
    }
  }
  out
}

# GSVA compatibility layer for both legacy and current APIs
run_gsva_ssgsea <- function(expr_gene_sample, gene_sets, min_gene_set = 10,
                            kcdf = "Gaussian", abs_ranking = TRUE, ssgsea_norm = TRUE) {
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= min_gene_set]
  if (length(gene_sets) == 0) return(matrix(numeric(), nrow = 0, ncol = ncol(expr_gene_sample)))

  # Prefer new API if available (GSVA >= 2 style)
  if (exists("ssgseaParam", where = asNamespace("GSVA"), mode = "function")) {
    res <- tryCatch({
      if (has_gseabase) {
        gsc <- GSEABase::GeneSetCollection(
          lapply(names(gene_sets), function(nm) GSEABase::GeneSet(geneIds = gene_sets[[nm]], setName = nm))
        )
        param <- GSVA::ssgseaParam(exprData = expr_gene_sample, geneSets = gsc,
                                   normalize = ssgsea_norm, alpha = 0.25)
      } else {
        param <- GSVA::ssgseaParam(exprData = expr_gene_sample, geneSets = gene_sets,
                                   normalize = ssgsea_norm, alpha = 0.25)
      }
      gsva_fun <- get("gsva", envir = asNamespace("GSVA"))
      # some versions accept verbose/BPPARAM, some do not
      tryCatch(
        suppressWarnings(gsva_fun(param, verbose = FALSE)),
        error = function(e) suppressWarnings(gsva_fun(param))
      )
    }, error = function(e) e)
    if (!inherits(res, "error")) return(as.matrix(res))
    warn_msg("New GSVA API path failed; trying legacy interface. Reason: {res$message}")
  }

  # Legacy API
  res2 <- tryCatch({
    suppressWarnings(
      GSVA::gsva(expr = expr_gene_sample, gset.idx.list = gene_sets,
                 method = "ssgsea", ssgsea.norm = ssgsea_norm,
                 kcdf = kcdf, abs.ranking = abs_ranking)
    )
  }, error = function(e) e)
  if (!inherits(res2, "error")) return(as.matrix(res2))

  # Final fallback
  warn_msg("Legacy GSVA API also failed; using internal ssGSEA fallback. Reason: {res2$message}")
  simple_ssgsea(expr_gene_sample, gene_sets, min_gene_set = min_gene_set)
}

compute_wilcox_effect <- function(df, value_col, group_col = "risk_group") {
  x <- as_tibble(df)
  vals <- suppressWarnings(as.numeric(x[[value_col]]))
  grp  <- as.character(x[[group_col]])
  ok <- is.finite(vals) & !is.na(grp)
  vals <- vals[ok]
  grp  <- grp[ok]
  if (length(unique(grp)) < 2) {
    return(c(p = NA_real_, delta = NA_real_, rbc = NA_real_, med_high = NA_real_, med_low = NA_real_, n_high = NA_real_, n_low = NA_real_))
  }
  hi <- vals[grp == "High"]
  lo <- vals[grp == "Low"]
  if (length(hi) < 3 || length(lo) < 3) {
    return(c(p = NA_real_, delta = median(hi, na.rm = TRUE) - median(lo, na.rm = TRUE), rbc = NA_real_, med_high = median(hi, na.rm = TRUE), med_low = median(lo, na.rm = TRUE), n_high = length(hi), n_low = length(lo)))
  }
  wt <- tryCatch(suppressWarnings(wilcox.test(hi, lo, exact = FALSE)), error = function(e) NULL)
  delta <- median(hi, na.rm = TRUE) - median(lo, na.rm = TRUE)
  rbc <- tryCatch({
    pooled <- rank(c(hi, lo), ties.method = "average")
    n1 <- length(hi); n2 <- length(lo)
    w <- sum(pooled[seq_along(hi)]) - n1 * (n1 + 1) / 2
    u <- w
    2 * u / (n1 * n2) - 1
  }, error = function(e) NA_real_)
  c(
    p = if (!is.null(wt)) wt$p.value else NA_real_,
    delta = delta,
    rbc = rbc,
    med_high = median(hi, na.rm = TRUE),
    med_low = median(lo, na.rm = TRUE),
    n_high = length(hi),
    n_low = length(lo)
  )
}

cohort_compare_long <- function(long_df, value_col = "score", feature_col = "feature") {
  bind_rows(lapply(unique(long_df$cohort), function(coh) {
    sub <- long_df %>% filter(.data$cohort == coh)
    bind_rows(lapply(unique(sub[[feature_col]]), function(ft) {
      cur <- sub[sub[[feature_col]] == ft, , drop = FALSE]
      eff <- compute_wilcox_effect(cur, value_col, "risk_group")
      tibble(
        cohort = coh,
        feature = ft,
        p_value = as.numeric(eff[["p"]]),
        delta_high_minus_low = as.numeric(eff[["delta"]]),
        rank_biserial = as.numeric(eff[["rbc"]]),
        median_high = as.numeric(eff[["med_high"]]),
        median_low = as.numeric(eff[["med_low"]]),
        n_high = as.numeric(eff[["n_high"]]),
        n_low = as.numeric(eff[["n_low"]])
      )
    }))
  }))
}

stouffer_meta <- function(p, direction, weights = NULL) {
  p <- suppressWarnings(as.numeric(p))
  direction <- suppressWarnings(as.numeric(direction))
  ok <- is.finite(p) & !is.na(p) & p > 0 & p <= 1 & is.finite(direction) & !is.na(direction)
  if (!any(ok)) return(c(z = NA_real_, p = NA_real_))
  p <- p[ok]
  direction <- direction[ok]
  if (is.null(weights)) weights <- rep(1, length(p)) else weights <- suppressWarnings(as.numeric(weights))[ok]
  weights[!is.finite(weights) | weights <= 0] <- 1
  z_each <- stats::qnorm(pmax(1e-300, 1 - p / 2)) * sign(direction)
  z_meta <- sum(weights * z_each) / sqrt(sum(weights^2))
  p_meta <- 2 * stats::pnorm(-abs(z_meta))
  c(z = z_meta, p = p_meta)
}



compute_cox_stats <- function(df, time_col = "time", event_col = "event", predictor_col = "risk_group") {
  x <- tibble::as_tibble(df)
  tm <- suppressWarnings(as.numeric(x[[time_col]]))
  ev <- suppressWarnings(as.numeric(x[[event_col]]))
  pr <- x[[predictor_col]]

  if (is.character(pr) || is.factor(pr)) {
    pr <- factor(as.character(pr), levels = c("Low", "High"))
  } else {
    pr <- suppressWarnings(as.numeric(pr))
  }

  ok <- is.finite(tm) & is.finite(ev) & !is.na(pr)
  tm <- tm[ok]; ev <- ev[ok]; pr <- pr[ok]
  if (length(tm) < 15 || length(unique(ev)) < 2) {
    return(tibble::tibble(hr = NA_real_, lower95 = NA_real_, upper95 = NA_real_, p_value = NA_real_, n = length(tm), n_events = sum(ev > 0, na.rm = TRUE)))
  }

  fit <- tryCatch({
    survival::coxph(survival::Surv(tm, ev) ~ pr)
  }, error = function(e) NULL)

  if (is.null(fit)) {
    return(tibble::tibble(hr = NA_real_, lower95 = NA_real_, upper95 = NA_real_, p_value = NA_real_, n = length(tm), n_events = sum(ev > 0, na.rm = TRUE)))
  }

  sm <- summary(fit)
  cf <- sm$coefficients
  ci <- sm$conf.int
  tibble::tibble(
    hr = unname(ci[1, "exp(coef)"]),
    lower95 = unname(ci[1, "lower .95"]),
    upper95 = unname(ci[1, "upper .95"]),
    p_value = unname(cf[1, "Pr(>|z|)"]),
    n = length(tm),
    n_events = sum(ev > 0, na.rm = TRUE)
  )
}

meta_loghr <- function(hr, lower95, upper95) {
  hr <- suppressWarnings(as.numeric(hr))
  lower95 <- suppressWarnings(as.numeric(lower95))
  upper95 <- suppressWarnings(as.numeric(upper95))
  ok <- is.finite(hr) & is.finite(lower95) & is.finite(upper95) & hr > 0 & lower95 > 0 & upper95 > 0
  if (!any(ok)) return(c(meta_hr = NA_real_, meta_lower95 = NA_real_, meta_upper95 = NA_real_, meta_p = NA_real_, k = 0))
  loghr <- log(hr[ok])
  se <- (log(upper95[ok]) - log(lower95[ok])) / (2 * 1.96)
  keep <- is.finite(se) & se > 0
  if (!any(keep)) return(c(meta_hr = NA_real_, meta_lower95 = NA_real_, meta_upper95 = NA_real_, meta_p = NA_real_, k = 0))
  loghr <- loghr[keep]
  se <- se[keep]
  w <- 1 / (se^2)
  mu <- sum(w * loghr) / sum(w)
  se_mu <- sqrt(1 / sum(w))
  z <- mu / se_mu
  p <- 2 * stats::pnorm(-abs(z))
  c(meta_hr = exp(mu), meta_lower95 = exp(mu - 1.96 * se_mu), meta_upper95 = exp(mu + 1.96 * se_mu), meta_p = p, k = sum(keep))
}

standardize_for_deconv <- function(expr_gene_sample) {
  expr <- expr_gene_sample
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  expr[is.na(expr)] <- 0
  mode(expr) <- "numeric"
  if (min(expr, na.rm = TRUE) < 0) expr <- expr - min(expr, na.rm = TRUE)
  expr
}

coerce_deconv_result <- function(obj) {
  if (is.null(obj)) return(NULL)
  if (inherits(obj, c("matrix", "data.frame"))) {
    return(as.data.frame(obj, check.names = FALSE))
  }
  if (is.list(obj)) {
    if (length(obj) == 1 && inherits(obj[[1]], c("matrix", "data.frame"))) {
      return(as.data.frame(obj[[1]], check.names = FALSE))
    }
    cand <- names(obj)[vapply(obj, function(x) inherits(x, c("matrix", "data.frame")), logical(1))]
    if (length(cand) >= 1) {
      return(as.data.frame(obj[[cand[[1]]]], check.names = FALSE))
    }
  }
  NULL
}

run_single_deconv_method <- function(expr_gene_sample, method_label) {
  stopifnot(is.matrix(expr_gene_sample) || is.data.frame(expr_gene_sample))

  if (identical(method_label, "quantiseq") && !requireNamespace("limSolve", quietly = TRUE)) {
    warn_msg("IOBR deconvolution skipped for method 'quantiseq' because dependency 'limSolve' is not installed.")
    return(NULL)
  }

  method_specific_args <- list(
    estimate   = list(list(eset = expr_gene_sample, method = method_label),
                      list(eset = expr_gene_sample, method = method_label, arrays = TRUE)),
    mcpcounter = list(list(eset = expr_gene_sample, method = method_label),
                      list(eset = expr_gene_sample, method = method_label, arrays = TRUE)),
    xcell      = list(list(eset = expr_gene_sample, method = method_label, arrays = TRUE),
                      list(eset = expr_gene_sample, method = method_label)),
    quantiseq  = list(list(eset = expr_gene_sample, method = method_label, arrays = FALSE),
                      list(eset = expr_gene_sample, method = method_label),
                      list(eset = expr_gene_sample, method = method_label, arrays = TRUE))
  )

  attempts <- method_specific_args[[method_label]]
  errs <- character(0)
  for (k in seq_along(attempts)) {
    arglist <- attempts[[k]]
    cur <- tryCatch(
      do.call(IOBR::deconvo_tme, arglist),
      error = function(e) {
        errs <<- c(errs, paste0("attempt", k, ": ", conditionMessage(e)))
        NULL
      }
    )
    cur_df <- coerce_deconv_result(cur)
    if (!is.null(cur_df) && nrow(cur_df) > 0 && ncol(cur_df) > 0) {
      msg("IOBR deconvolution succeeded for method '{method_label}' using attempt {k}.")
      return(cur_df)
    }
  }

  warn_msg("IOBR deconvolution failed for method '{method_label}'. {paste(errs, collapse = '; ')}")
  NULL
}

run_deconvolution_suite <- function(expr_gene_sample) {
  if (!has_iobr) return(list())
  methods <- c("estimate", "mcpcounter", "xcell", "quantiseq")
  out <- list()
  for (m in methods) {
    cur <- run_single_deconv_method(expr_gene_sample = expr_gene_sample, method_label = m)
    if (!is.null(cur)) out[[m]] <- cur
  }
  out
}

format_deconv_table <- function(df, method_label, cohort_label, sample_ids) {
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) return(NULL)
  d <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  sample_ids <- unique(as.character(sample_ids))
  sample_ids <- sample_ids[!is.na(sample_ids) & sample_ids != ""]

  rn <- rownames(d)
  cn <- colnames(d)
  first_col <- as.character(d[[1]])

  row_hit <- sum(rn %in% sample_ids)
  col_hit <- sum(cn %in% sample_ids)
  first_hit <- sum(first_col %in% sample_ids)

  # samples in rownames: rows=samples, cols=features
  if (row_hit >= max(2, col_hit, first_hit)) {
    d$sample_id <- rn
    long <- d %>%
      tibble::as_tibble() %>%
      dplyr::relocate(.data$sample_id) %>%
      tidyr::pivot_longer(cols = -sample_id, names_to = "feature", values_to = "score") %>%
      dplyr::mutate(method = method_label, cohort = cohort_label)
    return(long)
  }

  # samples in first column: rows=samples, cols=features
  if (first_hit >= max(2, row_hit, col_hit)) {
    names(d)[1] <- "sample_id"
    long <- d %>%
      tibble::as_tibble() %>%
      tidyr::pivot_longer(cols = -sample_id, names_to = "feature", values_to = "score") %>%
      dplyr::mutate(method = method_label, cohort = cohort_label)
    return(long)
  }

  # samples in column names: rows=features, cols=samples
  feature_vec <- if (!is.null(rn) && all(rn != "")) rn else first_col
  if (col_hit >= 2) {
    if (identical(feature_vec, first_col)) d <- d[, -1, drop = FALSE]
    d$feature <- feature_vec
    long <- d %>%
      tibble::as_tibble() %>%
      dplyr::relocate(.data$feature) %>%
      tidyr::pivot_longer(cols = -feature, names_to = "sample_id", values_to = "score") %>%
      dplyr::mutate(method = method_label, cohort = cohort_label)
    return(long)
  }

  warn_msg("Could not resolve deconvolution table orientation for cohort {cohort_label}, method {method_label}.")
  NULL
}

# -----------------------------
# Read ML outputs
# -----------------------------
ml_external <- read_flex(file.path(paths$out_v44, "external_ensemble_results.csv"))
ml_internal <- read_flex(file.path(paths$out_v44, "internal_loco_leaderboard.csv"))

possible_tbl <- as_tibble(read_flex(first_existing(c(
  file.path(paths$out_v44, "genes", "possible_prognostic_genes_FULLTRAIN.csv"),
  file.path(root, "luad_r_analysis_bundle", "possible_prognostic_genes_FULLTRAIN.csv")
))))
meta_tbl <- as_tibble(read_flex(file.path(paths$out_v44, "genes", "meta_univariate_top_genes__FULLTRAIN.csv")))
locked_tbl <- as_tibble(read_flex(file.path(paths$out_v44, "genes", "locked_signature_rna_genes.csv")))
selected_tbl <- as_tibble(read_flex(file.path(paths$out_v44, "genes", "selected_rna_genes_fulltrain.csv")))
wgcna_tbl <- as_tibble(read_flex(first_existing(c(
  file.path(paths$wgcna_dir, "tables", "wgcna_genes_for_ML.csv"),
  file.path(paths$wgcna_dir, "tables", "conservative_ferroptosis_genes_GS_MM.tsv"),
  file.path(paths$wgcna_dir, "tables", "module_hub_genes_MM_GS.tsv"),
  file.path(paths$wgcna_dir, "tables", "module_genes_GS_MM.tsv")
))))

ferr_all <- bind_rows(
  read_ferrdb(paths$ferr_driver, "FerrDb_driver"),
  read_ferrdb(paths$ferr_suppressor, "FerrDb_suppressor"),
  read_ferrdb(paths$ferr_marker, "FerrDb_marker")
) %>% distinct(gene, .keep_all = TRUE)

sig_sources <- bind_rows(
  extract_gene_tbl(possible_tbl, "possible"),
  extract_gene_tbl(meta_tbl, "meta"),
  extract_gene_tbl(locked_tbl, "locked"),
  extract_gene_tbl(selected_tbl, "selected"),
  extract_gene_tbl(wgcna_tbl, "wgcna"),
  ferr_all
) %>% distinct(gene, source)

master_signature <- sig_sources %>%
  group_by(gene) %>%
  summarise(sources = paste(sort(unique(source)), collapse = "; "), n_sources = n_distinct(source), .groups = "drop") %>%
  arrange(desc(n_sources), gene)
safe_write(master_signature, file.path(tab_dir, "table_ml_signature_master.tsv"))

locked_genes   <- extract_gene_tbl(locked_tbl, "locked")$gene
selected_genes <- extract_gene_tbl(selected_tbl, "selected")$gene
possible_genes <- extract_gene_tbl(possible_tbl, "possible")$gene
seed_genes <- unique(c(locked_genes, selected_genes, possible_genes))
if (length(seed_genes) < 10) seed_genes <- unique(master_signature$gene)
seed_tbl <- tibble(gene = seed_genes)
safe_write(seed_tbl, file.path(tab_dir, "table_ml_seed_genes.tsv"))

# -----------------------------
# Load LUAD GEO cohorts
# -----------------------------
geo_candidates <- c("GSE31210", "GSE136961", "GSE50081", "GSE72094", "GSE68465", "GSE30219")
cohort_data <- list()
for (coh in geo_candidates) {
  expr_path <- file.path(paths$geo_root, coh, "X_rna_symbol.csv")
  surv_path <- file.path(paths$geo_root, coh, "survival.csv")
  if (!file.exists(expr_path) || !file.exists(surv_path)) {
    warn_msg("Cohort {coh} skipped: missing expression or survival file.")
    next
  }
  expr <- read_flex(expr_path)
  surv <- read_flex(surv_path)
  if (is.null(expr) || is.null(surv)) {
    warn_msg("Cohort {coh} skipped: failed to read input tables.")
    next
  }

  sample_col <- names(expr)[1]
  expr_genes <- setdiff(names(expr), sample_col)
  expr_mat <- as.matrix(expr[, expr_genes, drop = FALSE])
  rownames(expr_mat) <- expr[[sample_col]]
  mode(expr_mat) <- "numeric"
  colnames(expr_mat) <- clean_expr_colnames(colnames(expr_mat))

  surv_sample <- names(surv)[1]
  time_idx <- which(tolower(names(surv)) %in% c("os.time", "time", "os_time", "days", "survival_time"))
  event_idx <- which(tolower(names(surv)) %in% c("os", "event", "status", "dead", "death"))
  if (length(time_idx) == 0 || length(event_idx) == 0) {
    warn_msg("Cohort {coh} skipped: survival time/event columns not found.")
    next
  }
  time_col <- names(surv)[time_idx[1]]
  event_col <- names(surv)[event_idx[1]]

  surv2 <- surv %>%
    transmute(
      sample_id = .data[[surv_sample]],
      time = suppressWarnings(as.numeric(.data[[time_col]])),
      event = suppressWarnings(as.numeric(.data[[event_col]]))
    ) %>%
    filter(!is.na(sample_id))

  common_ids <- intersect(rownames(expr_mat), surv2$sample_id)
  if (length(common_ids) < 10) {
    warn_msg("Cohort {coh} skipped: <10 matched samples after expression-survival intersection.")
    next
  }
  expr_mat <- expr_mat[common_ids, , drop = FALSE]
  surv2 <- surv2 %>% filter(sample_id %in% common_ids)
  cohort_data[[coh]] <- list(cohort = coh, expr = expr_mat, surv = surv2)
}

cohort_qc <- bind_rows(lapply(cohort_data, function(obj) {
  tibble(
    cohort = obj$cohort,
    n_samples = nrow(obj$expr),
    n_genes = ncol(obj$expr),
    n_events = sum(obj$surv$event > 0, na.rm = TRUE),
    median_time = median(obj$surv$time, na.rm = TRUE)
  )
}))
safe_write(cohort_qc, file.path(tab_dir, "table_cohort_qc.tsv"))

if (length(cohort_data) < 2) stop("At least two usable cohorts are required.")

# -----------------------------
# Signature scoring per cohort
# -----------------------------
common_genes <- Reduce(intersect, lapply(cohort_data, function(x) colnames(x$expr)))
common_genes <- common_genes[!is.na(common_genes) & common_genes != ""]

cohort_gene_overlap_qc <- bind_rows(lapply(cohort_data, function(obj) {
  sig_in <- intersect(seed_genes, colnames(obj$expr))
  fallback_in <- intersect(master_signature$gene, colnames(obj$expr))
  tibble(
    cohort = obj$cohort,
    n_genes_total = ncol(obj$expr),
    n_seed_genes_present = length(sig_in),
    n_master_signature_genes_present = length(fallback_in),
    n_common_genes_all_cohorts = length(common_genes)
  )
}))
safe_write(cohort_gene_overlap_qc, file.path(tab_dir, "table_cohort_gene_overlap_qc.tsv"))

anno_all <- bind_rows(lapply(cohort_data, function(obj) {
  mat <- obj$expr
  sig_in <- intersect(seed_genes, colnames(mat))
  if (length(sig_in) < 5) sig_in <- intersect(master_signature$gene, colnames(mat))

  if (length(sig_in) >= 5) {
    sig_mat <- mat[, sig_in, drop = FALSE]
    sig_mat_z <- apply(sig_mat, 2, robust_z)
    if (is.vector(sig_mat_z)) sig_mat_z <- matrix(sig_mat_z, ncol = 1, dimnames = list(rownames(sig_mat), sig_in[1]))
    rownames(sig_mat_z) <- rownames(sig_mat)
    colnames(sig_mat_z) <- sig_in
    sig_score_raw <- rowMeans(sig_mat, na.rm = TRUE)
    sig_score_z <- rowMeans(sig_mat_z, na.rm = TRUE)
    sig_score_pct <- rank01(sig_score_z)
  } else {
    sig_score_raw <- rep(NA_real_, nrow(mat))
    sig_score_z <- rep(NA_real_, nrow(mat))
    sig_score_pct <- rep(NA_real_, nrow(mat))
  }

  med_score <- suppressWarnings(median(sig_score_z, na.rm = TRUE))
  if (!is.finite(med_score)) med_score <- 0
  grp <- ifelse(sig_score_z >= med_score, "High", "Low")

  tibble(
    sample_id = rownames(mat),
    cohort = obj$cohort,
    n_signature_genes_used = length(sig_in),
    signature_score_raw = sig_score_raw,
    signature_score = sig_score_z,
    signature_score_percentile = sig_score_pct,
    risk_group = grp,
    time = obj$surv$time[match(rownames(mat), obj$surv$sample_id)],
    event = obj$surv$event[match(rownames(mat), obj$surv$sample_id)]
  )
}))

signature_score_tbl <- anno_all
score_summary <- signature_score_tbl %>%
  group_by(cohort) %>%
  summarise(
    n_samples = n(),
    median_signature_score_z = median(signature_score, na.rm = TRUE),
    mean_signature_score_z = mean(signature_score, na.rm = TRUE),
    median_signature_score_raw = median(signature_score_raw, na.rm = TRUE),
    mean_signature_score_raw = mean(signature_score_raw, na.rm = TRUE),
    median_signature_genes_used = median(n_signature_genes_used, na.rm = TRUE),
    high_fraction = mean(risk_group == "High", na.rm = TRUE),
    .groups = "drop"
  )
safe_write(signature_score_tbl, file.path(tab_dir, "table_signature_scores_by_sample.tsv"))
safe_write(score_summary, file.path(tab_dir, "table_signature_score_summary_by_cohort.tsv"))

signature_survival_tbl <- bind_rows(lapply(split(signature_score_tbl, signature_score_tbl$cohort), function(df) {
  cox_grp <- compute_cox_stats(df, predictor_col = "risk_group")
  cox_z <- compute_cox_stats(df, predictor_col = "signature_score")
  tibble::tibble(
    cohort = unique(df$cohort),
    hr_high_vs_low = cox_grp$hr,
    lower95_high_vs_low = cox_grp$lower95,
    upper95_high_vs_low = cox_grp$upper95,
    p_value_high_vs_low = cox_grp$p_value,
    hr_per_z = cox_z$hr,
    lower95_per_z = cox_z$lower95,
    upper95_per_z = cox_z$upper95,
    p_value_per_z = cox_z$p_value,
    n = cox_grp$n,
    n_events = cox_grp$n_events
  )
}))

sig_meta_grp <- meta_loghr(signature_survival_tbl$hr_high_vs_low, signature_survival_tbl$lower95_high_vs_low, signature_survival_tbl$upper95_high_vs_low)
sig_meta_z <- meta_loghr(signature_survival_tbl$hr_per_z, signature_survival_tbl$lower95_per_z, signature_survival_tbl$upper95_per_z)
signature_survival_meta <- tibble::tibble(
  endpoint = c("High_vs_Low", "Per_1_robust_z"),
  meta_hr = c(sig_meta_grp[["meta_hr"]], sig_meta_z[["meta_hr"]]),
  meta_lower95 = c(sig_meta_grp[["meta_lower95"]], sig_meta_z[["meta_lower95"]]),
  meta_upper95 = c(sig_meta_grp[["meta_upper95"]], sig_meta_z[["meta_upper95"]]),
  meta_p = c(sig_meta_grp[["meta_p"]], sig_meta_z[["meta_p"]]),
  k = c(sig_meta_grp[["k"]], sig_meta_z[["k"]])
)
signature_survival_meta$meta_fdr <- p.adjust(signature_survival_meta$meta_p, method = "BH")

safe_write(signature_survival_tbl, file.path(tab_dir, "table_signature_survival_by_cohort.tsv"))
safe_write(signature_survival_meta, file.path(tab_dir, "table_signature_survival_meta.tsv"))

# -----------------------------
# Package-derived TME program discovery
# -----------------------------
msig_tbl <- bind_rows(
  msigdbr::msigdbr(species = organism, collection = "H")  %>% mutate(program_source = "MSigDB_H"),
  msigdbr::msigdbr(species = organism, collection = "C7") %>% mutate(program_source = "MSigDB_C7"),
  msigdbr::msigdbr(species = organism, collection = "C8") %>% mutate(program_source = "MSigDB_C8")
) %>%
  transmute(gs_name, gene_symbol, gs_description, gs_collection, gs_subcollection = dplyr::coalesce(gs_subcollection, NA_character_), program_source)

term2gene <- msig_tbl %>% distinct(gs_name, gene_symbol)
term2name <- msig_tbl %>% distinct(gs_name, gs_description, gs_collection, gs_subcollection, program_source)

seed_in_msig <- intersect(seed_genes, unique(term2gene$gene_symbol))
if (length(seed_in_msig) < min_gene_set) {
  seed_in_msig <- intersect(master_signature$gene, unique(term2gene$gene_symbol))
}
if (length(seed_in_msig) < min_gene_set) stop("Too few ML genes mapped to MSigDB collections.")

seed_enrichment <- tryCatch({
  clusterProfiler::enricher(
    gene = seed_in_msig,
    TERM2GENE = term2gene,
    TERM2NAME = term2name %>% transmute(gs_name, term_label = gs_description),
    minGSSize = min_gene_set,
    pAdjustMethod = "BH"
  )
}, error = function(e) {
  err_msg("clusterProfiler::enricher failed: {e$message}")
  NULL
})

seed_enrichment_tbl <- if (!is.null(seed_enrichment) && nrow(as.data.frame(seed_enrichment)) > 0) {
  as_tibble(as.data.frame(seed_enrichment)) %>%
    rename(gs_name = ID, description = Description, gene_ratio = GeneRatio, bg_ratio = BgRatio, p_value = pvalue, p_adjust = p.adjust, q_value = qvalue, core_genes = geneID, overlap_size = Count) %>%
    left_join(term2name, by = "gs_name") %>%
    mutate(
      program_source = dplyr::coalesce(.data$program_source, "MSigDB"),
      gs_collection = dplyr::coalesce(.data$gs_collection, NA_character_),
      gs_subcollection = dplyr::coalesce(.data$gs_subcollection, NA_character_),
      gs_description = dplyr::coalesce(.data$gs_description, .data$description, .data$gs_name),
      tme_relevance_score = -log10(pmax(p_adjust, 1e-300)) * overlap_size,
      label = ifelse(is.na(gs_description) | gs_description == "", gs_name, gs_description)
    ) %>%
    arrange(p_adjust, desc(overlap_size), gs_name)
} else {
  tibble()
}

safe_write(seed_enrichment_tbl, file.path(tab_dir, "table_ml_seed_tme_enrichment.tsv"))
if (nrow(seed_enrichment_tbl) == 0) stop("No enriched package-derived TME programs identified from ML genes.")

n_program_sources <- if ("program_source" %in% names(seed_enrichment_tbl)) dplyr::n_distinct(seed_enrichment_tbl$program_source) else 1L
n_program_sources <- max(1L, as.integer(n_program_sources))
per_source_quota <- max(1L, ceiling(max_programs / n_program_sources))

selected_program_tbl <- seed_enrichment_tbl %>%
  filter(p_adjust < 0.10)

if (nrow(selected_program_tbl) > 0) {
  if (!("program_source" %in% names(selected_program_tbl))) selected_program_tbl$program_source <- "MSigDB"
  selected_program_tbl <- selected_program_tbl %>%
    arrange(program_source, p_adjust, desc(overlap_size)) %>%
    group_by(program_source) %>%
    group_modify(~ top_n_rows(.x, per_source_quota)) %>%
    ungroup() %>%
    arrange(p_adjust, desc(overlap_size)) %>%
    top_n_rows(max_programs)
}

if (nrow(selected_program_tbl) < 6) {
  selected_program_tbl <- seed_enrichment_tbl %>%
    arrange(p_adjust, desc(overlap_size)) %>%
    top_n_rows(max_programs)
}

program_catalog <- selected_program_tbl %>%
  mutate(
    display_label = compact_program_label(gs_name, gs_description, gs_collection, width = 56),
    short_label = compact_program_label(gs_name, gs_description, gs_collection, width = 34)
  ) %>%
  select(gs_name, display_label, short_label, label, gs_collection, gs_subcollection, program_source, overlap_size, p_adjust, core_genes)

safe_write(program_catalog, file.path(tab_dir, "table_data_driven_tme_program_catalog.tsv"))
safe_write(program_catalog %>% select(gs_name, display_label, short_label), file.path(tab_dir, "table_program_label_map.tsv"))

program_gene_sets <- split(
  term2gene$gene_symbol[term2gene$gs_name %in% program_catalog$gs_name],
  term2gene$gs_name[term2gene$gs_name %in% program_catalog$gs_name]
)
program_gene_sets <- lapply(program_gene_sets, unique)
program_gene_sets <- program_gene_sets[vapply(program_gene_sets, length, integer(1)) >= min_gene_set]

label_map <- program_catalog %>% select(gs_name, display_label, short_label)
program_labels <- setNames(label_map$display_label, label_map$gs_name)
program_short_labels <- setNames(label_map$short_label, label_map$gs_name)

program_long <- bind_rows(lapply(cohort_data, function(obj) {
  cohort_name <- obj$cohort
  expr_use <- t(obj$expr)
  expr_use <- cohort_row_z(expr_use)
  available_sets <- lapply(program_gene_sets, function(gs) intersect(gs, rownames(expr_use)))
  available_sets <- available_sets[vapply(available_sets, length, integer(1)) >= min_gene_set]
  if (length(available_sets) < 3) {
    warn_msg("Cohort {cohort_name}: <3 retained TME programs after cohort-specific gene intersection.")
    return(NULL)
  }

  msg("cohort {cohort_name}: retained package-derived TME programs = {length(available_sets)}")
  program_scores <- run_gsva_ssgsea(expr_gene_sample = expr_use, gene_sets = available_sets, min_gene_set = min_gene_set)
  if (is.null(program_scores) || nrow(program_scores) == 0) {
    warn_msg("Cohort {cohort_name}: no TME program scores generated.")
    return(NULL)
  }

  score_long <- as_tibble(as.data.frame(t(program_scores)), rownames = "sample_id") %>%
    pivot_longer(cols = -sample_id, names_to = "gs_name", values_to = "score") %>%
    mutate(
      feature = ifelse(.data$gs_name %in% names(program_labels), unname(program_labels[.data$gs_name]), clean_label(.data$gs_name, width = 56)),
      feature_short = ifelse(.data$gs_name %in% names(program_short_labels), unname(program_short_labels[.data$gs_name]), clean_label(.data$gs_name, width = 34)),
      cohort = cohort_name
    ) %>%
    select(sample_id, cohort, feature, feature_short, score, gs_name) %>%
    left_join(
      anno_all %>% select(sample_id, cohort, risk_group, signature_score, time, event),
      by = c("sample_id", "cohort")
    ) %>%
    relocate(sample_id, cohort, risk_group, signature_score, time, event, feature, score)
  score_long
}))

if (nrow(program_long) == 0) stop("No TME pathway scores were generated in the available cohorts.")

program_score_df <- program_long %>%
  select(sample_id, cohort, feature, score) %>%
  pivot_wider(names_from = feature, values_from = score) %>%
  left_join(anno_all, by = c("sample_id", "cohort")) %>%
  relocate(sample_id, cohort, signature_score, risk_group, time, event)

safe_write(program_score_df, file.path(tab_dir, "table_tme_scores_by_sample.tsv"))
safe_write(program_long, file.path(tab_dir, "table_tme_scores_long.tsv"))

program_effects <- cohort_compare_long(program_long, value_col = "score", feature_col = "feature") %>%
  group_by(feature) %>%
  mutate(fdr_feature = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(fdr_global = p.adjust(p_value, method = "BH"))
safe_write(program_effects, file.path(tab_dir, "table_tme_effects_by_cohort.tsv"))

program_meta <- bind_rows(lapply(split(program_effects, program_effects$feature), function(df) {
  mz <- stouffer_meta(df$p_value, df$delta_high_minus_low, weights = pmax(1, df$n_high + df$n_low))
  tibble(
    feature = unique(df$feature),
    n_cohorts = sum(is.finite(df$delta_high_minus_low)),
    mean_delta = ifelse(any(is.finite(df$delta_high_minus_low)), mean(df$delta_high_minus_low, na.rm = TRUE), NA_real_),
    median_delta = ifelse(any(is.finite(df$delta_high_minus_low)), median(df$delta_high_minus_low, na.rm = TRUE), NA_real_),
    mean_rank_biserial = ifelse(any(is.finite(df$rank_biserial)), mean(df$rank_biserial, na.rm = TRUE), NA_real_),
    min_fdr_feature = safe_min_num(df$fdr_feature),
    min_fdr_global = safe_min_num(df$fdr_global),
    stouffer_z = as.numeric(mz[["z"]]),
    stouffer_p = as.numeric(mz[["p"]])
  )
})) %>%
  mutate(stouffer_fdr = p.adjust(stouffer_p, method = "BH")) %>%
  arrange(stouffer_fdr, min_fdr_global, desc(abs(mean_delta)))
safe_write(program_meta, file.path(tab_dir, "table_tme_meta_summary.tsv"))

# -----------------------------
# IOBR-based deconvolution
# -----------------------------
deconv_long_all <- tibble()
deconv_effects <- tibble()
if (has_iobr) {
  for (obj in cohort_data) {
    cohort_name <- obj$cohort
    expr_raw_gene_sample <- t(obj$expr)
    expr_raw_gene_sample <- standardize_for_deconv(expr_raw_gene_sample)
    deconv_res <- run_deconvolution_suite(expr_raw_gene_sample)
    cur_long <- bind_rows(lapply(names(deconv_res), function(m) format_deconv_table(deconv_res[[m]], method_label = m, cohort_label = cohort_name, sample_ids = rownames(obj$expr))))
    if (nrow(cur_long) > 0) {
      cur_long <- cur_long %>% left_join(anno_all %>% select(sample_id, cohort, risk_group, signature_score), by = c("sample_id", "cohort"))
      deconv_long_all <- bind_rows(deconv_long_all, cur_long)
    }
  }
}

if (nrow(deconv_long_all) > 0) {
  deconv_long_all <- deconv_long_all %>%
    mutate(
      score = suppressWarnings(as.numeric(score)),
      feature = clean_label(feature, width = 42),
      sample_id = as.character(sample_id)
    ) %>%
    filter(is.finite(score)) %>%
    left_join(anno_all %>% select(sample_id, cohort, risk_group, signature_score), by = c("sample_id", "cohort"), suffix = c("", ".y")) %>%
    mutate(
      risk_group = dplyr::coalesce(.data$risk_group, .data$risk_group.y),
      signature_score = dplyr::coalesce(.data$signature_score, .data$signature_score.y)
    ) %>%
    select(-matches("\\.y$")) %>%
    filter(!is.na(risk_group))

  safe_write(deconv_long_all, file.path(tab_dir, "table_tme_deconvolution_scores_long.tsv"))

  deconv_effects <- deconv_long_all %>%
    mutate(feature_key = paste(method, feature, sep = "::")) %>%
    cohort_compare_long(value_col = "score", feature_col = "feature_key") %>%
    rename(feature_key = feature) %>%
    tidyr::separate(feature_key, into = c("method", "feature"), sep = "::", remove = FALSE) %>%
    group_by(method, feature) %>%
    mutate(fdr_feature = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(fdr_global = p.adjust(p_value, method = "BH"))
  safe_write(deconv_effects, file.path(tab_dir, "table_tme_deconvolution_effects_by_cohort.tsv"))

  deconv_meta <- bind_rows(lapply(split(deconv_effects, paste(deconv_effects$method, deconv_effects$feature, sep = "::")), function(df) {
    mz <- stouffer_meta(df$p_value, df$rank_biserial, weights = pmax(1, df$n_high + df$n_low))
    tibble(
      method = unique(df$method),
      feature = unique(df$feature),
      n_cohorts = sum(is.finite(df$delta_high_minus_low)),
      mean_delta = ifelse(any(is.finite(df$delta_high_minus_low)), mean(df$delta_high_minus_low, na.rm = TRUE), NA_real_),
      median_delta = ifelse(any(is.finite(df$delta_high_minus_low)), median(df$delta_high_minus_low, na.rm = TRUE), NA_real_),
      mean_rank_biserial = ifelse(any(is.finite(df$rank_biserial)), mean(df$rank_biserial, na.rm = TRUE), NA_real_),
      median_rank_biserial = ifelse(any(is.finite(df$rank_biserial)), median(df$rank_biserial, na.rm = TRUE), NA_real_),
      direction_consistency = ifelse(any(is.finite(df$rank_biserial)), mean(sign(df$rank_biserial[is.finite(df$rank_biserial)]) == sign(stats::median(df$rank_biserial, na.rm = TRUE))), NA_real_),
      min_fdr_feature = safe_min_num(df$fdr_feature),
      min_fdr_global = safe_min_num(df$fdr_global),
      stouffer_z = as.numeric(mz[["z"]]),
      stouffer_p = as.numeric(mz[["p"]])
    )
  })) %>%
    mutate(
      stouffer_fdr = p.adjust(stouffer_p, method = "BH"),
      method_label = dplyr::recode(method,
        estimate = "ESTIMATE",
        mcpcounter = "MCPcounter",
        quantiseq = "quanTIseq",
        xcell = "xCell",
        .default = method
      ),
      feature_short = clean_label(feature, width = 30),
      feature_label = clean_label(paste(method_label, feature_short, sep = " | "), width = 42)
    ) %>%
    arrange(stouffer_fdr, min_fdr_global, desc(abs(mean_rank_biserial)))
  safe_write(deconv_meta, file.path(tab_dir, "table_tme_deconvolution_meta_summary.tsv"))
} else {
  warn_msg("No deconvolution results available; IOBR absent or all methods failed.")
  deconv_meta <- tibble()
}

# -----------------------------
# Figures
# -----------------------------
# FIG01: Signature score distribution (within-cohort standardized)
p <- signature_score_tbl %>%
  mutate(
    cohort = factor(cohort, levels = cohort_qc$cohort),
    risk_group = factor(risk_group, levels = c("Low", "High"))
  ) %>%
  ggplot(aes(x = risk_group, y = signature_score, fill = risk_group)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.80, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95) +
  facet_wrap(~ cohort, ncol = 3, scales = "fixed") +
  labs(
    title = "Within-cohort standardized prognostic signature scores",
    x = NULL,
    y = "Signature score (robust z)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )
plot_save(p, file.path(fig_dir, "FIG01_signature_score_distribution.png"), 11, 7)

# FIG01B: Signature survival forest
sig_forest <- signature_survival_tbl %>%
  filter(is.finite(hr_high_vs_low), is.finite(lower95_high_vs_low), is.finite(upper95_high_vs_low)) %>%
  mutate(cohort = factor(cohort, levels = rev(cohort_qc$cohort)))
if (nrow(sig_forest) > 0) {
  p <- ggplot(sig_forest, aes(y = cohort, x = hr_high_vs_low, xmin = lower95_high_vs_low, xmax = upper95_high_vs_low)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
    geom_errorbar(aes(y = cohort, xmin = lower95_high_vs_low, xmax = upper95_high_vs_low),
                  orientation = "y", width = 0.18, color = "grey35") +
    geom_point(size = 2.6, color = "#B2182B") +
    scale_x_log10() +
    labs(title = "Prognostic effect of signature-defined high-risk group", x = "Hazard ratio (log scale)", y = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  plot_save(p, file.path(fig_dir, "FIG01B_signature_survival_forest.png"), 8.5, 4.5)
}

# FIG02: Heatmap of top pathway meta-signals
top_programs_heat <- program_meta %>%
  filter(is.finite(mean_delta), !is.na(stouffer_fdr)) %>%
  arrange(stouffer_fdr, desc(abs(mean_delta))) %>%
  top_n_rows(18) %>%
  pull(feature)

agg_program <- program_long %>%
  filter(feature %in% top_programs_heat) %>%
  group_by(feature, cohort, risk_group) %>%
  summarise(score = mean(as.numeric(score), na.rm = TRUE), .groups = "drop") %>%
  mutate(col_id = paste0(cohort, ".", risk_group))

if (nrow(agg_program) > 0) {
  mat_program <- agg_program %>%
    select(feature, col_id, score) %>%
    pivot_wider(names_from = col_id, values_from = score) %>%
    as.data.frame()
  rownames(mat_program) <- mat_program$feature
  mat_program$feature <- NULL
  mat_program <- as.matrix(mat_program)
  mat_program <- cohort_row_z(mat_program)
  ann_col <- data.frame(
    cohort = sub("\\..*$", "", colnames(mat_program)),
    risk_group = factor(sub("^.*\\.", "", colnames(mat_program)), levels = c("Low", "High")),
    stringsAsFactors = FALSE,
    row.names = colnames(mat_program)
  )
  save_pheatmap(
    mat_program,
    file.path(fig_dir, "FIG02_tme_program_heatmap.png"),
    width = 9, height = 8.5,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = ann_col,
    fontsize_row = 8,
    fontsize_col = 9,
    border_color = NA,
    main = "Top package-derived TME programs across cohorts"
  )
}

# FIG03: Pathway meta-summary
meta_plot_tbl <- program_meta %>%
  filter(is.finite(mean_delta), !is.na(stouffer_fdr)) %>%
  arrange(stouffer_fdr, desc(abs(mean_delta))) %>%
  top_n_rows(15) %>%
  mutate(
    feature_short = clean_label(feature, width = 38),
    feature_short = safe_factor_reorder(feature_short, mean_delta),
    neglog10_fdr = safe_neglog10(stouffer_fdr)
  )

if (nrow(meta_plot_tbl) > 0) {
  p <- ggplot(meta_plot_tbl, aes(x = mean_delta, y = feature_short, color = neglog10_fdr, size = n_cohorts)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_point() +
    scale_color_viridis_c(option = "C", na.value = "grey80") +
    labs(
      title = "Meta-summary of package-derived TME programs",
      x = "Mean delta (High - Low)",
      y = NULL,
      color = "-log10 FDR",
      size = "Cohorts"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  plot_save(p, file.path(fig_dir, "FIG03_tme_meta_summary.png"), 10, 6.5)
}

# FIG04: Cohort-wise forest-style panels
top_features <- program_meta %>%
  filter(is.finite(mean_delta), !is.na(stouffer_fdr)) %>%
  arrange(stouffer_fdr, desc(abs(mean_delta))) %>%
  top_n_rows(8) %>%
  pull(feature)

forest_tbl <- program_effects %>%
  filter(feature %in% top_features, is.finite(delta_high_minus_low)) %>%
  mutate(
    feature_short = clean_label(feature, width = 34),
    sig_flag = if_else(!is.na(fdr_global) & fdr_global < 0.05, "FDR<0.05", "NS")
  )

if (nrow(forest_tbl) > 0) {
  p <- ggplot(forest_tbl, aes(x = delta_high_minus_low, y = cohort, color = sig_flag)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
    geom_point(size = 2.8) +
    facet_wrap(~ feature_short, scales = "free_x") +
    scale_color_manual(values = c("FDR<0.05" = "#B2182B", "NS" = "grey45")) +
    labs(
      title = "Cohort-wise differences in top TME programs",
      x = "Delta (High - Low)",
      y = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  plot_save(p, file.path(fig_dir, "FIG04_tme_cohort_forest_panels.png"), 13, 9)
}

# FIG05: Enrichment overview with compact labels
bar_tbl <- selected_program_tbl %>%
  mutate(label_short = compact_program_label(gs_name, gs_description, gs_collection, width = 45)) %>%
  filter(is.finite(tme_relevance_score), !is.na(label_short), label_short != "") %>%
  arrange(tme_relevance_score) %>%
  tail_n_rows(15) %>%
  mutate(label_short = safe_factor_reorder(label_short, tme_relevance_score))

if (nrow(bar_tbl) > 0) {
  p <- ggplot(bar_tbl, aes(x = tme_relevance_score, y = label_short, fill = program_source)) +
    geom_col() +
    labs(
      title = "ML-derived genes enrich package-derived TME programs",
      x = "Relevance score (-log10 FDR × overlap)",
      y = NULL,
      fill = "Source"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
  plot_save(p, file.path(fig_dir, "FIG05_ml_seed_tme_enrichment_barplot.png"), 10.5, 7.5)
}

# FIG06: Deconvolution heatmap restricted to informative features
if (exists("deconv_meta") && nrow(deconv_meta) > 0 && nrow(deconv_long_all) > 0) {
  top_deconv <- deconv_meta %>%
    filter(is.finite(mean_delta), !is.na(stouffer_fdr)) %>%
    arrange(stouffer_fdr, desc(abs(mean_delta))) %>%
    top_n_rows(20)

  heat_df <- deconv_long_all %>%
    semi_join(top_deconv, by = c("method", "feature")) %>%
    group_by(method, feature, cohort, risk_group) %>%
    summarise(score = mean(as.numeric(score), na.rm = TRUE), .groups = "drop") %>%
    mutate(
      row_feature = clean_label(paste(method, feature, sep = " | "), width = 46),
      col_id = paste0(cohort, ".", risk_group)
    )

  if (nrow(heat_df) > 0) {
    mat <- heat_df %>%
      select(row_feature, col_id, score) %>%
      pivot_wider(names_from = col_id, values_from = score) %>%
      as.data.frame()
    rownames(mat) <- mat$row_feature
    mat$row_feature <- NULL
    mat <- as.matrix(mat)
    mat <- cohort_row_z(mat)
    ann_col <- data.frame(
      cohort = sub("\\..*$", "", colnames(mat)),
      risk_group = factor(sub("^.*\\.", "", colnames(mat)), levels = c("Low", "High")),
      stringsAsFactors = FALSE,
      row.names = colnames(mat)
    )
    save_pheatmap(
      mat,
      file.path(fig_dir, "FIG06_tme_deconvolution_heatmap.png"),
      width = 9.5, height = 10,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_col = ann_col,
      fontsize_row = 7,
      fontsize_col = 9,
      border_color = NA,
      main = "Top deconvolution-derived TME features"
    )
  }
}

# FIG07: Deconvolution meta-summary (method-stratified, scale-free effect)
if (exists("deconv_meta") && nrow(deconv_meta) > 0) {
  deconv_meta_plot <- deconv_meta %>%
    filter(is.finite(mean_rank_biserial), !is.na(stouffer_fdr), n_cohorts >= 3) %>%
    arrange(stouffer_fdr, desc(abs(mean_rank_biserial))) %>%
    group_by(method_label) %>%
    slice_head(n = 6) %>%
    ungroup() %>%
    mutate(
      method_label = factor(method_label, levels = c("ESTIMATE", "MCPcounter", "quanTIseq", "xCell")),
      feature_label = safe_factor_reorder(feature_label, mean_rank_biserial),
      neglog10_fdr = cap_neglog10(stouffer_fdr, cap = 20)
    )

  if (nrow(deconv_meta_plot) > 0) {
    p <- ggplot(deconv_meta_plot, aes(x = mean_rank_biserial, y = feature_label, color = neglog10_fdr, size = n_cohorts)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
      geom_point() +
      facet_wrap(~ method_label, scales = "free_y", ncol = 2) +
      scale_color_viridis_c(option = "C", na.value = "grey80", limits = c(0, 20)) +
      scale_x_continuous(limits = c(-1, 1)) +
      labs(
        title = "Meta-summary of deconvolution-derived TME features",
        subtitle = "Scale-free effect size shown as mean rank-biserial correlation",
        x = "Mean rank-biserial effect (High vs Low)",
        y = NULL,
        color = "-log10 FDR
(capped at 20)",
        size = "Cohorts"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "right",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold")
      )
    plot_save(p, file.path(fig_dir, "FIG07_tme_deconvolution_meta_summary.png"), 12, 8)
  }
}

# FIG08: Representative pathway distributions
sample_top_programs <- program_meta %>%
  filter(is.finite(mean_delta), !is.na(stouffer_fdr)) %>%
  arrange(stouffer_fdr, desc(abs(mean_delta))) %>%
  top_n_rows(9) %>%
  pull(feature)

violin_long <- program_long %>%
  filter(feature %in% sample_top_programs) %>%
  filter(!is.na(risk_group), !is.na(score)) %>%
  mutate(
    feature_short = clean_label(feature, width = 36),
    risk_group = factor(risk_group, levels = c("Low", "High"))
  )

if (nrow(violin_long) > 0) {
  p <- ggplot(violin_long, aes(x = risk_group, y = score, fill = risk_group)) +
    geom_violin(trim = FALSE, alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.95) +
    facet_wrap(~ feature_short, scales = "free_y", ncol = 3) +
    labs(
      title = "Representative package-derived TME programs by prognostic risk",
      x = NULL,
      y = "ssGSEA score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
  plot_save(p, file.path(fig_dir, "FIG08_tme_violin_panels.png"), 12, 10)
} else {
  warn_msg("FIG08 skipped because no pathway scores with non-missing risk_group were available.")
}

# -----------------------------
# Highlights / logs / workbook
# -----------------------------
result_highlights <- tibble(
  metric = c(
    "n_seed_genes",
    "n_common_genes_across_cohorts",
    "median_seed_genes_present_per_cohort",
    "n_package_derived_programs_scored",
    "n_iobr_deconvolution_methods_available",
    "n_deconvolution_features_scored",
    "best_internal_family",
    "best_internal_loco_mean_cindex",
    "signature_meta_hr_high_vs_low",
    "signature_meta_fdr_high_vs_low"
  ),
  value = c(
    length(seed_genes),
    length(common_genes),
    median(cohort_gene_overlap_qc$n_seed_genes_present, na.rm = TRUE),
    dplyr::n_distinct(program_long$feature),
    ifelse(nrow(deconv_long_all) > 0, n_distinct(deconv_long_all$method), 0),
    ifelse(nrow(deconv_long_all) > 0, n_distinct(paste(deconv_long_all$method, deconv_long_all$feature)), 0),
    if (!is.null(ml_internal) && nrow(ml_internal) > 0 && "family" %in% names(ml_internal) && "loco_mean_cindex" %in% names(ml_internal)) as.character(ml_internal$family[[which.max(ml_internal$loco_mean_cindex)]]) else NA,
    if (!is.null(ml_internal) && nrow(ml_internal) > 0 && "loco_mean_cindex" %in% names(ml_internal)) round(max(ml_internal$loco_mean_cindex, na.rm = TRUE), 4) else NA,
    if (exists("signature_survival_meta") && nrow(signature_survival_meta) > 0) round(signature_survival_meta$meta_hr[signature_survival_meta$endpoint == "High_vs_Low"], 4) else NA,
    if (exists("signature_survival_meta") && nrow(signature_survival_meta) > 0) signif(signature_survival_meta$meta_fdr[signature_survival_meta$endpoint == "High_vs_Low"], 3) else NA
  )
)
safe_write(result_highlights, file.path(tab_dir, "table_result_highlights.tsv"))

run_metadata <- tibble(
  field = c(
    "script_name", "script_stem", "project_root", "output_dir", "analysis_time",
    "install_missing", "min_gene_set", "max_programs", "organism", "seed",
    "has_iobr", "has_gseabase", "R_version", "GSVA_version", "IOBR_version"
  ),
  value = c(
    script_name, script_stem, root, out_dir, as.character(Sys.time()),
    as.character(install_missing), as.character(min_gene_set), as.character(max_programs), organism, as.character(seed_method),
    as.character(has_iobr), as.character(has_gseabase), R.version.string,
    as.character(utils::packageVersion("GSVA")),
    if (has_iobr) as.character(utils::packageVersion("IOBR")) else NA_character_
  )
)
safe_write(run_metadata, file.path(tab_dir, "table_run_metadata.tsv"))

warning_tbl <- tibble(log = unique(.warn_log))
error_tbl   <- tibble(log = unique(.err_log))
safe_write(warning_tbl, file.path(tab_dir, "table_warning_log.tsv"))
safe_write(error_tbl, file.path(tab_dir, "table_error_log.tsv"))
writeLines(unique(.warn_log), con = file.path(log_dir, "warnings.log"))
writeLines(unique(.err_log), con = file.path(log_dir, "errors.log"))

writeLines(c(
  "Tumor Microenvironment (TME) package-driven analysis summary",
  paste0("Project root: ", root),
  paste0("Generated at: ", Sys.time()),
  paste0("Output directory: ", out_dir),
  "",
  "Design:",
  "- ML-derived genes are used as the seed feature space.",
  "- No manually hard-coded checkpoint, HLA/APM, or TME marker lists are used.",
  "- TME programs are discovered from package-provided MSigDB collections using msigdbr + clusterProfiler.",
  "- Sample-level pathway activity is quantified with GSVA ssGSEA using a version-compatible wrapper.",
  "- Cell-contexture is quantified with package-native deconvolution through IOBR when available.",
  "- Program- and deconvolution-level meta-significance is summarized with Stouffer meta-analysis and BH control.",
  "- All outputs are exported as manuscript-style tables and figures."
), con = file.path(log_dir, "analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), con = file.path(log_dir, "sessionInfo.txt"))

if (has_openxlsx) {
  wb <- openxlsx::createWorkbook()
  add_sheet <- function(name, df) {
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, as.data.frame(df))
  }
  add_sheet("result_highlights", result_highlights)
  add_sheet("run_metadata", run_metadata)
  add_sheet("cohort_qc", cohort_qc)
  add_sheet("signature_master", master_signature)
  add_sheet("seed_genes", seed_tbl)
  add_sheet("signature_scores", signature_score_tbl)
  add_sheet("score_summary", score_summary)
  add_sheet("signature_survival", signature_survival_tbl)
  add_sheet("signature_survival_meta", signature_survival_meta)
  add_sheet("seed_enrichment", seed_enrichment_tbl)
  add_sheet("program_catalog", program_catalog)
  add_sheet("program_scores", program_score_df)
  add_sheet("program_effects", program_effects)
  add_sheet("program_meta", program_meta)
  if (nrow(deconv_long_all) > 0) add_sheet("deconv_scores", deconv_long_all)
  if (nrow(deconv_effects) > 0) add_sheet("deconv_effects", deconv_effects)
  if (exists("deconv_meta") && nrow(deconv_meta) > 0) add_sheet("deconv_meta", deconv_meta)
  add_sheet("warning_log", warning_tbl)
  add_sheet("error_log", error_tbl)
  openxlsx::saveWorkbook(wb, file.path(out_dir, paste0(script_stem, "_summary.xlsx")), overwrite = TRUE)
}

output_manifest <- write_output_manifest(out_dir)
safe_write(output_manifest, file.path(tab_dir, "table_output_manifest.tsv"))

msg("Outputs written to {out_dir}")
msg("Figures : {fig_dir}")
msg("Tables  : {tab_dir}")

