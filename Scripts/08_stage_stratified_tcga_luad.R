#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  req <- c("data.table","dplyr","tibble","readr","stringr","purrr","ggplot2","survival")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) stop("Missing required packages: ", paste(miss, collapse = ", "))
  invisible(lapply(req, library, character.only = TRUE))
})

`%||%` <- function(a,b) if (!is.null(a) && length(a) > 0 && !all(is.na(a)) && !identical(a, "")) a else b
msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[length(idx)] == length(args)) return(default)
  args[idx[length(idx)] + 1]
}
has_flag <- function(flag) flag %in% args

manifest_file <- get_arg("--manifest", "/mnt/data/full_paths.txt")
project_root  <- get_arg("--root", NA_character_)
genes_file    <- get_arg("--genes", "/mnt/data/possible_prognostic_genes_FULLTRAIN.csv")
outdir        <- get_arg("--outdir", file.path(getwd(), "stage_tcga_luad_q1_results_enhanced"))
tcga_only     <- has_flag("--tcga-only")
seed          <- as.integer(get_arg("--seed", "20260418"))
min_stage_n   <- as.integer(get_arg("--min-stage-n", "20"))
min_stage_e   <- as.integer(get_arg("--min-stage-events", "5"))
bootstrap_B   <- as.integer(get_arg("--bootstrap-B", "200"))
cal_B         <- as.integer(get_arg("--calibration-B", "200"))
time_horizons <- suppressWarnings(as.numeric(strsplit(get_arg("--time-horizons", "1,3,5"), ",")[[1]]))
install_optional <- has_flag("--install-optional")
min_overlap_genes <- as.integer(get_arg("--min-overlap-genes", "5"))

set.seed(seed)
options(stringsAsFactors = FALSE, scipen = 999)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(outdir, "figures")
tab_dir <- file.path(outdir, "tables")
log_dir <- file.path(outdir, "logs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# utilities
# ----------------------------
read_manifest <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x <- trimws(x)
  x[nzchar(x)]
}

infer_root_from_manifest <- function(paths) {
  if (length(paths) == 0) return(NA_character_)
  anchor <- paths[grepl("/(TCGA_GDC_download|GEO_PREP|OUT_v44|bioinf_scripts)/", paths)]
  if (length(anchor) > 0) {
    m <- regexpr("/(TCGA_GDC_download|GEO_PREP|OUT_v44|bioinf_scripts)/", anchor[1], perl = TRUE)
    if (m[1] > 1) return(substr(anchor[1], 1, m[1] - 1))
  }
  NA_character_
}

manifest_paths <- read_manifest(manifest_file)
if (is.na(project_root) || !nzchar(project_root)) project_root <- infer_root_from_manifest(manifest_paths)
if (is.na(project_root) || !nzchar(project_root)) project_root <- getwd()
project_root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
msg("Resolved project root: %s", project_root)

resolve_from_manifest <- function(patterns, manifest_paths, fallback = character()) {
  for (pt in patterns) {
    hit <- manifest_paths[grepl(pt, manifest_paths, ignore.case = TRUE, perl = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  for (x in fallback) {
    if (!is.na(x) && nzchar(x) && file.exists(x)) return(x)
  }
  NA_character_
}

sanitize_gene <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- gsub("\\.", "-", x)
  toupper(x)
}

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- gsub("_", "-", x)
  x <- gsub("\\.", "-", x)
  x
}

clean_tcga_patient <- function(x) {
  x <- clean_id(x)
  sub("^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}).*$", "\\1", x, perl = TRUE)
}

read_flex <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","txt")) return(data.table::fread(path, sep = "\t", data.table = FALSE))
  data.table::fread(path, data.table = FALSE)
}

pick_col <- function(df, patterns, required = TRUE) {
  nms <- names(df)
  low <- tolower(nms)
  for (p in patterns) {
    idx <- grep(p, low, perl = TRUE)
    if (length(idx) > 0) return(nms[idx[1]])
  }
  if (required) stop("Missing required column among patterns: ", paste(patterns, collapse = ", "))
  NA_character_
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

normalize_event <- function(x) {
  if (is.numeric(x)) return(ifelse(is.na(x), NA_real_, ifelse(x > 0, 1, 0)))
  z <- toupper(trimws(as.character(x)))
  out <- ifelse(grepl("DEAD|DECEASED|EVENT|^1$|TRUE|YES", z), 1,
                ifelse(grepl("ALIVE|LIVING|CENSOR|^0$|FALSE|NO", z), 0, suppressWarnings(as.numeric(z))))
  ifelse(is.na(out), NA_real_, ifelse(out > 0, 1, 0))
}

normalize_time_years <- function(x) {
  y <- safe_num(x)
  med <- suppressWarnings(median(y, na.rm = TRUE))
  if (!is.finite(med)) return(y)
  if (med > 100) return(y / 365.25)
  if (med > 15) return(y / 12)
  y
}

normalize_stage <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[z %in% c("", "NA", "N/A", "NULL", "NOT REPORTED", "[NOT AVAILABLE]")] <- NA_character_
  z <- gsub("STAGE\\s*", "", z)
  z <- gsub("[ABC]$", "", z)
  z <- gsub("[^IVX]", "", z)
  out <- dplyr::case_when(
    z == "I" ~ "I",
    z == "II" ~ "II",
    z == "III" ~ "III",
    z == "IV" ~ "IV",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("I","II","III","IV"), ordered = TRUE)
}

normalize_sex <- function(x) {
  z <- toupper(trimws(as.character(x)))
  out <- dplyr::case_when(
    z %in% c("MALE","M") ~ "Male",
    z %in% c("FEMALE","F") ~ "Female",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("Female","Male"))
}

normalize_tnm <- function(x, prefix = c("T","N","M")) {
  prefix <- match.arg(prefix)
  z <- toupper(trimws(as.character(x)))
  z[z %in% c("", "NA", "N/A", "NULL", "NOT REPORTED", "[NOT AVAILABLE]", "MX", "NX", "TX")] <- NA_character_
  z <- gsub(paste0("^.*(", prefix, "[0-4XISABC]+).*$"), "\\1", z, perl = TRUE)
  z <- gsub("[ABC]$", "", z)
  z[!grepl(paste0("^", prefix, "(IS|[0-4])$"), z)] <- NA_character_
  factor(z)
}

stage_group <- function(stage4) {
  dplyr::case_when(
    is.na(stage4) ~ NA_character_,
    stage4 %in% c("I","II") ~ "Early (I-II)",
    stage4 %in% c("III","IV") ~ "Late (III-IV)",
    TRUE ~ NA_character_
  )
}

complete_cases_for_formula <- function(df, formula_obj) {
  vars <- unique(all.vars(formula_obj))
  vars <- vars[vars %in% names(df)]
  vars <- unique(c("os_time", "os_event", vars))
  vars <- vars[vars %in% names(df)]
  if (length(vars) == 0) return(rep(FALSE, nrow(df)))
  complete.cases(as.data.frame(df[, vars, drop = FALSE]))
}

collapse_factor_levels <- function(x, max_levels = 4) {
  if (!is.factor(x)) x <- factor(x)
  tb <- sort(table(x), decreasing = TRUE)
  keep <- names(tb)[seq_len(min(length(tb), max_levels))]
  y <- as.character(x)
  y[!(y %in% keep)] <- "Other"
  factor(y)
}

usable_covariates <- function(df, covars, min_non_missing = 20, min_smallest_group = 5) {
  keep <- vapply(covars, function(v) {
    if (!(v %in% names(df))) return(FALSE)
    x <- df[[v]]
    x <- x[!is.na(x)]
    if (length(x) < min_non_missing) return(FALSE)
    if (is.numeric(x) || is.integer(x)) {
      ux <- unique(as.numeric(x))
      ux <- ux[is.finite(ux)]
      return(length(ux) >= 2 && stats::sd(as.numeric(x), na.rm = TRUE) > 0)
    }
    tab <- table(as.character(x))
    length(tab) >= 2 && min(tab) >= min_smallest_group
  }, logical(1))
  covars[keep]
}

safe_coxph <- function(formula_obj, data) {
  if (is.null(data) || nrow(data) < 20) return(NULL)
  fit <- tryCatch(
    suppressWarnings(survival::coxph(formula_obj, data = data, x = TRUE, singular.ok = TRUE, ties = "efron")),
    error = function(e) NULL
  )
  fit
}

build_formula <- function(include_risk = TRUE, covars = character()) {
  rhs <- covars
  if (include_risk) rhs <- c("risk_score", rhs)
  rhs <- unique(rhs)
  as.formula(paste("survival::Surv(os_time, os_event) ~", paste(rhs, collapse = " + ")))
}

safe_scalar_numeric <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  y <- suppressWarnings(as.numeric(x[1]))
  if (length(y) == 0 || !is.finite(y)) return(NA_real_)
  y
}

safe_write_tsv <- function(df, path) readr::write_tsv(tibble::as_tibble(df), path)

optional_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)
maybe_install_optional <- function(pkgs) {
  if (!install_optional) return(invisible(NULL))
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) install.packages(miss, repos = "https://cloud.r-project.org")
}

save_plot <- function(p, stem, width = 8, height = 6) {
  ggplot2::ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
}

km_plot_basic <- function(df, title) {
  sf <- survival::survfit(survival::Surv(os_time, os_event) ~ risk_group, data = df)
  sdat <- summary(sf)
  plt <- data.frame(time = sdat$time, surv = sdat$surv, strata = sdat$strata)
  ggplot2::ggplot(plt, ggplot2::aes(time, surv, color = strata)) +
    ggplot2::geom_step(linewidth = 0.95) +
    ggplot2::labs(title = title, x = "Time (years)", y = "Overall survival", color = NULL) +
    ggplot2::theme_bw(base_size = 12)
}

cox_extract <- function(fit, term_pattern) {
  sm <- summary(fit)
  rn <- rownames(sm$coefficients)
  idx <- grep(term_pattern, rn, perl = TRUE)
  if (length(idx) == 0) return(NULL)
  i <- idx[1]
  tibble::tibble(
    term = rn[i],
    beta = unname(sm$coefficients[i, "coef"]),
    HR = unname(sm$conf.int[i, "exp(coef)"]),
    lower95 = unname(sm$conf.int[i, "lower .95"]),
    upper95 = unname(sm$conf.int[i, "upper .95"]),
    p_value = unname(sm$coefficients[i, "Pr(>|z|)"])
  )
}

tidy_cox_table <- function(fit, model_name) {
  sm <- summary(fit)
  tibble::tibble(
    term = rownames(sm$coefficients),
    beta = sm$coefficients[, "coef"],
    HR = sm$conf.int[, "exp(coef)"],
    lower95 = sm$conf.int[, "lower .95"],
    upper95 = sm$conf.int[, "upper .95"],
    p_value = sm$coefficients[, "Pr(>|z|)"],
    model = model_name
  )
}

compute_cindex <- function(df, formula = survival::Surv(os_time, os_event) ~ risk_score) {
  fit <- tryCatch(survival::coxph(formula, data = df, x = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  lp <- predict(fit, type = "lp")
  cc <- tryCatch(survival::concordance(survival::Surv(df$os_time, df$os_event) ~ lp)$concordance, error = function(e) NA_real_)
  as.numeric(cc)
}

bootstrap_optimism <- function(df, formula, B = 200) {
  keep <- complete_cases_for_formula(df, formula)
  df <- df[keep, , drop = FALSE]
  if (nrow(df) < 80 || sum(df$os_event) < 20) {
    return(tibble::tibble(apparent_cindex = compute_cindex(df, formula), optimism = NA_real_, corrected_cindex = NA_real_))
  }
  app <- compute_cindex(df, formula)
  opt <- rep(NA_real_, B)
  n <- nrow(df)
  for (b in seq_len(B)) {
    ix <- sample.int(n, n, replace = TRUE)
    boot <- df[ix, , drop = FALSE]
    fit <- tryCatch(survival::coxph(formula, data = boot, x = TRUE), error = function(e) NULL)
    if (is.null(fit)) next
    lp_boot <- predict(fit, newdata = boot, type = "lp")
    lp_orig <- predict(fit, newdata = df, type = "lp")
    c_boot <- tryCatch(survival::concordance(survival::Surv(boot$os_time, boot$os_event) ~ lp_boot)$concordance, error = function(e) NA_real_)
    c_orig <- tryCatch(survival::concordance(survival::Surv(df$os_time, df$os_event) ~ lp_orig)$concordance, error = function(e) NA_real_)
    opt[b] <- c_boot - c_orig
  }
  optimism <- mean(opt, na.rm = TRUE)
  tibble::tibble(apparent_cindex = app, optimism = optimism, corrected_cindex = app - optimism)
}

meta_random_hr <- function(tbl) {
  d <- tbl %>% dplyr::filter(is.finite(beta), is.finite(lower95), is.finite(upper95))
  if (nrow(d) < 2) return(NULL)
  se <- (log(d$upper95) - log(d$lower95)) / (2 * 1.96)
  yi <- d$beta
  wi <- 1 / (se^2)
  mu_fixed <- sum(wi * yi) / sum(wi)
  Q <- sum(wi * (yi - mu_fixed)^2)
  C <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- max((Q - (nrow(d) - 1)) / C, 0)
  wi_re <- 1 / (se^2 + tau2)
  mu_re <- sum(wi_re * yi) / sum(wi_re)
  se_re <- sqrt(1 / sum(wi_re))
  z <- mu_re / se_re
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  i2 <- ifelse(Q <= 0, 0, max((Q - (nrow(d) - 1)) / Q, 0) * 100)
  tibble::tibble(
    k = nrow(d), pooled_beta = mu_re, pooled_HR = exp(mu_re),
    lower95 = exp(mu_re - 1.96 * se_re), upper95 = exp(mu_re + 1.96 * se_re),
    p_value = p, I2 = i2, tau2 = tau2
  )
}

compute_fixed_horizon_risk <- function(fit, newdata, time_year) {
  bh <- tryCatch(survival::basehaz(fit, centered = FALSE), error = function(e) NULL)
  if (is.null(bh) || nrow(bh) == 0) return(rep(NA_real_, nrow(newdata)))
  H0t <- stats::approx(bh$time, bh$hazard, xout = time_year, method = "linear", rule = 2)$y
  lp <- stats::predict(fit, newdata = newdata, type = "lp")
  risk <- 1 - exp(-H0t * exp(lp))
  as.numeric(risk)
}

calc_group_summary <- function(df, var) {

safe_group_pvalue <- function(df, var, group_var = "risk_group") {
  if (!(var %in% names(df)) || !(group_var %in% names(df))) return(NA_real_)
  x <- df[[var]]
  g <- df[[group_var]]
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]
  if (length(x) < 4 || length(unique(g)) < 2) return(NA_real_)
  if (is.numeric(x) || is.integer(x)) {
    p <- tryCatch(stats::wilcox.test(x ~ g)$p.value, error = function(e) NA_real_)
    return(as.numeric(p))
  }
  tb <- table(as.character(x), as.character(g))
  if (nrow(tb) < 2 || ncol(tb) < 2) return(NA_real_)
  expected <- tryCatch(suppressWarnings(stats::chisq.test(tb)$expected), error = function(e) NULL)
  if (is.null(expected)) return(NA_real_)
  if (any(expected < 5) || any(tb < 5)) {
    p <- tryCatch(stats::fisher.test(tb)$p.value, error = function(e) NA_real_)
  } else {
    p <- tryCatch(stats::chisq.test(tb, correct = FALSE)$p.value, error = function(e) NA_real_)
  }
  as.numeric(p)
}

  x <- df[[var]]
  if (is.numeric(x)) {
    high <- x[df$risk_group == "High"]
    low  <- x[df$risk_group == "Low"]
    pval <- tryCatch(stats::wilcox.test(high, low)$p.value, error = function(e) NA_real_)
    tibble::tibble(
      variable = var,
      level = c("High","Low"),
      n = c(sum(df$risk_group == "High" & !is.na(x)), sum(df$risk_group == "Low" & !is.na(x))),
      summary = c(sprintf("%.2f ± %.2f", mean(high, na.rm = TRUE), stats::sd(high, na.rm = TRUE)),
                  sprintf("%.2f ± %.2f", mean(low, na.rm = TRUE), stats::sd(low, na.rm = TRUE))),
      p_value = c(pval, pval)
    )
  } else {
    tb <- table(as.character(x), df$risk_group, useNA = "no")
    pval <- tryCatch(stats::chisq.test(tb)$p.value, error = function(e) NA_real_)
    tibble::as_tibble(as.data.frame.matrix(tb), rownames = "level") %>%
      dplyr::mutate(variable = var, summary = paste0("High=", High %||% 0, "; Low=", Low %||% 0), p_value = pval) %>%
      dplyr::select(variable, level, summary, p_value)
  }
}

standardize_score <- function(x) {
  x <- as.numeric(x)
  z <- (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  z[!is.finite(z)] <- 0
  z
}

time_roc_by_models <- function(df, model_list, use_times) {
  if (!optional_pkg("timeROC")) return(tibble::tibble())
  suppressPackageStartupMessages(library(timeROC))
  out <- list()
  for (nm in names(model_list)) {
    fit <- model_list[[nm]]
    if (is.null(fit)) next
    lp <- stats::predict(fit, newdata = df, type = "lp")
    roc <- tryCatch(timeROC::timeROC(T = df$os_time, delta = df$os_event, marker = lp, cause = 1, times = use_times, iid = TRUE), error = function(e) NULL)
    if (is.null(roc)) next
    out[[nm]] <- tibble::tibble(model = nm, time_year = use_times, AUC = roc$AUC)
  }
  dplyr::bind_rows(out)
}

extract_ph_table <- function(fit, model_name) {
  zph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  if (is.null(zph)) return(tibble::tibble())
  tb <- as.data.frame(zph$table)
  tb$term <- rownames(tb)
  rownames(tb) <- NULL
  tibble::as_tibble(tb) %>%
    dplyr::rename(rho = 1, chisq = 2, p_value = 3) %>%
    dplyr::mutate(model = model_name)
}

subgroup_cox <- function(df, subgroup_var, min_n = 20, min_events = 5) {
  if (!(subgroup_var %in% names(df))) return(tibble::tibble())
  levs <- unique(as.character(df[[subgroup_var]]))
  levs <- levs[!is.na(levs) & nzchar(levs)]
  out <- list()
  for (lv in levs) {
    d <- df[df[[subgroup_var]] == lv & !is.na(df[[subgroup_var]]), , drop = FALSE]
    if (nrow(d) < min_n || sum(d$os_event, na.rm = TRUE) < min_events) next
    fit <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_score, d)
    if (is.null(fit)) next
    tmp <- cox_extract(fit, "^risk_score$")
    if (is.null(tmp)) next
    out[[lv]] <- tmp %>% dplyr::mutate(subgroup = subgroup_var, level = lv, n = nrow(d), events = sum(d$os_event, na.rm = TRUE))
  }
  dplyr::bind_rows(out)
}

make_stage_assoc_table <- function(df) {
  out <- list()
  if ("stage" %in% names(df)) {
    p1 <- tryCatch(stats::kruskal.test(risk_score ~ stage, data = df)$p.value, error = function(e) NA_real_)
    out[[1]] <- tibble::tibble(test = "Kruskal-Wallis", variable = "stage", p_value = p1)
  }
  if ("stage_group" %in% names(df)) {
    p2 <- tryCatch(stats::wilcox.test(risk_score ~ stage_group, data = df)$p.value, error = function(e) NA_real_)
    out[[2]] <- tibble::tibble(test = "Wilcoxon", variable = "stage_group", p_value = p2)
  }
  dplyr::bind_rows(out)
}

# ----------------------------
# readers
# ----------------------------
read_gene_table <- function(path) {
  df <- read_flex(path)
  if (is.null(df) || nrow(df) == 0) stop("Gene file not found or empty: ", path)
  gene_col <- pick_col(df, c("^gene$", "gene.*symbol", "symbol", "feature"), required = FALSE)
  if (is.na(gene_col)) gene_col <- names(df)[1]
  z_col    <- pick_col(df, c("z_meta", "^z$", "meta.*z"), required = FALSE)
  bag_col  <- pick_col(df, c("bag_frac", "bag_freq", "freq"), required = FALSE)
  p_col    <- pick_col(df, c("p_meta", "p_value", "^p$", "fdr"), required = FALSE)
  out <- tibble::tibble(
    gene = sanitize_gene(df[[gene_col]]),
    z_meta = if (!is.na(z_col)) safe_num(df[[z_col]]) else NA_real_,
    bag_frac = if (!is.na(bag_col)) safe_num(df[[bag_col]]) else NA_real_,
    p_meta = if (!is.na(p_col)) safe_num(df[[p_col]]) else NA_real_
  ) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      z_meta = suppressWarnings(max(z_meta, na.rm = TRUE)),
      bag_frac = suppressWarnings(max(bag_frac, na.rm = TRUE)),
      p_meta = suppressWarnings(min(p_meta, na.rm = TRUE)),
      .groups = "drop"
    )
  out$z_meta[!is.finite(out$z_meta)] <- 0
  out$bag_frac[!is.finite(out$bag_frac)] <- 1
  out$p_meta[!is.finite(out$p_meta)] <- NA_real_
  out
}

read_expr_matrix <- function(path, cohort = "cohort") {
  df <- read_flex(path)
  if (is.null(df) || nrow(df) == 0) stop("Expression file not readable: ", path)
  id_col <- pick_col(df, c("sample_id", "sample", "patient", "submitter", "barcode", "geo_accession", "^id$"), required = FALSE)
  if (is.na(id_col)) id_col <- names(df)[1]
  gene_cols <- setdiff(names(df), id_col)
  if (length(gene_cols) < 5) stop("Expression file has too few gene columns: ", path)
  mat <- as.matrix(df[, gene_cols, drop = FALSE])
  mode(mat) <- "numeric"
  rownames(mat) <- clean_id(df[[id_col]])
  colnames(mat) <- sanitize_gene(colnames(mat))
  mat <- mat[!duplicated(rownames(mat)) & nzchar(rownames(mat)), , drop = FALSE]
  keep <- apply(mat, 2, function(v) any(is.finite(v)))
  mat <- mat[, keep, drop = FALSE]
  q99 <- suppressWarnings(as.numeric(stats::quantile(mat, 0.99, na.rm = TRUE)))
  if (is.finite(q99) && q99 > 30) {
    mat <- log2(mat + 1)
    msg("Applied log2(x+1) transform to %s expression", cohort)
  }
  mat
}

read_survival_clinical <- function(surv_path, clinical_path = NA_character_, tcga = FALSE) {
  surv <- read_flex(surv_path)
  clin <- read_flex(clinical_path)
  if (is.null(surv)) stop("Survival file not readable: ", surv_path)
  sid <- pick_col(surv, c("sample_id", "sample", "patient", "submitter", "barcode", "geo_accession", "^id$"), required = FALSE)
  if (is.na(sid)) sid <- names(surv)[1]
  time_col <- pick_col(surv, c("os.*time", "survival.*time", "futime", "time_to_event", "days", "time"), required = FALSE)
  event_col <- pick_col(surv, c("os.*event", "os.*status", "vital.*status", "status", "event", "fustat", "os$"), required = FALSE)
  if (is.na(time_col) && all(c("days_to_death","days_to_last_follow_up") %in% names(surv))) {
    surv$.__os_time__ <- dplyr::coalesce(safe_num(surv$days_to_death), safe_num(surv$days_to_last_follow_up))
    time_col <- ".__os_time__"
  }
  if (is.na(event_col) && "vital_status" %in% names(surv)) event_col <- "vital_status"
  if (is.na(time_col) || is.na(event_col)) stop("Could not identify survival columns in: ", surv_path)
  out <- tibble::tibble(
    sample_id = clean_id(surv[[sid]]),
    os_time = normalize_time_years(surv[[time_col]]),
    os_event = normalize_event(surv[[event_col]])
  )
  if (tcga) out$patient_id <- clean_tcga_patient(out$sample_id) else out$patient_id <- out$sample_id

  if (!is.null(clin)) {
    csid <- pick_col(clin, c("sample_id", "sample", "patient", "submitter", "barcode", "geo_accession", "^id$"), required = FALSE)
    if (is.na(csid)) csid <- names(clin)[1]
    clin2 <- tibble::tibble(sample_id = clean_id(clin[[csid]]))
    if (tcga) clin2$patient_id <- clean_tcga_patient(clin2$sample_id) else clin2$patient_id <- clin2$sample_id

    stg <- pick_col(clin, c("ajcc_pathologic_stage", "pathologic_stage", "tumor_stage", "clinical_stage", "stage"), required = FALSE)
    age_col <- pick_col(clin, c("^age$", "age_at.*diagnosis", "age_at_initial_pathologic_diagnosis", "diagnosis_age", "age_at_diagnosis", "days_to_birth"), required = FALSE)
    sex_col <- pick_col(clin, c("^sex$", "^gender$"), required = FALSE)
    t_col <- pick_col(clin, c("pathologic_t", "ajcc_pathologic_t", "^t_stage$", "^t$"), required = FALSE)
    n_col <- pick_col(clin, c("pathologic_n", "ajcc_pathologic_n", "^n_stage$", "^n$"), required = FALSE)
    m_col <- pick_col(clin, c("pathologic_m", "ajcc_pathologic_m", "^m_stage$", "^m$"), required = FALSE)
    smoke_col <- pick_col(clin, c("smok", "tobacco"), required = FALSE)

    if (!is.na(stg)) clin2$stage <- normalize_stage(clin[[stg]])
    if (!is.na(age_col)) {
      agev <- safe_num(clin[[age_col]])
      if (grepl("days_to_birth", tolower(age_col))) agev <- abs(agev) / 365.25
      clin2$age <- agev
    }
    if (!is.na(sex_col)) clin2$sex <- normalize_sex(clin[[sex_col]])
    if (!is.na(t_col)) clin2$t_stage <- normalize_tnm(clin[[t_col]], "T")
    if (!is.na(n_col)) clin2$n_stage <- normalize_tnm(clin[[n_col]], "N")
    if (!is.na(m_col)) clin2$m_stage <- normalize_tnm(clin[[m_col]], "M")
    if (!is.na(smoke_col)) clin2$smoking <- as.character(clin[[smoke_col]])

    clin2 <- clin2 %>% dplyr::distinct(patient_id, .keep_all = TRUE) %>% dplyr::select(-dplyr::any_of("sample_id"))
    out <- out %>% dplyr::left_join(clin2, by = "patient_id")
  }

  if (!"stage" %in% names(out)) out$stage <- NA
  if (!"age" %in% names(out)) out$age <- NA_real_
  if (!"sex" %in% names(out)) out$sex <- NA_character_
  if (!"t_stage" %in% names(out)) out$t_stage <- NA_character_
  if (!"n_stage" %in% names(out)) out$n_stage <- NA_character_
  if (!"m_stage" %in% names(out)) out$m_stage <- NA_character_
  if (!"smoking" %in% names(out)) out$smoking <- NA_character_

  out %>%
    dplyr::filter(!is.na(.data$sample_id), nzchar(.data$sample_id), is.finite(.data$os_time), !is.na(.data$os_event)) %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::summarise(
      sample_id = dplyr::first(.data$sample_id),
      os_time = max(.data$os_time, na.rm = TRUE),
      os_event = max(.data$os_event, na.rm = TRUE),
      stage = dplyr::first(.data$stage),
      age = dplyr::first(.data$age),
      sex = dplyr::first(.data$sex),
      t_stage = dplyr::first(.data$t_stage),
      n_stage = dplyr::first(.data$n_stage),
      m_stage = dplyr::first(.data$m_stage),
      smoking = dplyr::first(.data$smoking),
      .groups = "drop"
    ) %>%
    dplyr::mutate(stage_group = stage_group(as.character(.data$stage)))
}

align_tcga_expr_to_patients <- function(expr_mat) {
  pid <- clean_tcga_patient(rownames(expr_mat))
  split_idx <- split(seq_along(pid), pid)
  out <- t(vapply(split_idx, function(ix) colMeans(expr_mat[ix, , drop = FALSE], na.rm = TRUE), numeric(ncol(expr_mat))))
  rownames(out) <- names(split_idx)
  colnames(out) <- colnames(expr_mat)
  out
}

score_signature <- function(expr_mat, weights_tbl) {
  common <- intersect(colnames(expr_mat), weights_tbl$gene)
  if (length(common) < min_overlap_genes) stop("Too few overlapping genes to score signature")
  x <- expr_mat[, common, drop = FALSE]
  z <- scale(x)
  z[!is.finite(z)] <- 0
  w <- weights_tbl$weight[match(common, weights_tbl$gene)]
  score <- as.numeric(z %*% w) / sum(abs(w), na.rm = TRUE)
  names(score) <- rownames(expr_mat)
  score
}

resolve_cohort_paths <- function(cohort) {
  expr_pat <- paste0("/", cohort, "/X_rna_(symbol|tpm|tpm_patient)\\.csv$")
  surv_pat <- paste0("/", cohort, "/survival\\.csv$")
  clin_pat <- paste0("/", cohort, "/clinical\\.csv$")
  expr <- resolve_from_manifest(c(expr_pat), manifest_paths)
  surv <- resolve_from_manifest(c(surv_pat), manifest_paths)
  clin <- resolve_from_manifest(c(clin_pat), manifest_paths)
  list(expr = expr, surv = surv, clin = clin)
}

# ----------------------------
# resolve paths
# ----------------------------
tcga_paths <- resolve_cohort_paths("TCGA-LUAD_GDC")
if (is.na(tcga_paths$expr) || is.na(tcga_paths$surv)) stop("TCGA-LUAD full paths could not be resolved from manifest")

ext_cohorts <- c("GSE50081", "GSE68465", "GSE30219")
ext_paths <- lapply(ext_cohorts, resolve_cohort_paths)
names(ext_paths) <- ext_cohorts

if (is.na(genes_file) || !file.exists(genes_file)) {
  genes_file <- resolve_from_manifest(c("/OUT_v44/genes/possible_prognostic_genes_FULLTRAIN\\.csv$", "/possible_prognostic_genes_FULLTRAIN\\.csv$"), manifest_paths, fallback = c(genes_file))
}
if (is.na(genes_file) || !file.exists(genes_file)) stop("Gene file not found")
msg("Gene file: %s", genes_file)

# ----------------------------
# discovery data
# ----------------------------
msg("Loading TCGA-LUAD discovery cohort")
tcga_expr0 <- read_expr_matrix(tcga_paths$expr, "TCGA-LUAD")
tcga_expr <- align_tcga_expr_to_patients(tcga_expr0)
tcga_surv <- read_survival_clinical(tcga_paths$surv, tcga_paths$clin, tcga = TRUE)
common_tcga <- intersect(rownames(tcga_expr), tcga_surv$patient_id)
tcga_expr <- tcga_expr[common_tcga, , drop = FALSE]
tcga_surv <- tcga_surv %>% dplyr::filter(patient_id %in% common_tcga) %>% dplyr::arrange(match(patient_id, rownames(tcga_expr)))
stopifnot(identical(tcga_surv$patient_id, rownames(tcga_expr)))

cand <- read_gene_table(genes_file)
common_genes <- intersect(cand$gene, colnames(tcga_expr))
if (length(common_genes) < min_overlap_genes) stop("Too few candidate genes overlap with TCGA-LUAD expression")
msg("Candidate genes overlapping TCGA-LUAD: %d", length(common_genes))

# ----------------------------
# gene-level evidence
# ----------------------------
base_gene_tbl <- lapply(common_genes, function(g) {
  d <- tibble::tibble(os_time = tcga_surv$os_time, os_event = tcga_surv$os_event, expr = tcga_expr[, g], stage = tcga_surv$stage_group)
  fit_uni <- tryCatch(survival::coxph(survival::Surv(os_time, os_event) ~ expr, data = d), error = function(e) NULL)
  fit_adj <- tryCatch(survival::coxph(survival::Surv(os_time, os_event) ~ expr + stage, data = d), error = function(e) NULL)
  uni <- if (!is.null(fit_uni)) cox_extract(fit_uni, "^expr$") else NULL
  adj <- if (!is.null(fit_adj)) cox_extract(fit_adj, "^expr$") else NULL
  tibble::tibble(
    gene = g,
    beta_uni = uni$beta %||% NA_real_, p_uni = uni$p_value %||% NA_real_,
    beta_adj = adj$beta %||% NA_real_, p_adj = adj$p_value %||% NA_real_,
    hr_adj = adj$HR %||% NA_real_, lower95_adj = adj$lower95 %||% NA_real_, upper95_adj = adj$upper95 %||% NA_real_
  )
}) %>% dplyr::bind_rows() %>%
  dplyr::left_join(cand, by = "gene") %>%
  dplyr::mutate(
    dir_beta = ifelse(is.finite(beta_adj) & beta_adj != 0, sign(beta_adj), ifelse(z_meta == 0, 1, sign(z_meta))),
    evidence_weight = pmax(abs(z_meta), 0.25) * pmax(bag_frac, 0.25),
    weight = dir_beta * evidence_weight,
    abs_weight = abs(weight),
    fdr_adj = p.adjust(p_adj, method = "BH")
  ) %>%
  dplyr::arrange(dplyr::desc(abs_weight), p_adj, p_uni)

safe_write_tsv(base_gene_tbl, file.path(tab_dir, "table_tcga_gene_level_cox.tsv"))
safe_write_tsv(base_gene_tbl %>% dplyr::slice_head(n = 25), file.path(tab_dir, "table_top25_signature_genes.tsv"))

locked_tbl <- base_gene_tbl %>% dplyr::filter(is.finite(beta_adj), !is.na(p_adj), p_adj < 0.2)
if (nrow(locked_tbl) < 5) locked_tbl <- base_gene_tbl %>% dplyr::slice_head(n = min(15, n()))
locked_tbl <- locked_tbl %>% dplyr::mutate(weight = ifelse(weight == 0 | !is.finite(weight), sign(dplyr::coalesce(beta_adj, z_meta, 1)), weight))
safe_write_tsv(locked_tbl, file.path(tab_dir, "final_locked_signature_tcga.tsv"))

# ----------------------------
# scoring and baseline tables
# ----------------------------
risk_score_tcga <- score_signature(tcga_expr, locked_tbl %>% dplyr::select(gene, weight))
tcga_df <- tcga_surv %>%
  dplyr::mutate(
    risk_score = risk_score_tcga[patient_id],
    risk_score_z = standardize_score(risk_score),
    risk_group = factor(ifelse(risk_score >= stats::median(risk_score, na.rm = TRUE), "High", "Low"), levels = c("Low","High"))
  ) %>%
  dplyr::mutate(
    age = ifelse(is.finite(age), age, NA_real_),
    sex = if ("sex" %in% names(.)) sex else NA,
    t_stage = if ("t_stage" %in% names(.)) collapse_factor_levels(t_stage, 4) else NA,
    n_stage = if ("n_stage" %in% names(.)) collapse_factor_levels(n_stage, 4) else NA,
    m_stage = if ("m_stage" %in% names(.)) collapse_factor_levels(m_stage, 3) else NA
  )

safe_write_tsv(tcga_df, file.path(tab_dir, "tcga_patient_risk_scores.tsv"))

cohort_qc_tbl <- tibble::tibble(
  cohort = "TCGA-LUAD",
  n = nrow(tcga_df),
  events = sum(tcga_df$os_event, na.rm = TRUE),
  event_rate = mean(tcga_df$os_event, na.rm = TRUE),
  median_followup_years = stats::median(tcga_df$os_time, na.rm = TRUE),
  stage_I = sum(tcga_df$stage == "I", na.rm = TRUE),
  stage_II = sum(tcga_df$stage == "II", na.rm = TRUE),
  stage_III = sum(tcga_df$stage == "III", na.rm = TRUE),
  stage_IV = sum(tcga_df$stage == "IV", na.rm = TRUE),
  median_age = stats::median(tcga_df$age, na.rm = TRUE)
)
safe_write_tsv(cohort_qc_tbl, file.path(tab_dir, "table_tcga_cohort_qc.tsv"))

clin_vars <- intersect(c("age","sex","stage_group","t_stage","n_stage","m_stage"), names(tcga_df))
clin_summary_tbl <- dplyr::bind_rows(lapply(clin_vars, function(v) calc_group_summary(tcga_df, v)))
safe_write_tsv(clin_summary_tbl, file.path(tab_dir, "table_clinicopathologic_summary_by_risk.tsv"))
safe_write_tsv(make_stage_assoc_table(tcga_df), file.path(tab_dir, "table_risk_stage_association.tsv"))

# ----------------------------
# overall and stage models
# ----------------------------
fit_cont <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_score + stage_group, tcga_df)
fit_bin  <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_group + stage_group, tcga_df)
fit_int  <- tryCatch(survival::coxph(survival::Surv(os_time, os_event) ~ risk_score * stage_group, data = tcga_df), error = function(e) NULL)
fit_cont_z <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_score_z + stage_group, tcga_df)

overall_tbl <- dplyr::bind_rows(
  if (!is.null(fit_cont)) cox_extract(fit_cont, "^risk_score$") %>% dplyr::mutate(model = "Continuous + stage_group") else NULL,
  if (!is.null(fit_cont_z)) cox_extract(fit_cont_z, "^risk_score_z$") %>% dplyr::mutate(model = "Z-standardized continuous + stage_group") else NULL,
  if (!is.null(fit_bin)) cox_extract(fit_bin, "^risk_group") %>% dplyr::mutate(model = "Binary + stage_group") else NULL
) %>% dplyr::select(model, dplyr::everything())
if (!is.null(fit_int) && !is.null(fit_cont)) {
  an <- tryCatch(anova(fit_cont, fit_int, test = "LRT"), error = function(e) NULL)
  interaction_p <- safe_scalar_numeric(an[[ncol(an)]][2])
} else interaction_p <- NA_real_
safe_write_tsv(overall_tbl, file.path(tab_dir, "table_tcga_overall_models.tsv"))
safe_write_tsv(tibble::tibble(interaction_p = interaction_p), file.path(tab_dir, "table_tcga_stage_interaction.tsv"))

stage_levels <- c("I","II","III","IV")
stage_tbl <- list()
for (st in stage_levels) {
  d <- tcga_df %>% dplyr::filter(as.character(stage) == st)
  if (nrow(d) < min_stage_n || sum(d$os_event, na.rm = TRUE) < min_stage_e) next
  fit <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_score, d)
  fitb <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_group, d)
  stage_tbl[[st]] <- dplyr::bind_rows(
    if (!is.null(fit)) cox_extract(fit, "^risk_score$") %>% dplyr::mutate(stage_subset = st, type = "continuous", n = nrow(d), events = sum(d$os_event)) else NULL,
    if (!is.null(fitb)) cox_extract(fitb, "^risk_group") %>% dplyr::mutate(stage_subset = st, type = "binary", n = nrow(d), events = sum(d$os_event)) else NULL
  )
}
for (sg in c("Early (I-II)", "Late (III-IV)")) {
  d <- tcga_df %>% dplyr::filter(stage_group == sg)
  if (nrow(d) < min_stage_n || sum(d$os_event, na.rm = TRUE) < min_stage_e) next
  fit <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_score, d)
  fitb <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_group, d)
  stage_tbl[[sg]] <- dplyr::bind_rows(
    if (!is.null(fit)) cox_extract(fit, "^risk_score$") %>% dplyr::mutate(stage_subset = sg, type = "continuous", n = nrow(d), events = sum(d$os_event)) else NULL,
    if (!is.null(fitb)) cox_extract(fitb, "^risk_group") %>% dplyr::mutate(stage_subset = sg, type = "binary", n = nrow(d), events = sum(d$os_event)) else NULL
  )
}
stage_tbl <- dplyr::bind_rows(stage_tbl)
safe_write_tsv(stage_tbl, file.path(tab_dir, "table_tcga_stage_stratified.tsv"))

# ----------------------------
# multivariable models
# ----------------------------
candidate_covars <- c("age","sex","stage_group","t_stage","n_stage","m_stage")
available_covars <- candidate_covars[vapply(candidate_covars, function(v) v %in% names(tcga_df) && sum(!is.na(tcga_df[[v]])) >= max(20, ceiling(0.25 * nrow(tcga_df))), logical(1))]
clinic_formula <- build_formula(FALSE, available_covars)
full_formula <- build_formula(TRUE, available_covars)

clinic_df <- tcga_df[complete_cases_for_formula(tcga_df, clinic_formula), , drop = FALSE]
full_df <- tcga_df[complete_cases_for_formula(tcga_df, full_formula), , drop = FALSE]

clinic_covars2 <- usable_covariates(clinic_df, setdiff(all.vars(clinic_formula), c("os_time", "os_event")))
full_covars2 <- usable_covariates(full_df, setdiff(all.vars(full_formula), c("os_time", "os_event", "risk_score")))
clinic_formula <- build_formula(FALSE, clinic_covars2)
full_formula <- build_formula(TRUE, full_covars2)
clinic_df <- clinic_df[complete_cases_for_formula(clinic_df, clinic_formula), , drop = FALSE]
full_df <- full_df[complete_cases_for_formula(full_df, full_formula), , drop = FALSE]

fit_clinic <- safe_coxph(clinic_formula, clinic_df)
fit_full <- safe_coxph(full_formula, full_df)
fit_risk_only <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_score, tcga_df)

multi_tbl <- dplyr::bind_rows(
  if (!is.null(fit_risk_only)) tidy_cox_table(fit_risk_only, "Risk only") else NULL,
  if (!is.null(fit_clinic)) tidy_cox_table(fit_clinic, "Clinicopathologic only") else NULL,
  if (!is.null(fit_full)) tidy_cox_table(fit_full, "Risk + clinicopathologic") else NULL
)
safe_write_tsv(multi_tbl, file.path(tab_dir, "table_multivariable_clinicopathologic_models.tsv"))

compare_df <- full_df
fit_risk_cmp <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_score, compare_df)
fit_clinic_cmp <- safe_coxph(clinic_formula, compare_df)
fit_full_cmp <- safe_coxph(full_formula, compare_df)

model_compare_tbl <- tibble::tibble(
  model = c("Risk only", "Clinicopathologic only", "Risk + clinicopathologic"),
  n = c(if (!is.null(fit_risk_cmp)) nobs(fit_risk_cmp) else NA_integer_,
        if (!is.null(fit_clinic_cmp)) nobs(fit_clinic_cmp) else NA_integer_,
        if (!is.null(fit_full_cmp)) nobs(fit_full_cmp) else NA_integer_),
  cindex = c(if (!is.null(fit_risk_cmp)) compute_cindex(compare_df, formula(fit_risk_cmp)) else NA_real_,
             if (!is.null(fit_clinic_cmp)) compute_cindex(compare_df, formula(fit_clinic_cmp)) else NA_real_,
             if (!is.null(fit_full_cmp)) compute_cindex(compare_df, formula(fit_full_cmp)) else NA_real_),
  AIC = c(if (!is.null(fit_risk_cmp)) AIC(fit_risk_cmp) else NA_real_,
          if (!is.null(fit_clinic_cmp)) AIC(fit_clinic_cmp) else NA_real_,
          if (!is.null(fit_full_cmp)) AIC(fit_full_cmp) else NA_real_),
  BIC = c(if (!is.null(fit_risk_cmp)) BIC(fit_risk_cmp) else NA_real_,
          if (!is.null(fit_clinic_cmp)) BIC(fit_clinic_cmp) else NA_real_,
          if (!is.null(fit_full_cmp)) BIC(fit_full_cmp) else NA_real_),
  lrtest_vs_clinic = c(NA_real_, NA_real_, NA_real_)
)
if (!is.null(fit_clinic_cmp) && !is.null(fit_full_cmp)) {
  lrt <- tryCatch(anova(fit_clinic_cmp, fit_full_cmp, test = "LRT"), error = function(e) NULL)
  improve_p <- NA_real_
  if (!is.null(lrt) && nrow(lrt) >= 2) {
    pcol <- grep("Pr(>|Chi|)", names(lrt), fixed = TRUE, value = TRUE)
    if (length(pcol) > 0) improve_p <- safe_scalar_numeric(lrt[[pcol[1]]][2])
  }
  idx_full <- which(model_compare_tbl$model == "Risk + clinicopathologic")
  if (length(idx_full) == 1) model_compare_tbl$lrtest_vs_clinic[idx_full] <- improve_p
}
safe_write_tsv(model_compare_tbl, file.path(tab_dir, "table_model_comparison.tsv"))

# ----------------------------
# bootstrap and PH assumptions
# ----------------------------
boot_risk <- bootstrap_optimism(tcga_df, survival::Surv(os_time, os_event) ~ risk_score, B = bootstrap_B)
boot_full <- if (!is.null(fit_full)) bootstrap_optimism(full_df, formula(fit_full), B = bootstrap_B) else tibble::tibble(apparent_cindex=NA_real_, optimism=NA_real_, corrected_cindex=NA_real_)
boot_tbl <- dplyr::bind_rows(
  boot_risk %>% dplyr::mutate(model = "Risk only"),
  boot_full %>% dplyr::mutate(model = "Risk + clinicopathologic")
)
safe_write_tsv(boot_tbl, file.path(tab_dir, "table_bootstrap_optimism.tsv"))

ph_tbl <- dplyr::bind_rows(
  if (!is.null(fit_risk_only)) extract_ph_table(fit_risk_only, "Risk only") else NULL,
  if (!is.null(fit_full)) extract_ph_table(fit_full, "Risk + clinicopathologic") else NULL
)
safe_write_tsv(ph_tbl, file.path(tab_dir, "table_proportional_hazards_checks.tsv"))

# ----------------------------
# subgroup analyses
# ----------------------------
subgroup_tbl <- dplyr::bind_rows(
  subgroup_cox(tcga_df, "sex", min_stage_n, min_stage_e),
  subgroup_cox(tcga_df, "stage_group", min_stage_n, min_stage_e),
  subgroup_cox(tcga_df, "t_stage", min_stage_n, min_stage_e),
  subgroup_cox(tcga_df, "n_stage", min_stage_n, min_stage_e)
)
if ("age" %in% names(tcga_df) && sum(!is.na(tcga_df$age)) >= 30) {
  tcga_df$age_group <- factor(ifelse(tcga_df$age >= stats::median(tcga_df$age, na.rm = TRUE), "Older", "Younger"))
  subgroup_tbl <- dplyr::bind_rows(subgroup_tbl, subgroup_cox(tcga_df, "age_group", min_stage_n, min_stage_e))
}
safe_write_tsv(subgroup_tbl, file.path(tab_dir, "table_tcga_subgroup_forest.tsv"))

# ----------------------------
# optional modules
# ----------------------------
maybe_install_optional(c("rms","rmda","timeROC"))
nomogram_tbl <- tibble::tibble()
calibration_tbl <- tibble::tibble()
dca_tbl <- tibble::tibble()
optional_notes <- c()

if (!is.null(fit_full) && optional_pkg("rms")) {
  suppressPackageStartupMessages(library(rms))
  dd <- datadist(full_df)
  oldopt <- options(datadist = "dd")
  on.exit(options(oldopt), add = TRUE)

  if (requireNamespace("Hmisc", quietly = TRUE)) {
    try(Hmisc::units(full_df$os_time) <- "Year", silent = TRUE)
  } else {
    attr(full_df$os_time, "units") <- "Year"
  }
  cph_fit <- tryCatch(rms::cph(formula(fit_full), data = full_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 3), error = function(e) NULL)
  if (!is.null(cph_fit)) {
    surv_fun <- rms::Survival(cph_fit)
    nom <- tryCatch(rms::nomogram(
      cph_fit,
      fun = list(function(lp) surv_fun(1, lp), function(lp) surv_fun(3, lp), function(lp) surv_fun(5, lp)),
      funlabel = c("1-year OS","3-year OS","5-year OS")
    ), error = function(e) NULL)
    if (!is.null(nom)) {
      grDevices::pdf(file.path(fig_dir, "FIG04_nomogram_tcga_full_model.pdf"), width = 12, height = 8)
      plot(nom, xfrac = 0.45)
      grDevices::dev.off()
      grDevices::png(file.path(fig_dir, "FIG04_nomogram_tcga_full_model.png"), width = 1800, height = 1200, res = 220)
      plot(nom, xfrac = 0.45)
      grDevices::dev.off()
      nomogram_tbl <- tibble::tibble(model = "Risk + clinicopathologic", horizons = "1,3,5 years", note = "Nomogram generated with rms::nomogram")
    }
    cal_list <- list()
    for (u in time_horizons) {
      cal <- tryCatch(rms::calibrate(cph_fit, method = "boot", u = u, B = cal_B), error = function(e) NULL)
      if (is.null(cal)) next
      cdf <- as.data.frame(cal)
      names(cdf) <- make.names(names(cdf))
      xcol <- names(cdf)[grep("predy|mean.predicted", names(cdf), ignore.case = TRUE)][1]
      ycol <- names(cdf)[grep("km|calibrated.corrected|mean.observed", names(cdf), ignore.case = TRUE)][1]
      if (!is.na(xcol) && !is.na(ycol)) {
        cal_out <- tibble::tibble(time_year = u, predicted = cdf[[xcol]], observed = cdf[[ycol]])
        cal_list[[as.character(u)]] <- cal_out
        pcal <- ggplot2::ggplot(cal_out, ggplot2::aes(predicted, observed)) +
          ggplot2::geom_point() + ggplot2::geom_line() +
          ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "firebrick") +
          ggplot2::coord_equal() + ggplot2::theme_bw(base_size = 12) +
          ggplot2::labs(title = paste0("Calibration at ", u, "-year OS"), x = "Predicted survival", y = "Observed survival")
        save_plot(pcal, file.path(fig_dir, paste0("FIG05_calibration_", gsub("\\.", "p", as.character(u)), "yr")), 6, 6)
      }
    }
    calibration_tbl <- dplyr::bind_rows(cal_list)
  }
} else {
  optional_notes <- c(optional_notes, "rms package unavailable: nomogram/calibration skipped.")
}

if (!is.null(fit_full) && optional_pkg("rmda")) {
  suppressPackageStartupMessages(library(rmda))
  dca_list <- list()
  for (u in time_horizons) {
    d <- full_df %>% dplyr::mutate(event_by_u = ifelse(os_event == 1 & os_time <= u, 1, 0))
    if (sum(d$event_by_u, na.rm = TRUE) < 10) next
    fit_r <- stats::glm(event_by_u ~ risk_score, data = d, family = stats::binomial())
    d$risk_only_prob <- stats::predict(fit_r, type = "response")
    fit_f <- stats::glm(as.formula(paste("event_by_u ~", paste(setdiff(all.vars(full_formula), c("os_time","os_event")), collapse = " + "))), data = d, family = stats::binomial())
    d$full_model_prob <- stats::predict(fit_f, type = "response")
    formulas <- c("event_by_u ~ risk_only_prob", "event_by_u ~ full_model_prob")
    labels <- c("Risk only", "Risk + clinicopathologic")
    if (!is.null(fit_clinic)) {
      fit_c <- tryCatch(stats::glm(as.formula(paste("event_by_u ~", paste(setdiff(all.vars(clinic_formula), c("os_time","os_event")), collapse = " + "))), data = d, family = stats::binomial()), error = function(e) NULL)
      if (!is.null(fit_c)) {
        d$clinic_prob <- stats::predict(fit_c, type = "response")
        formulas <- c("event_by_u ~ risk_only_prob", "event_by_u ~ clinic_prob", "event_by_u ~ full_model_prob")
        labels <- c("Risk only", "Clinicopathologic only", "Risk + clinicopathologic")
      }
    }
    dc <- tryCatch(rmda::decision_curve(formulas = formulas, data = d, fitted.risk = TRUE, confidence.intervals = FALSE, thresholds = seq(0.01, 0.60, by = 0.01), study.design = "cohort"), error = function(e) NULL)
    if (!is.null(dc)) {
      grDevices::pdf(file.path(fig_dir, paste0("FIG06_dca_", gsub("\\.", "p", as.character(u)), "yr.pdf")), width = 8, height = 6)
      plot_decision_curve(dc, curve.names = labels, xlim = c(0, 0.6), legend.position = "bottomright", standardize = FALSE)
      grDevices::dev.off()
      grDevices::png(file.path(fig_dir, paste0("FIG06_dca_", gsub("\\.", "p", as.character(u)), "yr.png")), width = 1600, height = 1200, res = 220)
      plot_decision_curve(dc, curve.names = labels, xlim = c(0, 0.6), legend.position = "bottomright", standardize = FALSE)
      grDevices::dev.off()
      if (!is.null(dc$derived.data)) {
        tmp <- tibble::as_tibble(dc$derived.data)
        tmp$time_year <- u
        dca_list[[as.character(u)]] <- tmp
      }
    }
  }
  dca_tbl <- dplyr::bind_rows(dca_list)
} else {
  optional_notes <- c(optional_notes, "rmda package unavailable: DCA skipped.")
}

safe_write_tsv(nomogram_tbl, file.path(tab_dir, "table_nomogram_summary.tsv"))
safe_write_tsv(calibration_tbl, file.path(tab_dir, "table_calibration_curves.tsv"))
safe_write_tsv(dca_tbl, file.path(tab_dir, "table_dca_net_benefit.tsv"))

# ----------------------------
# time-dependent AUC
# ----------------------------
time_auc_tbl <- tibble::tibble()
if (optional_pkg("timeROC") && !is.null(fit_full)) {
  suppressPackageStartupMessages(library(timeROC))
  d <- full_df
  use_times <- time_horizons[time_horizons < max(d$os_time, na.rm = TRUE)]
  if (length(use_times) > 0) {
    time_auc_tbl <- time_roc_by_models(d, list(
      "Risk only" = fit_risk_cmp,
      "Clinicopathologic only" = fit_clinic_cmp,
      "Risk + clinicopathologic" = fit_full_cmp
    ), use_times)
    safe_write_tsv(time_auc_tbl, file.path(tab_dir, "table_time_dependent_auc_models.tsv"))
    if (nrow(time_auc_tbl) > 0) {
      p_auc <- ggplot2::ggplot(time_auc_tbl, ggplot2::aes(time_year, AUC, color = model)) +
        ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2) +
        ggplot2::ylim(0.5, 1) + ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(title = "Time-dependent AUC across prognostic models", x = "Time (years)", y = "AUC", color = NULL)
      save_plot(p_auc, file.path(fig_dir, "FIG07_time_dependent_auc_models"), 7, 5)
    }
  }
} else {
  optional_notes <- c(optional_notes, "timeROC package unavailable: time-dependent AUC skipped.")
}

# ----------------------------
# external validation
# ----------------------------
external_results <- list()
external_auc_tbl <- list()
external_stage_assoc <- list()

for (coh in names(ext_paths)) {
  if (tcga_only) next
  p <- ext_paths[[coh]]
  if (is.na(p$expr) || is.na(p$surv)) next
  msg("Loading external cohort: %s", coh)
  expr <- tryCatch(read_expr_matrix(p$expr, coh), error = function(e) NULL)
  surv <- tryCatch(read_survival_clinical(p$surv, p$clin, tcga = FALSE), error = function(e) NULL)
  if (is.null(expr) || is.null(surv)) next
  common <- intersect(rownames(expr), surv$sample_id)
  if (length(common) < 20) next
  expr <- expr[common, , drop = FALSE]
  surv <- surv %>% dplyr::filter(sample_id %in% common) %>% dplyr::arrange(match(sample_id, rownames(expr)))
  score <- tryCatch(score_signature(expr, locked_tbl %>% dplyr::select(gene, weight)), error = function(e) NULL)
  if (is.null(score)) next
  df <- surv %>%
    dplyr::mutate(risk_score = score[sample_id],
                  risk_score_z = standardize_score(risk_score),
                  risk_group = factor(ifelse(risk_score >= stats::median(risk_score, na.rm = TRUE), "High", "Low"), levels = c("Low","High")))
  fitc <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_score, df)
  fitb <- safe_coxph(survival::Surv(os_time, os_event) ~ risk_group, df)
  res <- dplyr::bind_rows(
    if (!is.null(fitc)) cox_extract(fitc, "^risk_score$") %>% dplyr::mutate(cohort = coh, type = "continuous", n = nrow(df), events = sum(df$os_event, na.rm = TRUE), cindex = compute_cindex(df, survival::Surv(os_time, os_event) ~ risk_score)) else NULL,
    if (!is.null(fitb)) cox_extract(fitb, "^risk_group") %>% dplyr::mutate(cohort = coh, type = "binary", n = nrow(df), events = sum(df$os_event, na.rm = TRUE), cindex = compute_cindex(df, survival::Surv(os_time, os_event) ~ risk_score)) else NULL
  )
  external_results[[coh]] <- res
  safe_write_tsv(df, file.path(tab_dir, paste0("external_", coh, "_risk_scores.tsv")))
  external_stage_assoc[[coh]] <- make_stage_assoc_table(df) %>% dplyr::mutate(cohort = coh)
  if (optional_pkg("timeROC") && !is.null(fitc)) {
    use_times <- time_horizons[time_horizons < max(df$os_time, na.rm = TRUE)]
    if (length(use_times) > 0 && sum(df$os_event, na.rm = TRUE) >= 10) {
      roc <- tryCatch(timeROC::timeROC(T = df$os_time, delta = df$os_event, marker = df$risk_score, cause = 1, times = use_times, iid = TRUE), error = function(e) NULL)
      if (!is.null(roc)) external_auc_tbl[[coh]] <- tibble::tibble(cohort = coh, time_year = use_times, AUC = roc$AUC)
    }
  }
}
external_tbl <- dplyr::bind_rows(external_results)
safe_write_tsv(external_tbl, file.path(tab_dir, "table_external_validation.tsv"))
safe_write_tsv(dplyr::bind_rows(external_auc_tbl), file.path(tab_dir, "table_external_time_auc.tsv"))
safe_write_tsv(dplyr::bind_rows(external_stage_assoc), file.path(tab_dir, "table_external_stage_association.tsv"))

meta_tbl <- if (nrow(external_tbl) > 0) {
  dplyr::bind_rows(
    external_tbl %>% dplyr::filter(type == "continuous") %>% meta_random_hr() %>% dplyr::mutate(type = "continuous"),
    external_tbl %>% dplyr::filter(type == "binary") %>% meta_random_hr() %>% dplyr::mutate(type = "binary")
  )
} else tibble::tibble()
safe_write_tsv(meta_tbl, file.path(tab_dir, "table_external_meta_analysis.tsv"))

if (nrow(external_tbl %>% dplyr::filter(type == "continuous")) >= 3) {
  loo <- lapply(unique(external_tbl$cohort[external_tbl$type == "continuous"]), function(drop_coh) {
    tmp <- external_tbl %>% dplyr::filter(type == "continuous", cohort != drop_coh)
    mt <- meta_random_hr(tmp)
    if (is.null(mt)) return(NULL)
    mt$left_out <- drop_coh
    mt
  }) %>% dplyr::bind_rows()
  safe_write_tsv(loo, file.path(tab_dir, "table_leave_one_cohort_out_meta.tsv"))
}

# ----------------------------
# plots
# ----------------------------
if (nrow(stage_tbl) > 0) {
  p1 <- ggplot2::ggplot(stage_tbl %>% dplyr::filter(type == "continuous"),
                        ggplot2::aes(x = stage_subset, y = HR, ymin = lower95, ymax = upper95)) +
    ggplot2::geom_pointrange() +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "firebrick") +
    ggplot2::coord_flip() + ggplot2::scale_y_log10() + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "TCGA-LUAD stage-stratified validation", x = NULL, y = "Hazard ratio (log scale)")
  save_plot(p1, file.path(fig_dir, "FIG01_tcga_stage_stratified_forest"), 8, 5)
}

p2 <- km_plot_basic(tcga_df, "TCGA-LUAD overall survival by locked signature")
save_plot(p2, file.path(fig_dir, "FIG02_tcga_km_locked_signature"), 8, 6)

if (nrow(external_tbl) > 0) {
  p3 <- ggplot2::ggplot(external_tbl %>% dplyr::filter(type == "continuous"),
                        ggplot2::aes(x = cohort, y = HR, ymin = lower95, ymax = upper95)) +
    ggplot2::geom_pointrange() +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "firebrick") +
    ggplot2::coord_flip() + ggplot2::scale_y_log10() + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "External GEO validation of TCGA-locked signature", x = NULL, y = "Hazard ratio (log scale)")
  save_plot(p3, file.path(fig_dir, "FIG03_external_validation_forest"), 8, 4.8)
}

if (nrow(model_compare_tbl) > 0) {
  p_model <- ggplot2::ggplot(model_compare_tbl, ggplot2::aes(model, cindex)) +
    ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "Discrimination of prognostic models", x = NULL, y = "C-index")
  save_plot(p_model, file.path(fig_dir, "FIG08_model_cindex_comparison"), 7, 4.5)
}

p_gene <- base_gene_tbl %>% dplyr::slice_head(n = 20) %>%
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>%
  ggplot2::ggplot(ggplot2::aes(gene, abs_weight, fill = dir_beta > 0)) +
  ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(title = "Top 20 weighted signature genes", x = NULL, y = "Absolute weight", fill = "Positive beta")
save_plot(p_gene, file.path(fig_dir, "FIG09_top_signature_gene_weights"), 7.5, 6)

p_risk_stage <- ggplot2::ggplot(tcga_df, ggplot2::aes(stage_group, risk_score, fill = stage_group)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  ggplot2::theme_bw(base_size = 12) + ggplot2::labs(title = "Risk score distribution across stage groups", x = NULL, y = "Risk score") +
  ggplot2::guides(fill = "none")
save_plot(p_risk_stage, file.path(fig_dir, "FIG10_risk_score_by_stage_group"), 7, 5)

for (st in c("I","II","III","IV")) {
  d <- tcga_df %>% dplyr::filter(as.character(stage) == st)
  if (nrow(d) < min_stage_n || sum(d$os_event) < min_stage_e || length(unique(d$risk_group)) < 2) next
  pk <- km_plot_basic(d, paste0("TCGA-LUAD stage ", st, " survival by locked signature"))
  save_plot(pk, file.path(fig_dir, paste0("FIG11_KM_stage_", st)), 8, 6)
}

if (nrow(subgroup_tbl) > 0) {
  p_sub <- ggplot2::ggplot(subgroup_tbl, ggplot2::aes(x = paste(subgroup, level, sep = ": "), y = HR, ymin = lower95, ymax = upper95)) +
    ggplot2::geom_pointrange() + ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "firebrick") +
    ggplot2::coord_flip() + ggplot2::scale_y_log10() + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "Subgroup consistency of the continuous risk score", x = NULL, y = "Hazard ratio (log scale)")
  save_plot(p_sub, file.path(fig_dir, "FIG12_tcga_subgroup_forest"), 8, 6)
}

if (nrow(ph_tbl) > 0) {
  p_ph <- ggplot2::ggplot(ph_tbl %>% dplyr::filter(term != "GLOBAL"), ggplot2::aes(term, -log10(p_value), fill = model)) +
    ggplot2::geom_col(position = "dodge") + ggplot2::coord_flip() + ggplot2::theme_bw(base_size = 12) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, color = "firebrick") +
    ggplot2::labs(title = "Proportional-hazards diagnostics", x = NULL, y = expression(-log[10](p)))
  save_plot(p_ph, file.path(fig_dir, "FIG13_ph_diagnostics"), 8, 5)
}

# ----------------------------
# workbook export
# ----------------------------
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  sheets <- list(
    tcga_cohort_qc = cohort_qc_tbl,
    tcga_gene_level = base_gene_tbl,
    locked_signature = locked_tbl,
    tcga_risk_scores = tcga_df,
    clin_by_risk = clin_summary_tbl,
    stage_assoc = make_stage_assoc_table(tcga_df),
    tcga_overall = overall_tbl,
    tcga_stage_stratified = stage_tbl,
    multivariable_models = multi_tbl,
    model_comparison = model_compare_tbl,
    bootstrap = boot_tbl,
    ph_checks = ph_tbl,
    subgroup_forest = subgroup_tbl,
    nomogram = nomogram_tbl,
    calibration = calibration_tbl,
    dca = dca_tbl,
    time_auc = time_auc_tbl,
    external_validation = external_tbl,
    external_time_auc = dplyr::bind_rows(external_auc_tbl),
    external_meta = meta_tbl
  )
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, substr(gsub("[^A-Za-z0-9_]", "_", nm), 1, 31))
    openxlsx::writeData(wb, substr(gsub("[^A-Za-z0-9_]", "_", nm), 1, 31), as.data.frame(sheets[[nm]]))
  }
  openxlsx::saveWorkbook(wb, file.path(tab_dir, "stage_stratified_tcga_luad_q1_enhanced.xlsx"), overwrite = TRUE)
}

# ----------------------------
# output manifest and summary
# ----------------------------
out_files <- list.files(outdir, recursive = TRUE, full.names = TRUE)
out_files <- out_files[file.info(out_files)$isdir %in% FALSE]
out_manifest <- tibble::tibble(
  relative_path = sub(paste0("^", normalizePath(outdir, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(out_files, winslash = "/", mustWork = FALSE)),
  bytes = file.info(out_files)$size
) %>% dplyr::arrange(relative_path)
safe_write_tsv(out_manifest, file.path(log_dir, "output_manifest.tsv"))

summary_lines <- c(
  "Enhanced stage-stratified TCGA-LUAD validation and multi-cohort transportability analysis",
  sprintf("Project root: %s", project_root),
  sprintf("Discovery cohort: TCGA-LUAD only (%d patients, %d events)", nrow(tcga_df), sum(tcga_df$os_event, na.rm = TRUE)),
  sprintf("Candidate gene file: %s", genes_file),
  sprintf("Locked signature size: %d genes", nrow(locked_tbl)),
  sprintf("Clinical covariates retained: %s", paste(available_covars, collapse = ", ")),
  sprintf("Stage interaction p-value: %s", format(interaction_p, digits = 4)),
  sprintf("Bootstrap optimism-corrected c-index (risk only): %s", format(boot_tbl$corrected_cindex[boot_tbl$model == "Risk only"][1], digits = 4)),
  sprintf("Bootstrap optimism-corrected c-index (full model): %s", format(boot_tbl$corrected_cindex[boot_tbl$model == "Risk + clinicopathologic"][1], digits = 4)),
  sprintf("External validation enabled: %s", ifelse(tcga_only, "No (--tcga-only)", "Yes")),
  sprintf("External cohorts resolved: %s", paste(names(ext_paths)[vapply(ext_paths, function(z) !is.na(z$expr) && !is.na(z$surv), logical(1))], collapse = ", ")),
  "Added Q1-grade features: PH diagnostics, subgroup forest, cross-model time-dependent AUC, stage-association testing, workbook export, and stricter multivariable comparison on a common analysis set.",
  "GSE81089 intentionally excluded because it belongs to WGCNA discovery, not prognostic validation in this analysis.",
  if (length(optional_notes) > 0) paste(optional_notes, collapse = " ") else "All optional modules executed successfully."
)
writeLines(summary_lines, file.path(log_dir, "analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

msg("Analysis completed. Results written to: %s", outdir)
