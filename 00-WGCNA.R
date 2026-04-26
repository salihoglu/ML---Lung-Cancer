#!/usr/bin/env Rscript

############################################################
#Rana Salihoglu
#WGCNA Analysis 
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
  library(matrixStats)
})
suppressWarnings({
  if (requireNamespace("survival", quietly = TRUE)) library(survival)
})

enableWGCNAThreads()
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
set.seed(2025)

# -------------------- PATHS --------------------
BASE_DIR   <- "GEO_PREP/GSE81089"
EXPR_FILE  <- file.path(BASE_DIR, "X_rna_symbol.csv")
SURV_FILE  <- file.path(BASE_DIR, "survival.csv")
PHENO_FILE <- file.path(BASE_DIR, "GSE81089_pheno_expanded.csv")

SURV_SAMPLE_COL <- "sample"
SURV_TIME_COL   <- "OS.time"
SURV_EVENT_COL  <- "OS"

FERR_DRIVER <- "ferrdb_driver.csv"
FERR_SUPP   <- "ferrdb_suppressor.csv"
FERR_MARKER <- "ferrdb_marker.csv"   # optional

OUTDIR <- file.path(BASE_DIR, "WGCNA_GSE81089_out_v12_signature_fix")
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")
CYTODIR <- file.path(OUTDIR, "cytoscape")
LOGDIR <- file.path(OUTDIR, "logs")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CYTODIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

# -------------------- SETTINGS --------------------
MIN_SAMPLES_EXPR <- 20
MAX_GENES_MAD <- 8000

NETWORK_TYPE <- "signed"
TOM_TYPE <- "signed"
COR_TYPE <- "bicor"          # v13 FIX: bicor (outlier-robust); v12'de pearson idi
POWER_RANGE <- 1:20
R2_TARGET <- 0.80            # v13 FIX: 0.85 → 0.80 (signed network için daha güvenli)
FORCE_POWER <- NA_integer_   # örn: 6 — elle belirlemek için
MIN_MODULE_SIZE <- 30
MERGE_CUT_HEIGHT <- 0.25

# ---- v13: dual-branch kontrolleri (artık gerçekten kullanılıyor) ----
ANALYSIS_MODE <- "both"      # "conservative" | "discovery" | "both"
EXPORT_MODULE_TOPN <- 200    # discovery branch için top N gen (MM+GS skoruna göre)
MM_MIN_ABS <- 0.30           # discovery: |Module Membership| minimum eşiği
GS_MIN_ABS <- 0.20           # discovery: |Gene Significance to FERRO_TRAIT| minimum eşiği
UNIV_P_CUTOFF <- 0.05        # discovery: univariate Cox p-değeri eşiği
MIN_FERRO_IN_MODULE <- 3     # v13 YENİ: modülde minimum FerrDb geni (biyolojik güvence)
MAX_MODULE_GENES_ML <- 300   # v13 YENİ: ML pipeline'a gönderilecek max gen sayısı
                              # Büyük modüllerde (>300 gen) hub score sıralamasıyla daralt
HUB_SCORE_TOPN <- 150        # v13 YENİ: conservative branch'te top N hub gen (MM×GS skoru)

DEEP_SPLIT <- 2

DO_SAMPLE_OUTLIER_QC <- TRUE
SAMPLE_CUT_HEIGHT <- NULL

MODULE_PICK_TRAIT <- "FERRO_TRAIT"
MODULE_RHO_MIN <- 0.20
MODULE_FDR_MAX <- 0.10

ABS_MM_MIN <- 0.60
ABS_GS_MIN <- 0.30
GS_FDR_MAX <- 0.05
MM_P_MAX <- 0.05

MIN_SAMPLES_COX <- 30
P_FLOOR <- 1e-300
TOP_GENES_FOR_ML <- 200

CYTO_MAX_EDGES <- 30000
CYTO_WEIGHT_Q <- 0.95

MIN_GENES_FOR_WGCNA <- 200
MIN_GENES_FOR_MODULES <- 60
MIN_SAMPLES_FOR_WGCNA <- 20

# -------------------- HELPERS --------------------
p_floor <- function(p, floor = P_FLOOR) {
  p <- as.numeric(p); p[!is.finite(p)] <- 1; pmax(p, floor)
}

safe_fwrite <- function(x, file, sep = "\t") {
  if (is.null(x)) x <- data.frame()
  if (!is.data.frame(x)) x <- as.data.frame(x)
  fwrite(x, file = file, sep = sep)
}

save_pdf_png <- function(stem, width, height, plot_fun) {
  pdf(file.path(FIGDIR, paste0(stem, ".pdf")), width = width, height = height)
  plot_fun(); dev.off()
  png(file.path(FIGDIR, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300)
  plot_fun(); dev.off()
}

save_placeholder_figure <- function(stem, msg, width = 8, height = 5) {
  save_pdf_png(stem, width, height, function() {
    plot.new(); title(main = stem); text(0.5, 0.55, msg, cex = 1)
  })
}

normalize_sample_id <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s[s %in% c("", "na", "n/a", "nan", "none")] <- NA_character_
  s <- sub("_[0-9]+$", "", s)
  s <- gsub("rna[- _\\.]*seq", "", s)
  s <- gsub("[^a-z0-9]+", "", s)
  s <- sub("^x(?=[0-9])", "", s, perl = TRUE)
  s
}

parse_stage_token <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s[s %in% c("", "NA", "N/A", "NONE", "NULL", "UNKNOWN", "--", "-", "NOT AVAILABLE")] <- NA
  s <- gsub("STAGE", "", s)
  s <- gsub("[^A-Z0-9]", "", s)
  
  out <- rep(NA_real_, length(s))
  num <- regmatches(s, regexpr("^[1-4]", s))
  has_num <- nchar(num) > 0
  out[has_num] <- as.numeric(num[has_num])
  
  idx <- which(is.na(out) & !is.na(s))
  if (length(idx)) {
    ss <- s[idx]
    out[idx[grepl("^IV", ss)]]  <- 4
    out[idx[grepl("^III", ss)]] <- 3
    out[idx[grepl("^II", ss)]]  <- 2
    out[idx[grepl("^I", ss)]]   <- 1
    out[idx[grepl("^T1", ss)]]  <- 1
    out[idx[grepl("^T2", ss)]]  <- 2
    out[idx[grepl("^T3", ss)]]  <- 3
    out[idx[grepl("^T4", ss)]]  <- 4
  }
  out
}

safe_cor_test_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20) return(list(rho = NA_real_, p = 1, n = sum(ok)))
  rho <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  p <- tryCatch(cor.test(x[ok], y[ok], method = "spearman")$p.value, error = function(e) 1)
  list(rho = as.numeric(rho), p = p_floor(p), n = sum(ok))
}

safe_cor_test_pearson <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20) return(list(r = NA_real_, p = 1, n = sum(ok)))
  r <- suppressWarnings(cor(x[ok], y[ok], method = "pearson"))
  p <- tryCatch(cor.test(x[ok], y[ok], method = "pearson")$p.value, error = function(e) 1)
  list(r = as.numeric(r), p = p_floor(p), n = sum(ok))
}

median_impute_by_gene <- function(mat) {
  if (!anyNA(mat)) return(mat)
  for (j in seq_len(ncol(mat))) {
    x <- mat[, j]
    if (anyNA(x)) {
      md <- median(x, na.rm = TRUE); if (!is.finite(md)) md <- 0
      x[is.na(x)] <- md
      mat[, j] <- x
    }
  }
  mat
}

make_cor_options <- function(cor_type = "pearson") {
  if (tolower(cor_type) == "pearson") {
    list(use = "pairwise.complete.obs")
  } else {
    list(use = "pairwise.complete.obs", maxPOutliers = 0.1)
  }
}

# ---- CRITICAL: signature-compatible wrappers ----
call_pickSoftThreshold_compat <- function(datExpr, powerVector, networkType, cor_type = "pearson", verbose = 5) {
  f <- get("pickSoftThreshold", asNamespace("WGCNA"))
  fn <- names(formals(f))
  corFnc <- if (tolower(cor_type) == "pearson") "cor" else "bicor"
  corOptions <- make_cor_options(cor_type)
  
  # A) positional first arg + corFnc/corOptions if supported
  argsA <- list(datExpr, powerVector = powerVector, networkType = networkType, verbose = verbose)
  if ("corFnc" %in% fn) argsA$corFnc <- corFnc
  if ("corOptions" %in% fn) argsA$corOptions <- corOptions
  
  outA <- tryCatch(do.call(f, argsA), error = function(e) e)
  if (!inherits(outA, "error")) return(outA)
  message("[WARN] pickSoftThreshold mode-A failed: ", outA$message)
  
  # B) try corType if supported
  argsB <- list(datExpr, powerVector = powerVector, networkType = networkType, verbose = verbose)
  if ("corType" %in% fn) argsB$corType <- tolower(cor_type)
  outB <- tryCatch(do.call(f, argsB), error = function(e) e)
  if (!inherits(outB, "error")) return(outB)
  message("[WARN] pickSoftThreshold mode-B failed: ", outB$message)
  
  # C) minimal fallback
  argsC <- list(datExpr, powerVector = powerVector, networkType = networkType, verbose = verbose)
  outC <- tryCatch(do.call(f, argsC), error = function(e) e)
  if (!inherits(outC, "error")) return(outC)
  
  stop("pickSoftThreshold failed in all compatibility modes: ", outC$message)
}

call_adjacency_compat <- function(datExpr, power, networkType, cor_type = "pearson") {
  f <- get("adjacency", asNamespace("WGCNA"))
  fn <- names(formals(f))
  corFnc <- if (tolower(cor_type) == "pearson") "cor" else "bicor"
  corOptions <- make_cor_options(cor_type)
  
  argsA <- list(datExpr, power = power, type = networkType)
  if ("corFnc" %in% fn) argsA$corFnc <- corFnc
  if ("corOptions" %in% fn) argsA$corOptions <- corOptions
  outA <- tryCatch(do.call(f, argsA), error = function(e) e)
  if (!inherits(outA, "error")) return(outA)
  message("[WARN] adjacency mode-A failed: ", outA$message)
  
  argsB <- list(datExpr, power = power, type = networkType)
  outB <- tryCatch(do.call(f, argsB), error = function(e) e)
  if (!inherits(outB, "error")) return(outB)
  
  stop("adjacency failed in compatibility modes: ", outB$message)
}

ssgsea_safe <- function(expr_gene_by_sample, genesets) {
  if (requireNamespace("GSVA", quietly = TRUE)) {
    GSVA <- asNamespace("GSVA")
    if (exists("ssgseaParam", envir = GSVA, inherits = FALSE)) {
      ssgseaParam <- get("ssgseaParam", envir = GSVA)
      gsva <- get("gsva", envir = GSVA)
      return(as.matrix(gsva(ssgseaParam(expr = expr_gene_by_sample, geneSets = genesets))))
    }
    if (exists("gsva", envir = GSVA, inherits = FALSE)) {
      gsva <- get("gsva", envir = GSVA)
      return(as.matrix(gsva(expr = expr_gene_by_sample, gset.idx.list = genesets, method = "ssgsea", verbose = FALSE)))
    }
  }
  X <- expr_gene_by_sample
  rnks <- apply(X, 2, rank, ties.method = "average") / nrow(X)
  out <- matrix(NA_real_, nrow = length(genesets), ncol = ncol(X),
                dimnames = list(names(genesets), colnames(X)))
  for (k in seq_along(genesets)) {
    gs <- intersect(genesets[[k]], rownames(X))
    if (length(gs) < 5) next
    out[k, ] <- colMeans(rnks[gs, , drop = FALSE])
  }
  out
}

run_go_bp <- function(genes, universe) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(data.frame(note = "clusterProfiler/org.Hs.eg.db not installed"))
  }
  eg <- tryCatch(
    clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db),
    error = function(e) NULL
  )
  eg_bg <- tryCatch(
    clusterProfiler::bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db),
    error = function(e) NULL
  )
  if (is.null(eg) || nrow(eg) < 5 || is.null(eg_bg) || nrow(eg_bg) < 20) {
    return(data.frame(note = "Insufficient mapped ENTREZ IDs for GO"))
  }
  ego <- tryCatch(
    clusterProfiler::enrichGO(
      gene = unique(eg$ENTREZID),
      universe = unique(eg_bg$ENTREZID),
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.20, readable = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(ego)) return(data.frame(note = "enrichGO failed"))
  as.data.frame(ego)
}

# -------------------- LOADERS --------------------
load_expr_geo <- function(path_expr, min_samples_expr, max_genes_mad) {
  if (!file.exists(path_expr)) stop("Missing expression file: ", path_expr)
  expr <- fread(path_expr, data.table = FALSE)
  colnames(expr)[1] <- "sample"
  expr$sample <- as.character(expr$sample)
  
  mat <- as.matrix(expr[, -1, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  rownames(mat) <- expr$sample
  
  mat <- mat[, colSums(is.finite(mat)) > 0, drop = FALSE]
  mat <- log2(mat + 1)
  
  ns <- nrow(mat)
  ladder <- unique(pmax(3, c(min_samples_expr, floor(0.8*ns), floor(0.6*ns), floor(0.4*ns), floor(0.25*ns), 3)))
  ladder <- sort(ladder, decreasing = TRUE)
  chosen_thr <- 0; best <- mat
  for (k in ladder) {
    keep <- colSums(mat > 0, na.rm = TRUE) >= k
    tmp <- mat[, keep, drop = FALSE]
    if (ncol(tmp) >= MIN_GENES_FOR_WGCNA) { best <- tmp; chosen_thr <- k; break }
  }
  mat <- best
  
  vars <- matrixStats::colVars(mat, na.rm = TRUE)
  mat <- mat[, is.finite(vars) & vars > 1e-12, drop = FALSE]
  colnames(mat) <- make.unique(colnames(mat), sep = "_")
  
  if (ncol(mat) > max_genes_mad) {
    mads <- apply(mat, 2, mad, na.rm = TRUE)
    mat <- mat[, order(mads, decreasing = TRUE)[seq_len(max_genes_mad)], drop = FALSE]
  }
  
  mat <- median_impute_by_gene(mat)
  mat <- mat[rowSums(is.finite(mat)) == ncol(mat), colSums(is.finite(mat)) == nrow(mat), drop = FALSE]
  
  gsg <- tryCatch(goodSamplesGenes(mat, verbose = 2), error = function(e) NULL)
  if (!is.null(gsg) && (!all(gsg$goodSamples) || !all(gsg$goodGenes))) {
    mat <- mat[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }
  
  if (nrow(mat) < MIN_SAMPLES_FOR_WGCNA) stop("Too few samples after QC: ", nrow(mat))
  if (ncol(mat) < MIN_GENES_FOR_MODULES) stop("Too few genes after QC: ", ncol(mat))
  
  list(mat = mat, threshold_used = chosen_thr)
}

load_survival_geo <- function(path, sample_ids_keep) {
  if (!file.exists(path)) return(NULL)
  s <- fread(path, data.table = FALSE)
  need <- c(SURV_SAMPLE_COL, SURV_TIME_COL, SURV_EVENT_COL)
  if (!all(need %in% colnames(s))) return(NULL)
  
  s <- s[, need, drop = FALSE]
  colnames(s) <- c("sample", "time", "event")
  s$sample <- as.character(s$sample)
  s$time <- suppressWarnings(as.numeric(s$time))
  s$event <- suppressWarnings(as.numeric(s$event))
  s <- s[is.finite(s$time) & is.finite(s$event), , drop = FALSE]
  s$event <- ifelse(s$event > 0, 1, 0)
  
  s2 <- s[s$sample %in% sample_ids_keep, , drop = FALSE]
  if (nrow(s2) == 0) {
    idx <- match(normalize_sample_id(sample_ids_keep), normalize_sample_id(s$sample))
    ok <- which(!is.na(idx))
    if (!length(ok)) return(NULL)
    s2 <- s[idx[ok], , drop = FALSE]
    s2$sample <- sample_ids_keep[ok]
  } else {
    s2 <- s2[match(sample_ids_keep, s2$sample), , drop = FALSE]
  }
  
  s2 <- s2[is.finite(s2$time) & is.finite(s2$event), , drop = FALSE]
  rownames(s2) <- s2$sample
  s2
}

load_stage_geo <- function(pheno_file, sample_ids_keep) {
  if (!file.exists(pheno_file)) return(NULL)
  ph <- fread(pheno_file, data.table = FALSE)
  if (!("sample" %in% colnames(ph))) colnames(ph)[1] <- "sample"
  
  stage_candidates <- c(
    "pathological stage:ch1","pathological.stage",
    "pstage iorii:ch1","pstage.iorii",
    "stage:ch1","stage","Stage:ch1"
  )
  st_col <- stage_candidates[stage_candidates %in% colnames(ph)][1]
  if (is.na(st_col)) {
    stage_like <- colnames(ph)[grepl("stage", tolower(colnames(ph)))]
    st_col <- ifelse(length(stage_like), stage_like[1], NA)
  }
  if (is.na(st_col)) return(NULL)
  
  st <- parse_stage_token(ph[[st_col]])
  x <- data.frame(sample = as.character(ph$sample), stage_num = st, stringsAsFactors = FALSE)
  
  stage <- x$stage_num[match(sample_ids_keep, x$sample)]
  if (sum(is.finite(stage)) < 5) {
    stage <- x$stage_num[match(normalize_sample_id(sample_ids_keep), normalize_sample_id(x$sample))]
  }
  
  med <- median(stage, na.rm = TRUE); if (!is.finite(med)) med <- 2
  imp <- stage; imp[!is.finite(imp)] <- med
  
  data.frame(sample = sample_ids_keep, stage_num = as.numeric(imp),
             stage_missing = as.integer(!is.finite(stage)), stringsAsFactors = FALSE)
}

# -------------------- START --------------------
param_tbl <- data.frame(
  key = c("EXPR_FILE","SURV_FILE","PHENO_FILE","NETWORK_TYPE","TOM_TYPE","COR_TYPE","R2_TARGET","MIN_MODULE_SIZE","MERGE_CUT_HEIGHT","MAX_GENES_MAD"),
  value = c(EXPR_FILE,SURV_FILE,PHENO_FILE,NETWORK_TYPE,TOM_TYPE,COR_TYPE,R2_TARGET,MIN_MODULE_SIZE,MERGE_CUT_HEIGHT,MAX_GENES_MAD)
)
safe_fwrite(param_tbl, file.path(LOGDIR, "run_parameters.tsv"), sep = "\t")

expr_pack <- load_expr_geo(EXPR_FILE, MIN_SAMPLES_EXPR, MAX_GENES_MAD)
mat <- expr_pack$mat

safe_fwrite(data.frame(sample = rownames(mat), mat, check.names = FALSE),
            file.path(TABDIR, "datExpr_used.tsv"), sep = "\t")

qc_tbl <- data.frame(
  metric = c("n_samples", "n_genes", "expression_threshold_used"),
  value = c(nrow(mat), ncol(mat), expr_pack$threshold_used)
)
safe_fwrite(qc_tbl, file.path(TABDIR, "qc_summary.tsv"), sep = "\t")

sampleTree <- hclust(dist(mat), method = "average")
save_pdf_png("FIG_sample_clustering_dendrogram", 12, 6, function() {
  par(cex = 0.7)
  plot(sampleTree, main = "Sample clustering (outlier QC)", xlab = "", sub = "")
  if (!is.null(SAMPLE_CUT_HEIGHT)) abline(h = SAMPLE_CUT_HEIGHT, col = "red")
})

if (DO_SAMPLE_OUTLIER_QC && !is.null(SAMPLE_CUT_HEIGHT)) {
  cl <- cutreeStatic(sampleTree, cutHeight = SAMPLE_CUT_HEIGHT, minSize = 10)
  mat <- mat[cl == 1, , drop = FALSE]
}

if (!file.exists(FERR_DRIVER) || !file.exists(FERR_SUPP)) {
  stop("Missing FerrDb files (ferrdb_driver.csv / ferrdb_suppressor.csv).")
}
fd <- fread(FERR_DRIVER, data.table = FALSE)
fs <- fread(FERR_SUPP, data.table = FALSE)
fm <- if (file.exists(FERR_MARKER)) fread(FERR_MARKER, data.table = FALSE) else NULL

dc <- ifelse("Symbol" %in% colnames(fd), "Symbol", colnames(fd)[1])
sc <- ifelse("Symbol" %in% colnames(fs), "Symbol", colnames(fs)[1])
mc <- if (!is.null(fm)) ifelse("Symbol" %in% colnames(fm), "Symbol", colnames(fm)[1]) else NA_character_

driver_genes <- unique(as.character(fd[[dc]]))
suppressor_genes <- unique(as.character(fs[[sc]]))
marker_genes <- if (!is.null(fm)) unique(as.character(fm[[mc]])) else character(0)

genesets <- list(
  FERR_DRIVER = intersect(driver_genes, colnames(mat)),
  FERR_SUPPRESSOR = intersect(suppressor_genes, colnames(mat))
)
if (length(marker_genes)) genesets$FERR_MARKER <- intersect(marker_genes, colnames(mat))

safe_fwrite(data.frame(geneset = names(genesets), n = sapply(genesets, length)),
            file.path(TABDIR, "ferro_gene_set_sizes.tsv"), sep = "\t")

expr_gene_by_sample <- t(mat)
rownames(expr_gene_by_sample) <- colnames(mat)
colnames(expr_gene_by_sample) <- rownames(mat)

ss <- ssgsea_safe(expr_gene_by_sample, genesets)
ss_df <- as.data.frame(t(ss))
ss_df$sample <- rownames(ss_df)

if (!("FERR_DRIVER" %in% colnames(ss_df))) ss_df$FERR_DRIVER <- NA_real_
if (!("FERR_SUPPRESSOR" %in% colnames(ss_df))) ss_df$FERR_SUPPRESSOR <- NA_real_
if (!("FERR_MARKER" %in% colnames(ss_df))) ss_df$FERR_MARKER <- NA_real_

ss_df$FERRO_DRIVER <- ss_df$FERR_DRIVER
ss_df$FERRO_SUPPRESSOR <- ss_df$FERR_SUPPRESSOR
ss_df$FERRO_MARKER <- ss_df$FERR_MARKER
ss_df$FERRO_TOTAL <- rowMeans(ss_df[, c("FERRO_DRIVER", "FERRO_SUPPRESSOR"), drop = FALSE], na.rm = TRUE)
ss_df$FERRO_TRAIT <- ss_df$FERRO_DRIVER - ss_df$FERRO_SUPPRESSOR
safe_fwrite(ss_df, file.path(TABDIR, "ferroptosis_scores.tsv"), sep = "\t")

traits <- data.frame(
  FERRO_DRIVER = ss_df$FERRO_DRIVER[match(rownames(mat), ss_df$sample)],
  FERRO_SUPPRESSOR = ss_df$FERRO_SUPPRESSOR[match(rownames(mat), ss_df$sample)],
  FERRO_TOTAL = ss_df$FERRO_TOTAL[match(rownames(mat), ss_df$sample)],
  FERRO_TRAIT = ss_df$FERRO_TRAIT[match(rownames(mat), ss_df$sample)],
  check.names = FALSE
)
if ("FERRO_MARKER" %in% colnames(ss_df)) {
  traits$FERRO_MARKER <- ss_df$FERRO_MARKER[match(rownames(mat), ss_df$sample)]
}
rownames(traits) <- rownames(mat)

# ---- FIXED pickSoftThreshold call ----
sft <- call_pickSoftThreshold_compat(
  datExpr = mat,
  powerVector = POWER_RANGE,
  networkType = NETWORK_TYPE,
  cor_type = COR_TYPE,
  verbose = 5
)

fit_tbl <- as.data.frame(sft$fitIndices)
safe_fwrite(fit_tbl, file.path(TABDIR, "soft_threshold_fit.tsv"), sep = "\t")

signedR2 <- -sign(fit_tbl[, 3]) * fit_tbl[, 2]
powers <- fit_tbl[, 1]
meanK <- fit_tbl[, 5]
ok <- which(signedR2 >= R2_TARGET)
chosen_power <- if (length(ok)) as.integer(powers[min(ok)]) else as.integer(powers[which.max(signedR2)])
if (is.finite(FORCE_POWER)) chosen_power <- as.integer(FORCE_POWER)

save_pdf_png("FIG_soft_threshold_selection", 11, 5.5, function() {
  par(mfrow = c(1,2))
  plot(powers, signedR2, xlab="Soft Threshold (power)", ylab="Scale-free fit (signed R^2)",
       type="n", main="Scale independence")
  text(powers, signedR2, labels=powers, col="red", cex=0.9)
  abline(h=R2_TARGET, col="blue", lty=2)
  abline(v=chosen_power, col="darkgreen", lty=2)
  
  plot(powers, meanK, xlab="Soft Threshold (power)", ylab="Mean connectivity",
       type="n", main="Mean connectivity")
  text(powers, meanK, labels=powers, col="red", cex=0.9)
  abline(v=chosen_power, col="darkgreen", lty=2)
})

run_blockwise <- function(datExpr, power, networkType, TOMType, corType, minModuleSize, mergeCutHeight, deepSplit) {
  blockwiseModules(
    datExpr,
    power = power,
    networkType = networkType,
    TOMType = TOMType,
    minModuleSize = minModuleSize,
    mergeCutHeight = mergeCutHeight,
    deepSplit = deepSplit,
    numericLabels = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    maxBlockSize = ncol(datExpr),
    reassignThreshold = 0,
    minKMEtoStay = 0,
    verbose = 5,
    corType = corType
  )
}

net <- tryCatch(
  run_blockwise(mat, chosen_power, NETWORK_TYPE, TOM_TYPE, COR_TYPE, MIN_MODULE_SIZE, MERGE_CUT_HEIGHT, DEEP_SPLIT),
  error = function(e) {
    message("[WARN] primary blockwiseModules failed, fallback running...")
    run_blockwise(
      mat,
      power = max(3, min(chosen_power, 8)),
      networkType = "unsigned",
      TOMType = "unsigned",
      corType = "pearson",
      minModuleSize = max(20, min(30, floor(ncol(mat)/20))),
      mergeCutHeight = 0.30,
      deepSplit = 1
    )
  }
)

module_colors <- net$colors
MEs <- orderMEs(net$MEs)

safe_fwrite(data.frame(gene = colnames(mat), module_color = module_colors),
            file.path(TABDIR, "gene_modules.tsv"), sep = "\t")
safe_fwrite(data.frame(sample = rownames(mat), MEs, check.names = FALSE),
            file.path(TABDIR, "module_eigengenes.tsv"), sep = "\t")

if (length(net$dendrograms) >= 1 && length(net$blockGenes) >= 1) {
  save_pdf_png("FIG_gene_dendrogram_and_module_colors", 12, 7, function() {
    plotDendroAndColors(
      net$dendrograms[[1]],
      module_colors[net$blockGenes[[1]]],
      "Module colors",
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05,
      main = "Gene dendrogram and module colors"
    )
  })
} else {
  save_placeholder_figure("FIG_gene_dendrogram_and_module_colors", "No dendrogram available from WGCNA run.")
}

# Module-trait
ME_mat <- as.data.frame(MEs)
mts <- data.frame()
for (me in colnames(ME_mat)) {
  for (tn in colnames(traits)) {
    sp <- safe_cor_test_spearman(ME_mat[[me]], traits[[tn]])
    pr <- safe_cor_test_pearson(ME_mat[[me]], traits[[tn]])
    mts <- rbind(mts, data.frame(
      module = me, trait = tn,
      rho_spearman = sp$rho, p_spearman = sp$p,
      r_pearson = pr$r, p_pearson = pr$p,
      n = min(sp$n, pr$n)
    ))
  }
}
mts$FDR_spearman <- p_floor(p.adjust(mts$p_spearman, method = "BH"))
mts$FDR_pearson  <- p_floor(p.adjust(mts$p_pearson, method = "BH"))
safe_fwrite(mts, file.path(TABDIR, "module_trait_stats.tsv"), sep = "\t")

mods <- unique(mts$module); trs <- unique(mts$trait)
mat_r <- matrix(NA_real_, nrow = length(mods), ncol = length(trs), dimnames = list(mods, trs))
mat_fdr <- mat_r
for (i in seq_len(nrow(mts))) {
  rr <- mts[i, ]
  mat_r[rr$module, rr$trait] <- rr$r_pearson
  mat_fdr[rr$module, rr$trait] <- rr$FDR_pearson
}
txt <- matrix("", nrow = nrow(mat_r), ncol = ncol(mat_r), dimnames = dimnames(mat_r))
for (i in seq_len(nrow(txt))) for (j in seq_len(ncol(txt))) {
  txt[i, j] <- sprintf("%.2f\nFDR=%.2g", mat_r[i, j], mat_fdr[i, j])
}
save_pdf_png("FIG_module_trait_heatmap", 9, 6.5, function() {
  labeledHeatmap(
    Matrix = mat_r, xLabels = colnames(mat_r), yLabels = rownames(mat_r), ySymbols = rownames(mat_r),
    colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt,
    setStdMargins = FALSE, cex.text = 0.68, zlim = c(-1, 1),
    main = "Module-trait relationships (Pearson)"
  )
})

# ── v13 FIX-2 + FIX-5: FerrDb-aware modül seçimi ─────────────────────────────
# Tüm FerrDb genleri (driver ∪ suppressor ∪ marker)
ferrdb_all <- unique(c(driver_genes, suppressor_genes, marker_genes))

pick <- mts[mts$trait == MODULE_PICK_TRAIT & is.finite(mts$r_pearson), , drop = FALSE]
pick_candidates <- pick[pick$r_pearson >= MODULE_RHO_MIN & pick$FDR_pearson <= MODULE_FDR_MAX, , drop = FALSE]
if (nrow(pick_candidates) == 0) {
  # Fallback: FDR şartını kaldır, sadece r kriterine bak
  pick_candidates <- pick[pick$r_pearson > 0, , drop = FALSE]
  if (nrow(pick_candidates) == 0) stop("No positive module for MODULE_PICK_TRAIT.")
  message("[WARN] No module met r>=", MODULE_RHO_MIN, " AND FDR<=", MODULE_FDR_MAX,
          "; relaxed to r>0 fallback.")
}
pick_candidates <- pick_candidates[order(-pick_candidates$r_pearson, pick_candidates$FDR_pearson), ]

# v13 FIX-5: MIN_FERRO_IN_MODULE garantisi
# Adaylar içinden yeterli FerrDb geni olan modülü seç
best <- NULL
for (i in seq_len(nrow(pick_candidates))) {
  cand_mod   <- as.character(pick_candidates$module[i])
  cand_color <- sub("^ME", "", cand_mod)
  cand_genes <- unique(colnames(mat)[module_colors == cand_color])
  n_ferro    <- length(intersect(cand_genes, ferrdb_all))
  message(sprintf("  Module candidate: %s | r=%.3f | FDR=%.3g | FerrDb genes in module: %d",
                  cand_color,
                  pick_candidates$r_pearson[i],
                  pick_candidates$FDR_pearson[i],
                  n_ferro))
  if (n_ferro >= MIN_FERRO_IN_MODULE) {
    best <- pick_candidates[i, , drop = FALSE]
    message(sprintf("  → Selected: %s (first module with >= %d FerrDb genes)",
                    cand_color, MIN_FERRO_IN_MODULE))
    break
  }
}
if (is.null(best)) {
  # Fallback: FerrDb garantisi olmadan en yüksek r'ye sahip modülü seç
  best <- pick_candidates[1, , drop = FALSE]
  message("[WARN] No module had >= ", MIN_FERRO_IN_MODULE,
          " FerrDb genes. Falling back to highest-r module: ",
          sub("^ME", "", best$module))
}

best_module <- as.character(best$module)
best_color  <- sub("^ME", "", best_module)
sel_genes   <- unique(colnames(mat)[module_colors == best_color])

message(sprintf("[MODULE] Selected: %s | n_genes=%d | r_pearson=%.3f | FDR=%.3g",
                best_color, length(sel_genes),
                best$r_pearson, best$FDR_pearson))

# Log modül seçim gerekçesi
module_select_log <- data.frame(
  selected_module        = best_color,
  r_pearson_FERRO_TRAIT  = round(best$r_pearson, 4),
  FDR_pearson            = signif(best$FDR_pearson, 3),
  n_total_module_genes   = length(sel_genes),
  n_ferrdb_in_module     = length(intersect(sel_genes, ferrdb_all)),
  n_ferrdb_driver        = length(intersect(sel_genes, driver_genes)),
  n_ferrdb_suppressor    = length(intersect(sel_genes, suppressor_genes)),
  min_ferro_threshold    = MIN_FERRO_IN_MODULE,
  stringsAsFactors       = FALSE
)
safe_fwrite(module_select_log, file.path(LOGDIR, "module_selection_log.tsv"), sep = "\t")

safe_fwrite(data.frame(sample = rownames(mat), ME_selected = ME_mat[[best_module]]),
            file.path(TABDIR, "wgcna_module_scores.csv"), sep = "\t")

# GS/MM
me_vec <- ME_mat[[best_module]]
trait_vec <- traits[[MODULE_PICK_TRAIT]]
gsmm <- data.frame(gene = sel_genes, MM = NA_real_, MM_p = 1, GS = NA_real_, GS_p = 1)
for (i in seq_len(nrow(gsmm))) {
  g <- gsmm$gene[i]
  x <- mat[, g]
  a <- safe_cor_test_spearman(x, me_vec)
  b <- safe_cor_test_spearman(x, trait_vec)
  gsmm$MM[i] <- a$rho; gsmm$MM_p[i] <- a$p
  gsmm$GS[i] <- b$rho; gsmm$GS_p[i] <- b$p
}
gsmm$MM_FDR <- p_floor(p.adjust(gsmm$MM_p, "BH"))
gsmm$GS_FDR <- p_floor(p.adjust(gsmm$GS_p, "BH"))
gsmm$absMM <- abs(gsmm$MM); gsmm$absGS <- abs(gsmm$GS)
safe_fwrite(gsmm, file.path(TABDIR, "module_genes_GS_MM.tsv"), sep = "\t")

hub <- gsmm[
  is.finite(gsmm$absMM) & is.finite(gsmm$absGS) &
    gsmm$absMM >= ABS_MM_MIN & gsmm$absGS >= ABS_GS_MIN &
    gsmm$MM_p <= MM_P_MAX & gsmm$GS_FDR <= GS_FDR_MAX, , drop = FALSE
]
hub <- hub[order(-hub$absMM, -hub$absGS, hub$GS_FDR), , drop = FALSE]
safe_fwrite(hub, file.path(TABDIR, "module_hub_genes_MM_GS.tsv"), sep = "\t")

# ══════════════════════════════════════════════════════════════════════════════
# v13 FIX-1: DUAL-BRANCH İMPLEMENTASYONU
# conservative: sel_genes ∩ ferrdb_all  (biyolojik olarak doğrulanmış)
# discovery:    sel_genes \ ferrdb_all  (yeni aday genler, yüksek MM+GS)
# ══════════════════════════════════════════════════════════════════════════════

# ── Conservative branch ────────────────────────────────────────────────────────
conservative_genes <- intersect(sel_genes, ferrdb_all)
conservative_df <- gsmm[gsmm$gene %in% conservative_genes, , drop = FALSE]
conservative_df$ferrdb_role <- dplyr::case_when(
  conservative_df$gene %in% driver_genes     ~ "driver",
  conservative_df$gene %in% suppressor_genes ~ "suppressor",
  conservative_df$gene %in% marker_genes     ~ "marker",
  TRUE                                        ~ "other"
) %||% ifelse(conservative_df$gene %in% driver_genes, "driver",
              ifelse(conservative_df$gene %in% suppressor_genes, "suppressor",
                     ifelse(conservative_df$gene %in% marker_genes, "marker", "other")))
# Base R uyumu için %||% yerine:
conservative_df$ferrdb_role <- ifelse(
  conservative_df$gene %in% driver_genes,     "driver",
  ifelse(conservative_df$gene %in% suppressor_genes, "suppressor",
         ifelse(conservative_df$gene %in% marker_genes, "marker", "other")))
conservative_df <- conservative_df[order(-conservative_df$absMM, -conservative_df$absGS), , drop = FALSE]

message(sprintf("[CONSERVATIVE] %d FerrDb genes in module (driver=%d, suppressor=%d, marker=%d)",
                length(conservative_genes),
                sum(conservative_df$ferrdb_role == "driver",   na.rm = TRUE),
                sum(conservative_df$ferrdb_role == "suppressor", na.rm = TRUE),
                sum(conservative_df$ferrdb_role == "marker",   na.rm = TRUE)))

# ── v13 YENİ: Büyük modül uyarısı ve hub-score sıralaması ─────────────────────
if (length(sel_genes) > MAX_MODULE_GENES_ML) {
  message(sprintf(
    "[WARN] Modül çok büyük: %d gen (MAX_MODULE_GENES_ML=%d). \n" ,
    length(sel_genes), MAX_MODULE_GENES_ML,
    "  Conservative branch hub-score (MM×|GS|) sıralamasıyla daraltılıyor."))

  # Hub score: MM × |GS| — yüksek hem module membership hem ferroptosis ilişkisi
  conservative_df$hub_score <- conservative_df$absMM * conservative_df$absGS
  conservative_df <- conservative_df[order(-conservative_df$hub_score), , drop = FALSE]

  # Top HUB_SCORE_TOPN ile daralt
  if (nrow(conservative_df) > HUB_SCORE_TOPN) {
    conservative_df_full <- conservative_df   # tam listeyi sakla
    conservative_df      <- head(conservative_df, HUB_SCORE_TOPN)
    conservative_genes   <- conservative_df$gene
    message(sprintf("  [HUB-SCORE] Conservative branch %d → %d gen (hub_score top%d)",
                    nrow(conservative_df_full), length(conservative_genes), HUB_SCORE_TOPN))
    safe_fwrite(conservative_df_full,
                file.path(TABDIR, "conservative_genes_full_unfiltered.tsv"), sep = "\t")
  }
} else {
  conservative_df$hub_score <- conservative_df$absMM * conservative_df$absGS
  conservative_df <- conservative_df[order(-conservative_df$hub_score), , drop = FALSE]
}

# ── Discovery branch ──────────────────────────────────────────────────────────
novel_candidates_all <- setdiff(sel_genes, ferrdb_all)
discovery_df <- gsmm[gsmm$gene %in% novel_candidates_all, , drop = FALSE]

# Discovery filtre: yüksek MM + yüksek GS (FERRO_TRAIT)
discovery_filtered <- discovery_df[
  is.finite(discovery_df$absMM) & is.finite(discovery_df$absGS) &
    discovery_df$absMM >= MM_MIN_ABS &
    discovery_df$absGS >= GS_MIN_ABS, , drop = FALSE
]
discovery_filtered <- discovery_filtered[
  order(-discovery_filtered$absMM, -discovery_filtered$absGS), , drop = FALSE
]

# Univariate Cox filtresi (survival verisi varsa)
cox_surv_df_tmp <- load_survival_geo(SURV_FILE, rownames(mat))
if (!is.null(cox_surv_df_tmp) && nrow(discovery_filtered) > 0 &&
    requireNamespace("survival", quietly = TRUE)) {
  common_tmp <- intersect(rownames(mat), rownames(cox_surv_df_tmp))
  if (length(common_tmp) >= MIN_SAMPLES_COX) {
    cox_p_discovery <- sapply(discovery_filtered$gene, function(g) {
      if (!g %in% colnames(mat)) return(1)
      x   <- mat[common_tmp, g]
      ss  <- cox_surv_df_tmp[common_tmp, , drop = FALSE]
      ok  <- is.finite(x) & is.finite(ss$time) & is.finite(ss$event)
      if (sum(ok) < MIN_SAMPLES_COX) return(1)
      fu <- tryCatch(
        survival::coxph(survival::Surv(time, event) ~ x,
                        data = data.frame(time  = ss$time[ok],
                                          event = ss$event[ok],
                                          x     = x[ok])),
        error = function(e) NULL)
      if (is.null(fu)) return(1)
      as.numeric(summary(fu)$coefficients[1, "Pr(>|z|)"])
    })
    discovery_filtered$cox_p_univariate <- cox_p_discovery
    discovery_filtered$cox_p_adj        <- p.adjust(cox_p_discovery, "BH")
    # Sadece Cox p < UNIV_P_CUTOFF olanları tut
    discovery_filtered <- discovery_filtered[
      is.finite(discovery_filtered$cox_p_univariate) &
        discovery_filtered$cox_p_univariate < UNIV_P_CUTOFF, , drop = FALSE
    ]
    message(sprintf("[DISCOVERY] After Cox filter (p<%.2f): %d novel candidates",
                    UNIV_P_CUTOFF, nrow(discovery_filtered)))
  }
} else {
  discovery_filtered$cox_p_univariate <- NA_real_
  discovery_filtered$cox_p_adj        <- NA_real_
}

# Top N discovery (EXPORT_MODULE_TOPN)
discovery_genes <- head(discovery_filtered$gene, EXPORT_MODULE_TOPN)
message(sprintf("[DISCOVERY] Final novel candidates: %d genes", length(discovery_genes)))

# ── Çıktılar ────────────────────────────────────────────────────────────────
if (ANALYSIS_MODE %in% c("conservative", "both")) {
  # ferroptosis_related_genes.csv = conservative (biyolojik filtreli)
  safe_fwrite(data.frame(gene = conservative_genes),
              file.path(TABDIR, "ferroptosis_related_genes.csv"), sep = ",")
  safe_fwrite(conservative_df,
              file.path(TABDIR, "conservative_ferroptosis_genes_GS_MM.tsv"), sep = "\t")
  message("[CONSERVATIVE] Saved: ferroptosis_related_genes.csv (",
          length(conservative_genes), " genes)")
}

if (ANALYSIS_MODE %in% c("discovery", "both")) {
  safe_fwrite(data.frame(gene = discovery_genes),
              file.path(TABDIR, "discovery_candidates.csv"), sep = ",")
  safe_fwrite(discovery_filtered,
              file.path(TABDIR, "discovery_candidates_full_stats.tsv"), sep = "\t")
  message("[DISCOVERY] Saved: discovery_candidates.csv (",
          length(discovery_genes), " genes)")
}

# ML için birleşik gen listesi: conservative ∪ discovery
ml_gene_union <- unique(c(conservative_genes, discovery_genes))
safe_fwrite(data.frame(
  gene   = ml_gene_union,
  branch = c(rep("conservative", length(conservative_genes)),
             rep("discovery",    length(discovery_genes)))
), file.path(TABDIR, "wgcna_genes_for_ML.csv"), sep = ",")

# Dual-branch özet tablo
dual_summary <- data.frame(
  branch                    = c("conservative", "discovery", "union"),
  n_genes                   = c(length(conservative_genes),
                                 length(discovery_genes),
                                 length(ml_gene_union)),
  description               = c(
    paste0("FerrDb genes in module (driver/suppressor/marker intersect)"),
    paste0("Novel module genes: MM>=", MM_MIN_ABS,
           ", GS>=", GS_MIN_ABS,
           ", Cox p<", UNIV_P_CUTOFF),
    "Union for ML pipeline"
  ),
  stringsAsFactors = FALSE
)
safe_fwrite(dual_summary, file.path(TABDIR, "dual_branch_summary.tsv"), sep = "\t")
message("[DUAL-BRANCH] Summary:")
print(dual_summary)

# Eski sel_genes tabanlı tüm modül genleri (referans için koru)
safe_fwrite(data.frame(gene = sel_genes),
            file.path(TABDIR, "all_module_genes_unfiltered.csv"), sep = ",")

# Conservative + Discovery scatter plot
if (length(conservative_genes) >= 2 || length(discovery_genes) >= 2) {
  save_pdf_png("FIG_dual_branch_GS_MM", 8, 6.5, function() {
    xlim_r <- range(gsmm$MM, na.rm = TRUE)
    ylim_r <- range(gsmm$GS, na.rm = TRUE)
    plot(gsmm$MM, gsmm$GS, pch = 16, cex = 0.55, col = "grey70",
         xlim = xlim_r, ylim = ylim_r,
         xlab = paste0("Module Membership (", best_module, ")"),
         ylab = paste0("Gene Significance (", MODULE_PICK_TRAIT, ")"),
         main = "Dual-Branch Gene Classification")
    abline(v = c(-MM_MIN_ABS, MM_MIN_ABS), lty = 2, col = "grey40")
    abline(h = c(-GS_MIN_ABS, GS_MIN_ABS), lty = 2, col = "grey40")
    # Conservative (FerrDb)
    cf <- gsmm[gsmm$gene %in% conservative_genes, , drop = FALSE]
    if (nrow(cf)) points(cf$MM, cf$GS, pch = 17, cex = 0.9, col = "#D85A30")
    # Discovery (novel)
    df2 <- gsmm[gsmm$gene %in% discovery_genes, , drop = FALSE]
    if (nrow(df2)) points(df2$MM, df2$GS, pch = 15, cex = 0.8, col = "#378ADD")
    legend("topright",
           legend = c(paste0("Conservative (FerrDb, n=", length(conservative_genes), ")"),
                      paste0("Discovery (novel, n=", length(discovery_genes), ")"),
                      "Other module genes"),
           col = c("#D85A30", "#378ADD", "grey70"),
           pch = c(17, 15, 16), pt.cex = c(0.9, 0.8, 0.55), cex = 0.8)
  })
}

save_pdf_png("FIG_GS_vs_MM_selected_module", 6.5, 5.5, function() {
  plot(gsmm$MM, gsmm$GS, pch = 16, cex = 0.7,
       xlab = paste0("MM (cor gene, ", best_module, ")"),
       ylab = paste0("GS (cor gene, ", MODULE_PICK_TRAIT, ")"),
       main = paste("Selected module:", best_color))
  abline(v = c(-ABS_MM_MIN, ABS_MM_MIN), lty = 2, col = "grey40")
  abline(h = c(-ABS_GS_MIN, ABS_GS_MIN), lty = 2, col = "grey40")
  if (nrow(hub) > 0) points(hub$MM, hub$GS, pch = 16, cex = 0.8, col = "red")
})

# Survival/Cox
surv_df <- load_survival_geo(SURV_FILE, rownames(mat))
stage_df <- load_stage_geo(PHENO_FILE, rownames(mat))

cox_uni <- data.frame(); cox_adj <- data.frame()

if (!is.null(surv_df) && requireNamespace("survival", quietly = TRUE)) {
  common <- intersect(rownames(mat), rownames(surv_df))
  Xs <- mat[common, sel_genes, drop = FALSE]
  Ss <- surv_df[common, , drop = FALSE]
  st <- if (!is.null(stage_df)) stage_df$stage_num[match(common, stage_df$sample)] else rep(NA_real_, length(common))
  
  for (g in colnames(Xs)) {
    x <- Xs[, g]
    ok <- is.finite(x) & is.finite(Ss$time) & is.finite(Ss$event)
    if (sum(ok) < MIN_SAMPLES_COX) next
    
    fu <- tryCatch(
      survival::coxph(survival::Surv(time, event) ~ x,
                      data = data.frame(time = Ss$time[ok], event = Ss$event[ok], x = x[ok])),
      error = function(e) NULL
    )
    if (!is.null(fu)) {
      su <- summary(fu)
      cox_uni <- rbind(cox_uni, data.frame(
        gene = g, beta = su$coefficients[1, "coef"], HR = su$coefficients[1, "exp(coef)"],
        z = su$coefficients[1, "z"], p = p_floor(su$coefficients[1, "Pr(>|z|)"]), n = sum(ok)
      ))
    }
    
    if (sum(is.finite(st)) >= MIN_SAMPLES_COX) {
      ok2 <- ok & is.finite(st)
      if (sum(ok2) >= MIN_SAMPLES_COX) {
        fa <- tryCatch(
          survival::coxph(survival::Surv(time, event) ~ x + stage_num,
                          data = data.frame(time = Ss$time[ok2], event = Ss$event[ok2], x = x[ok2], stage_num = st[ok2])),
          error = function(e) NULL
        )
        if (!is.null(fa)) {
          sa <- summary(fa)
          if ("x" %in% rownames(sa$coefficients)) {
            cox_adj <- rbind(cox_adj, data.frame(
              gene = g, beta = sa$coefficients["x", "coef"], HR = sa$coefficients["x", "exp(coef)"],
              z = sa$coefficients["x", "z"], p = p_floor(sa$coefficients["x", "Pr(>|z|)"]), n = sum(ok2)
            ))
          }
        }
      }
    }
  }
}

if (nrow(cox_uni) > 0) {
  cox_uni$FDR <- p_floor(p.adjust(cox_uni$p, "BH"))
  cox_uni <- cox_uni[order(cox_uni$FDR, cox_uni$p), , drop = FALSE]
}
if (nrow(cox_adj) > 0) {
  cox_adj$FDR <- p_floor(p.adjust(cox_adj$p, "BH"))
  cox_adj <- cox_adj[order(cox_adj$FDR, cox_adj$p), , drop = FALSE]
}

safe_fwrite(cox_uni, file.path(TABDIR, "univariate_cox_module_genes.tsv"), sep = "\t")
safe_fwrite(cox_adj, file.path(TABDIR, "stage_adjusted_cox_module_genes.tsv"), sep = "\t")
safe_fwrite(if (nrow(cox_uni)) cox_uni[cox_uni$FDR < 0.05, , drop = FALSE] else data.frame(),
            file.path(TABDIR, "prognostic_genes_univariate_FDR0.05.csv"), sep = ",")
safe_fwrite(if (nrow(cox_adj)) cox_adj[cox_adj$FDR < 0.05, , drop = FALSE] else data.frame(),
            file.path(TABDIR, "prognostic_genes_stage_adjusted_FDR0.05.csv"), sep = ",")
safe_fwrite(if (nrow(cox_uni)) cox_uni[seq_len(min(TOP_GENES_FOR_ML, nrow(cox_uni))), , drop = FALSE] else data.frame(),
            file.path(TABDIR, "top_genes_for_ML_univariate.tsv"), sep = "\t")
safe_fwrite(if (nrow(cox_adj)) cox_adj[seq_len(min(TOP_GENES_FOR_ML, nrow(cox_adj))), , drop = FALSE] else data.frame(),
            file.path(TABDIR, "top_genes_for_ML_stage_adjusted.tsv"), sep = "\t")

# v13: possible_prognostic_genes_consensus = conservative ∪ discovery ∪ hub ∪ Cox-sig
possible <- sort(unique(c(
  conservative_genes,                                                            # FerrDb ∩ module
  discovery_genes,                                                               # novel candidates
  if (nrow(hub)) hub$gene else character(0),                                    # yüksek MM+GS
  if (nrow(cox_uni)) cox_uni$gene[cox_uni$FDR < 0.05] else character(0),       # Cox univariate sig
  if (nrow(cox_adj)) cox_adj$gene[cox_adj$FDR < 0.05] else character(0)        # Cox stage-adj sig
)))
safe_fwrite(data.frame(gene = possible),
            file.path(TABDIR, "possible_prognostic_genes_consensus.csv"), sep = ",")
message(sprintf("[CONSENSUS] %d prognostic candidate genes (conservative=%d, discovery=%d, hub=%d)",
                length(possible), length(conservative_genes), length(discovery_genes), nrow(hub)))

if (nrow(cox_uni) > 0) {
  save_pdf_png("FIG_univariate_cox_volcano", 7, 6, function() {
    x <- log2(cox_uni$HR); y <- -log10(p_floor(cox_uni$FDR))
    plot(x, y, pch = 16, cex = 0.6, xlab = "log2(HR)", ylab = "-log10(FDR)",
         main = "Univariate Cox (selected module genes)")
    abline(h = -log10(0.05), lty = 2, col = "grey40")
  })
} else {
  save_placeholder_figure("FIG_univariate_cox_volcano", "No univariate Cox result.")
}

if (nrow(cox_adj) > 0) {
  save_pdf_png("FIG_stage_adjusted_cox_volcano", 7, 6, function() {
    x <- log2(cox_adj$HR); y <- -log10(p_floor(cox_adj$FDR))
    plot(x, y, pch = 16, cex = 0.6, xlab = "log2(HR)", ylab = "-log10(FDR)",
         main = "Stage-adjusted Cox (selected module genes)")
    abline(h = -log10(0.05), lty = 2, col = "grey40")
  })
} else {
  save_placeholder_figure("FIG_stage_adjusted_cox_volcano", "No stage-adjusted Cox result.")
}

# ---- FIXED adjacency ----
if (length(sel_genes) >= 5) {
  Xmod <- mat[, sel_genes, drop = FALSE]
  adj <- tryCatch(call_adjacency_compat(Xmod, power = chosen_power, networkType = NETWORK_TYPE, cor_type = COR_TYPE),
                  error = function(e) NULL)
  
  if (!is.null(adj)) {
    thr <- as.numeric(quantile(adj[upper.tri(adj)], probs = CYTO_WEIGHT_Q, na.rm = TRUE))
    if (!is.finite(thr)) thr <- 0.1
    idx <- which(adj >= thr & upper.tri(adj), arr.ind = TRUE)
    if (nrow(idx) > CYTO_MAX_EDGES) {
      ord <- order(adj[idx], decreasing = TRUE)
      idx <- idx[ord[seq_len(CYTO_MAX_EDGES)], , drop = FALSE]
    }
    
    edges <- data.frame(from = colnames(adj)[idx[, 1]], to = colnames(adj)[idx[, 2]], weight = as.numeric(adj[idx]))
    node <- data.frame(gene = sel_genes, module = best_color, stringsAsFactors = FALSE)
    node <- merge(node, gsmm[, c("gene","MM","GS","absMM","absGS","GS_FDR","MM_FDR"), drop = FALSE],
                  by = "gene", all.x = TRUE, sort = FALSE)
    
    safe_fwrite(edges, file.path(CYTODIR, "cytoscape_edges_selected_module.tsv"), sep = "\t")
    safe_fwrite(node, file.path(CYTODIR, "cytoscape_nodes_selected_module.tsv"), sep = "\t")
  } else {
    safe_fwrite(data.frame(), file.path(CYTODIR, "cytoscape_edges_selected_module.tsv"), sep = "\t")
    safe_fwrite(data.frame(), file.path(CYTODIR, "cytoscape_nodes_selected_module.tsv"), sep = "\t")
  }
} else {
  safe_fwrite(data.frame(), file.path(CYTODIR, "cytoscape_edges_selected_module.tsv"), sep = "\t")
  safe_fwrite(data.frame(), file.path(CYTODIR, "cytoscape_nodes_selected_module.tsv"), sep = "\t")
}

mod_tab <- sort(table(module_colors), decreasing = TRUE)
safe_fwrite(data.frame(module = names(mod_tab), n_genes = as.integer(mod_tab)),
            file.path(TABDIR, "module_sizes.tsv"), sep = "\t")
save_pdf_png("FIG_module_size_barplot", 9, 5, function() {
  barplot(mod_tab, las = 2, cex.names = 0.8, main = "Module sizes", ylab = "Number of genes")
})

summary_tbl <- data.frame(
  script_version           = "v13",
  expr_used                = basename(EXPR_FILE),
  n_samples                = nrow(mat),
  n_genes                  = ncol(mat),
  chosen_power             = chosen_power,
  r2_target                = R2_TARGET,
  network_type             = NETWORK_TYPE,
  tom_type                 = TOM_TYPE,
  cor_type                 = COR_TYPE,              # v13: bicor
  module_pick_trait        = MODULE_PICK_TRAIT,
  selected_module          = best_module,
  selected_color           = best_color,
  n_selected_module_genes  = length(sel_genes),
  n_ferrdb_in_module       = length(intersect(sel_genes, ferrdb_all)),
  n_conservative_genes     = length(conservative_genes),
  n_conservative_hub_filtered = nrow(conservative_df),
  n_discovery_genes        = length(discovery_genes),
  n_ml_union_genes         = length(ml_gene_union),
  module_too_large         = length(sel_genes) > MAX_MODULE_GENES_ML,
  max_module_genes_ml      = MAX_MODULE_GENES_ML,
  hub_score_topn           = HUB_SCORE_TOPN,
  n_hub_genes              = nrow(hub),
  analysis_mode            = ANALYSIS_MODE,
  min_ferro_in_module      = MIN_FERRO_IN_MODULE,
  survival_used            = !is.null(surv_df),
  stage_used               = !is.null(stage_df),
  stringsAsFactors         = FALSE
)
safe_fwrite(summary_tbl, file.path(OUTDIR, "wgcna_run_summary.tsv"), sep = "\t")

sink(file.path(LOGDIR, "sessionInfo.txt")); print(sessionInfo()); sink()

cat("\n=====================================================\n")
cat("[DONE] WGCNA v13 pipeline completed.\n")
cat("[OUTDIR]  ", OUTDIR, "\n", sep = "")
cat("[CONSERVATIVE] ferroptosis_related_genes.csv\n")
cat("[DISCOVERY]    discovery_candidates.csv\n")
cat("[ML UNION]     wgcna_genes_for_ML.csv\n")
cat("[SUMMARY]      dual_branch_summary.tsv\n")
cat("=====================================================\n")
