#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE, warn = 1)
})

# ============================================================
# Rana Salihoglu 
# ============================================================

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
  "07_signature_reversal_connectivity_q1_package_driven"
}

project_root <- path.expand(arg_value("--project_root", getwd()))
gene_file <- path.expand(arg_value("--gene_file", "possible_prognostic_genes_FULLTRAIN.csv"))
path_manifest <- path.expand(arg_value("--path_manifest", "full_paths(2).txt"))
outdir_name <- arg_value("--outdir_name", script_basename)
outdir_name <- gsub("[^A-Za-z0-9._-]+", "_", outdir_name)
if (!nzchar(outdir_name)) outdir_name <- script_basename
install_missing <- tolower(arg_value("--install_missing", "false")) %in% c("true", "1", "yes", "y")
min_gs_size <- suppressWarnings(as.integer(arg_value("--min_gs_size", "10")))
max_gs_size <- suppressWarnings(as.integer(arg_value("--max_gs_size", "500")))
max_reverse_terms <- suppressWarnings(as.integer(arg_value("--max_reverse_terms", "25")))
max_hallmark_terms <- suppressWarnings(as.integer(arg_value("--max_hallmark_terms", "20")))
max_canonical_terms <- suppressWarnings(as.integer(arg_value("--max_canonical_terms", "20")))
max_overlap_terms <- suppressWarnings(as.integer(arg_value("--max_overlap_terms", "12")))
shortlist_n <- suppressWarnings(as.integer(arg_value("--shortlist_n", "12")))
top_leading_edge <- suppressWarnings(as.integer(arg_value("--top_leading_edge", "10")))
druglike_only_shortlist <- tolower(arg_value("--druglike_only_shortlist", "false")) %in% c("true", "1", "yes", "y")
use_reactome_enrich <- tolower(arg_value("--use_reactome_enrich", "true")) %in% c("true", "1", "yes", "y")

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

read_manifest <- function(path_manifest) {
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

manifest_paths <- read_manifest(path_manifest)

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
    looks_like_root <- file.exists(file.path(cand_n, "OUT_v44")) ||
      file.exists(file.path(cand_n, "GEO_PREP")) ||
      file.exists(file.path(cand_n, "ferrdb_driver.csv"))
    if (looks_like_root) return(cand_n)
  }
  if (length(manifest_paths) > 0) {
    root_hits <- unique(dirname(dirname(manifest_paths[grepl("/(OUT_v44|GEO_PREP|bioinf_scripts)/", manifest_paths)])))
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
  file.path(project_root, "bioinf_scripts", "scripts", basename(gene_file)),
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

cran_pkgs <- c(
  "data.table", "dplyr", "tibble", "tidyr", "stringr", "purrr", "readr",
  "ggplot2", "forcats", "scales", "jsonlite", "glue", "cli", "digest"
)
bioc_pkgs <- c(
  "msigdbr", "clusterProfiler", "enrichplot", "org.Hs.eg.db", "AnnotationDbi",
  "ReactomePA", "reactome.db"
)
optional_pkgs <- c("openxlsx", "patchwork", "KEGGREST")

ensure_pkg <- function(pkg, bioc = FALSE, install_missing = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  if (!install_missing) return(FALSE)
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (bioc || pkg %in% bioc_pkgs) {
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
library(jsonlite)
library(glue)
library(cli)
library(digest)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ReactomePA)
library(reactome.db)

has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
has_keggrest <- requireNamespace("KEGGREST", quietly = TRUE)

anno_select <- function(x, keys, columns, keytype) {
  suppressMessages(suppressWarnings(
    AnnotationDbi::select(x, keys = keys, columns = columns, keytype = keytype)
  ))
}

root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
out_dir <- file.path(root, outdir_name)
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
cache_dir <- file.path(out_dir, "cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

read_flex <- function(path, ...) {
  if (is.null(path) || is.na(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) fread(path, sep = "\t", data.table = FALSE, ...)
  else fread(path, data.table = FALSE, ...)
}

safe_write <- function(df, path) {
  if (is.null(df)) return(invisible(NULL))
  readr::write_tsv(as_tibble(df), path)
}

plot_save <- function(p, file_stub, width = 10, height = 7) {
  ggplot2::ggsave(paste0(file_stub, ".png"), plot = p, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(paste0(file_stub, ".pdf"), plot = p, width = width, height = height, bg = "white")
}

q1_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, color = "#4D4D4D", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.25, linetype = "dashed"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#EFEFEF", color = "grey60"),
      strip.text = element_text(face = "bold")
    )
}

top_n_df <- function(df, n, arrange_by = NULL, decreasing = TRUE) {
  if (is.null(df) || nrow(df) == 0) return(df)
  out <- df
  if (!is.null(arrange_by) && arrange_by %in% names(out)) {
    ord <- order(out[[arrange_by]], decreasing = decreasing, na.last = TRUE)
    out <- out[ord, , drop = FALSE]
  }
  utils::head(out, min(n, nrow(out)))
}

sanitize_gene <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- sub("^.*__", "", x)
  x <- gsub("\\.", "-", x)
  toupper(x)
}

wrap_label <- function(x, width = 40) stringr::str_wrap(x, width = width)

safe_rescale <- function(x, to = c(0, 1)) {
  x <- suppressWarnings(as.numeric(x))
  if (length(unique(x[is.finite(x)])) <= 1) return(rep(mean(to), length(x)))
  scales::rescale(x, to = to, from = range(x, na.rm = TRUE))
}

read_gene_table <- function(path) {
  df <- read_flex(path)
  if (is.null(df) || nrow(df) == 0) stop("Gene file not found or empty: ", path)
  nms <- names(df)
  low <- tolower(nms)
  pick <- function(patterns, required = FALSE) {
    for (p in patterns) {
      idx <- grep(p, low, perl = TRUE)
      if (length(idx) > 0) return(nms[idx[1]])
    }
    if (required) stop("Missing required column: ", paste(patterns, collapse = ", "))
    NA_character_
  }
  gene_col <- pick(c("^gene$", "symbol", "feature"), required = FALSE)
  if (is.na(gene_col)) gene_col <- nms[1]
  z_col <- pick(c("z_meta", "^z$", "meta.*z"), required = FALSE)
  bag_col <- pick(c("bag_frac", "bag_freq", "freq"), required = FALSE)

  tibble(
    gene = sanitize_gene(df[[gene_col]]),
    z_meta = if (!is.na(z_col)) suppressWarnings(as.numeric(df[[z_col]])) else 0,
    bag_frac = if (!is.na(bag_col)) suppressWarnings(as.numeric(df[[bag_col]])) else 1
  ) %>%
    filter(!is.na(gene), gene != "") %>%
    group_by(gene) %>%
    summarise(
      z_meta = suppressWarnings(max(z_meta, na.rm = TRUE)),
      bag_frac = suppressWarnings(max(bag_frac, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      z_meta = ifelse(is.finite(z_meta), z_meta, 0),
      bag_frac = ifelse(is.finite(bag_frac), bag_frac, 1),
      direction = ifelse(z_meta > 0, "Risk-promoting", "Protective"),
      abs_z = abs(z_meta),
      stability_tier = case_when(
        bag_frac >= 1.00 ~ "Core",
        bag_frac >= 0.75 ~ "Stable",
        TRUE ~ "Variable"
      )
    ) %>%
    arrange(desc(abs_z), desc(bag_frac), gene)
}

gene_tbl <- read_gene_table(gene_file)

id_map <- anno_select(
  org.Hs.eg.db,
  keys = unique(gene_tbl$gene),
  columns = c("ENTREZID", "SYMBOL", "GENENAME", "ALIAS"),
  keytype = "SYMBOL"
) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  filter(!is.na(ENTREZID))

gene_tbl <- gene_tbl %>% left_join(id_map, by = c("gene" = "SYMBOL"))
ranked <- gene_tbl %>% filter(!is.na(ENTREZID)) %>% arrange(desc(z_meta))
rank_vec <- ranked$z_meta
names(rank_vec) <- ranked$ENTREZID
rank_vec <- sort(rank_vec, decreasing = TRUE)

if (length(rank_vec) < min_gs_size) {
  stop("Too few ranked genes after Entrez mapping: ", length(rank_vec))
}

msig_get <- function(collection, subcollection = NULL, species = "Homo sapiens") {
  args <- list(species = species, collection = collection)
  x <- tryCatch(do.call(msigdbr::msigdbr, args), error = function(e) NULL)
  if (is.null(x)) stop("msigdbr query failed for collection: ", collection)
  if (!is.null(subcollection)) {
    subcol_name <- if ("gs_subcollection" %in% names(x)) "gs_subcollection" else if ("gs_subcat" %in% names(x)) "gs_subcat" else NA_character_
    if (!is.na(subcol_name)) x <- x %>% filter(.data[[subcol_name]] %in% subcollection)
  }
  x
}

term_gene_from_msig <- function(msig_df) {
  gene_col <- if ("entrez_gene" %in% names(msig_df)) "entrez_gene" else if ("ncbi_gene" %in% names(msig_df)) "ncbi_gene" else NA_character_
  if (is.na(gene_col)) stop("No Entrez-compatible gene column found in msigdbr result")
  msig_df %>%
    filter(!is.na(.data[[gene_col]])) %>%
    transmute(gs_name = as.character(gs_name), entrez_gene = as.character(.data[[gene_col]])) %>%
    distinct()
}

term_name_from_msig <- function(msig_df) {
  desc_col <- if ("gs_description" %in% names(msig_df)) "gs_description" else "gs_name"
  subcol_name <- if ("gs_subcollection" %in% names(msig_df)) "gs_subcollection" else if ("gs_subcat" %in% names(msig_df)) "gs_subcat" else NA_character_
  msig_df %>%
    transmute(
      gs_name = as.character(gs_name),
      gs_description = as.character(.data[[desc_col]]),
      gs_subcollection = if (!is.na(subcol_name)) as.character(.data[[subcol_name]]) else NA_character_
    ) %>%
    distinct()
}

msig_cgp <- msig_get("C2", subcollection = c("CGP"))
msig_cp  <- msig_get("C2", subcollection = c("CP:REACTOME", "CP:KEGG", "CP:BIOCARTA", "CP:PID", "CP:WIKIPATHWAYS"))
msig_h   <- msig_get("H")

term2gene_cgp <- term_gene_from_msig(msig_cgp)
term2name_cgp <- term_name_from_msig(msig_cgp)
term2gene_cp  <- term_gene_from_msig(msig_cp)
term2name_cp  <- term_name_from_msig(msig_cp)
term2gene_h   <- term_gene_from_msig(msig_h)
term2name_h   <- term_name_from_msig(msig_h)

run_gsea <- function(rank_vec, term2gene, minGSSize = 10, maxGSSize = 500) {
  clusterProfiler::GSEA(
    geneList = rank_vec,
    TERM2GENE = term2gene,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    verbose = FALSE,
    eps = 1e-10,
    seed = TRUE,
    by = "fgsea"
  )
}

run_enricher <- function(genes, term2gene, universe = NULL) {
  if (length(genes) < min_gs_size) return(NULL)
  clusterProfiler::enricher(
    gene = genes,
    TERM2GENE = term2gene,
    universe = universe,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size
  )
}

extract_tbl <- function(res) {
  if (is.null(res)) return(tibble())
  as_tibble(res@result)
}

# -----------------------------
# package-driven lexicons
# -----------------------------

gsub_multi <- function(x, patterns, replacement = " ") {
  out <- x
  for (p in patterns) out <- gsub(p, replacement, out, perl = TRUE)
  out
}

normalize_text <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub_multi(x, c("[\\/\\|,:;()\\[\\]{}+-]", "[_]", "\\s+"), " ")
  trimws(x)
}

split_tokens <- function(x) {
  toks <- unique(unlist(strsplit(normalize_text(x), " ", fixed = TRUE)))
  toks <- toks[!is.na(toks) & nzchar(toks)]
  toks[nchar(toks) > 1]
}

extract_alpha_tokens <- function(x) {
  toks <- split_tokens(x)
  toks[grepl("^[A-Z0-9-]+$", toks)]
}

safe_save_rds <- function(object, path) {
  tryCatch(saveRDS(object, path), error = function(e) NULL)
}

build_gene_dictionary <- function(cache_dir) {
  cache_file <- file.path(cache_dir, "gene_dictionary.rds")
  if (file.exists(cache_file)) return(readRDS(cache_file))

  gene_df <- bind_rows(
    anno_select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "SYMBOL"), columns = c("SYMBOL"), keytype = "SYMBOL") %>%
      transmute(token = toupper(SYMBOL), source = "SYMBOL"),
    anno_select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "ALIAS"), columns = c("ALIAS"), keytype = "ALIAS") %>%
      transmute(token = toupper(ALIAS), source = "ALIAS")
  ) %>%
    filter(!is.na(token), token != "") %>%
    mutate(token = trimws(token)) %>%
    filter(grepl("^[A-Z0-9-]+$", token), nchar(token) >= 2, nchar(token) <= 15) %>%
    distinct(token, .keep_all = TRUE)

  safe_save_rds(gene_df, cache_file)
  gene_df
}

build_pathway_dictionary <- function(msig_cp, msig_h, cache_dir) {
  cache_file <- file.path(cache_dir, "pathway_dictionary.rds")
  if (file.exists(cache_file)) return(readRDS(cache_file))

  cp_terms <- unique(c(msig_cp$gs_name, msig_cp$gs_description, msig_h$gs_name, msig_h$gs_description))
  cp_terms <- cp_terms[!is.na(cp_terms) & nzchar(cp_terms)]

  reactome_terms <- unique(as.character(AnnotationDbi::toTable(reactome.db::reactomePATHID2NAME)$path_name))
  reactome_terms <- reactome_terms[!is.na(reactome_terms) & nzchar(reactome_terms)]

  kegg_terms <- character()
  if (has_keggrest) {
    kegg_terms <- tryCatch({
      kk <- KEGGREST::keggList("pathway", "hsa")
      unname(sub(" - Homo sapiens \\(human\\)$", "", kk))
    }, error = function(e) character())
  }

  raw_terms_tbl <- tibble(raw_term = unique(c(cp_terms, reactome_terms, kegg_terms))) %>%
    mutate(normalized = normalize_text(raw_term)) %>%
    filter(normalized != "")

  pathway_token_vec <- raw_terms_tbl$normalized %>%
    purrr::map(extract_alpha_tokens) %>%
    unlist(use.names = FALSE)

  pathway_df <- tibble(token = pathway_token_vec) %>%
    filter(!is.na(token), token != "", nchar(token) >= 3) %>%
    count(token, sort = TRUE, name = "n_sources") %>%
    mutate(is_pathway_token = n_sources >= 2 | token %in% c("PATHWAY", "SIGNALING", "SIGNALLING", "RESPONSE", "METABOLISM", "CHECKPOINT", "APOPTOSIS", "AUTOPHAGY", "FERROPTOSIS"))

  safe_save_rds(pathway_df, cache_file)
  pathway_df
}

build_drug_dictionary <- function(cache_dir) {
  cache_file <- file.path(cache_dir, "drug_dictionary.rds")
  if (file.exists(cache_file)) return(readRDS(cache_file))

  drug_df <- tibble(token = character(), source = character())
  if (has_keggrest) {
    drug_df <- tryCatch({
      kd <- KEGGREST::keggList("drug")
      nm <- unname(kd)
      nm <- gsub(";.*$", "", nm)
      tibble(raw_name = nm) %>%
        mutate(norm = normalize_text(raw_name)) %>%
        mutate(token = gsub("[^A-Z0-9 ]", " ", norm)) %>%
        separate_rows(token, sep = "\\s+") %>%
        filter(token != "", nchar(token) >= 4, !grepl("^[0-9]+$", token)) %>%
        transmute(token, source = "KEGG_DRUG") %>%
        distinct()
    }, error = function(e) tibble(token = character(), source = character()))
  }

  safe_save_rds(drug_df, cache_file)
  drug_df
}

gene_dict <- build_gene_dictionary(cache_dir)
pathway_dict <- build_pathway_dictionary(msig_cp, msig_h, cache_dir)
drug_dict <- build_drug_dictionary(cache_dir)

safe_write(gene_dict, file.path(tab_dir, "dictionary_gene_tokens.tsv"))
safe_write(pathway_dict, file.path(tab_dir, "dictionary_pathway_tokens.tsv"))
safe_write(drug_dict, file.path(tab_dir, "dictionary_drug_tokens.tsv"))

# -----------------------------
# CGP term parsing + evidence-backed classification
# -----------------------------

extract_perturbagen_phrase <- function(gs_name) {
  x <- toupper(gs_name)
  x <- gsub("_UP$|_DN$|_DOWN$", "", x)
  x <- gsub("_VS_.*$", "", x)
  x <- gsub("^KEGG_|^REACTOME_|^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  trimws(x)
}

extract_contrast <- function(gs_name) {
  x <- toupper(gs_name)
  if (grepl("_UP$", x)) return("UP gene set")
  if (grepl("_DN$|_DOWN$", x)) return("DOWN gene set")
  "Mixed"
}

count_token_hits <- function(tokens, dictionary_tokens) {
  sum(tokens %in% dictionary_tokens)
}

first_hits <- function(tokens, dictionary_tokens, n = 5) {
  hits <- unique(tokens[tokens %in% dictionary_tokens])
  paste(utils::head(hits, n), collapse = "; ")
}

classify_cgp_term <- function(term_label, drug_tokens, gene_tokens, pathway_tokens) {
  tokens <- extract_alpha_tokens(term_label)
  token_string <- paste(tokens, collapse = " ")

  drug_hits <- count_token_hits(tokens, drug_tokens)
  gene_hits <- count_token_hits(tokens, gene_tokens)
  pathway_hits <- count_token_hits(tokens, pathway_tokens)

  has_genetic_action <- grepl("\\b(KNOCKDOWN|KNOCKOUT|OVEREXPRESSION|DEPLETION|SILENCED|DEFICIENT|MUTANT|WILDTYPE|WT|SHRNA|SIRNA|CRISPR|KD|KO)\\b", token_string)
  has_exposure_action <- grepl("\\b(TREATED|TREATMENT|EXPOSURE|STIMULATED|INHIBITED|AFTER|BEFORE|RESISTANT|SENSITIVE)\\b", token_string)
  has_cell_context <- grepl("\\b(CELL|CELLS|LINE|LINES|PATIENT|TUMOR|TUMOUR|XENOGRAFT|ORGANOID|TISSUE|A549|H1975|H1299)\\b", token_string)

  context_family <- case_when(
    drug_hits >= 1 ~ "Compound / drug-response",
    has_genetic_action && gene_hits >= 1 ~ "Genetic / regulatory perturbation",
    pathway_hits >= 2 ~ "Pathway / program perturbation",
    has_cell_context ~ "Disease / lineage context",
    has_exposure_action ~ "Exposure / treatment context",
    TRUE ~ "Other / unresolved"
  )

  actionability_class <- case_when(
    context_family == "Compound / drug-response" ~ "Directly drug-like",
    context_family == "Genetic / regulatory perturbation" ~ "Mechanistically actionable",
    context_family == "Pathway / program perturbation" ~ "Biology-supportive",
    context_family == "Disease / lineage context" ~ "Phenotype-supportive",
    context_family == "Exposure / treatment context" ~ "Contextual only",
    TRUE ~ "Exploratory"
  )

  tibble(
    matched_drug_tokens_n = drug_hits,
    matched_gene_tokens_n = gene_hits,
    matched_pathway_tokens_n = pathway_hits,
    matched_drug_tokens = first_hits(tokens, drug_tokens),
    matched_gene_tokens = first_hits(tokens, gene_tokens),
    matched_pathway_tokens = first_hits(tokens, pathway_tokens),
    has_genetic_action = has_genetic_action,
    has_exposure_action = has_exposure_action,
    has_cell_context = has_cell_context,
    context_family = context_family,
    actionability_class = actionability_class,
    is_druglike = context_family == "Compound / drug-response"
  )
}

infer_reversal <- function(df, drug_dict, gene_dict, pathway_dict) {
  if (nrow(df) == 0) return(df)

  drug_tokens <- unique(drug_dict$token)
  gene_tokens <- unique(gene_dict$token)
  pathway_tokens <- unique(pathway_dict$token[pathway_dict$is_pathway_token])

  parsed <- purrr::map_dfr(df$ID, ~classify_cgp_term(extract_perturbagen_phrase(.x), drug_tokens, gene_tokens, pathway_tokens))
  bind_cols(df, parsed) %>%
    mutate(
      perturbagen = vapply(ID, extract_perturbagen_phrase, character(1)),
      perturbation_set_direction = vapply(ID, extract_contrast, character(1)),
      reversal_class = case_when(
        NES < 0 ~ "Candidate reversal",
        NES > 0 ~ "Concordant / phenocopying",
        TRUE ~ "Indeterminate"
      ),
      absNES = abs(NES),
      neglog10_fdr = -log10(pmax(p.adjust, 1e-300))
    ) %>%
    arrange(p.adjust, desc(absNES))
}

# -----------------------------
# enrichment analyses
# -----------------------------

res_cgp <- tryCatch(run_gsea(rank_vec, term2gene_cgp, minGSSize = min_gs_size, maxGSSize = max_gs_size), error = function(e) NULL)
res_cp  <- tryCatch(run_gsea(rank_vec, term2gene_cp,  minGSSize = min_gs_size, maxGSSize = max_gs_size), error = function(e) NULL)
res_h   <- tryCatch(run_gsea(rank_vec, term2gene_h,   minGSSize = min_gs_size, maxGSSize = max_gs_size), error = function(e) NULL)

tbl_cgp <- extract_tbl(res_cgp) %>% left_join(term2name_cgp, by = c("ID" = "gs_name"))
tbl_cp  <- extract_tbl(res_cp)  %>% left_join(term2name_cp,  by = c("ID" = "gs_name"))
tbl_h   <- extract_tbl(res_h)   %>% left_join(term2name_h,   by = c("ID" = "gs_name"))

tbl_cgp <- infer_reversal(tbl_cgp, drug_dict = drug_dict, gene_dict = gene_dict, pathway_dict = pathway_dict)

tbl_cp <- tbl_cp %>%
  mutate(
    absNES = abs(NES),
    neglog10_fdr = -log10(pmax(p.adjust, 1e-300)),
    canonical_direction = ifelse(NES > 0, "Activated in risk-high signature", "Suppressed in risk-high signature")
  ) %>%
  arrange(p.adjust, desc(absNES))

tbl_h <- tbl_h %>%
  mutate(
    absNES = abs(NES),
    neglog10_fdr = -log10(pmax(p.adjust, 1e-300)),
    hallmark_direction = ifelse(NES > 0, "Enriched in risk-high signature", "Depleted in risk-high signature")
  ) %>%
  arrange(p.adjust, desc(absNES))

risk_entrez <- gene_tbl %>% filter(direction == "Risk-promoting", !is.na(ENTREZID)) %>% pull(ENTREZID) %>% as.character()
protective_entrez <- gene_tbl %>% filter(direction == "Protective", !is.na(ENTREZID)) %>% pull(ENTREZID) %>% as.character()
universe_entrez <- unique(as.character(gene_tbl$ENTREZID[!is.na(gene_tbl$ENTREZID)]))

ora_risk_cgp <- tryCatch(run_enricher(risk_entrez, term2gene_cgp, universe_entrez), error = function(e) NULL)
ora_protective_cgp <- tryCatch(run_enricher(protective_entrez, term2gene_cgp, universe_entrez), error = function(e) NULL)
ora_risk_cp <- tryCatch(run_enricher(risk_entrez, term2gene_cp, universe_entrez), error = function(e) NULL)
ora_protective_cp <- tryCatch(run_enricher(protective_entrez, term2gene_cp, universe_entrez), error = function(e) NULL)

reactome_gsea <- NULL
reactome_risk <- NULL
reactome_protective <- NULL
if (use_reactome_enrich) {
  reactome_gsea <- tryCatch(
    ReactomePA::gsePathway(
      geneList = rank_vec,
      organism = "human",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      eps = 1e-10,
      verbose = FALSE,
      seed = TRUE
    ),
    error = function(e) NULL
  )
  reactome_risk <- tryCatch(ReactomePA::enrichPathway(gene = risk_entrez, organism = "human", universe = universe_entrez, pAdjustMethod = "BH", pvalueCutoff = 1, minGSSize = min_gs_size, maxGSSize = max_gs_size, readable = TRUE), error = function(e) NULL)
  reactome_protective <- tryCatch(ReactomePA::enrichPathway(gene = protective_entrez, organism = "human", universe = universe_entrez, pAdjustMethod = "BH", pvalueCutoff = 1, minGSSize = min_gs_size, maxGSSize = max_gs_size, readable = TRUE), error = function(e) NULL)
}

tbl_reactome_gsea <- extract_tbl(reactome_gsea) %>%
  mutate(absNES = abs(NES), neglog10_fdr = -log10(pmax(p.adjust, 1e-300))) %>%
  arrange(p.adjust, desc(absNES))

tbl_reactome_risk <- extract_tbl(reactome_risk) %>% mutate(signature_side = "Risk-promoting", neglog10_fdr = -log10(pmax(p.adjust, 1e-300)))
tbl_reactome_protective <- extract_tbl(reactome_protective) %>% mutate(signature_side = "Protective", neglog10_fdr = -log10(pmax(p.adjust, 1e-300)))

tbl_ora_risk <- extract_tbl(ora_risk_cgp) %>% left_join(term2name_cgp, by = c("ID" = "gs_name")) %>% mutate(signature_side = "Risk-promoting")
tbl_ora_protective <- extract_tbl(ora_protective_cgp) %>% left_join(term2name_cgp, by = c("ID" = "gs_name")) %>% mutate(signature_side = "Protective")
tbl_ora_risk_cp <- extract_tbl(ora_risk_cp) %>% left_join(term2name_cp, by = c("ID" = "gs_name")) %>% mutate(signature_side = "Risk-promoting")
tbl_ora_protective_cp <- extract_tbl(ora_protective_cp) %>% left_join(term2name_cp, by = c("ID" = "gs_name")) %>% mutate(signature_side = "Protective")

priority_tier_fun <- function(score) {
  case_when(
    score >= 0.80 ~ "Tier 1",
    score >= 0.60 ~ "Tier 2",
    score >= 0.40 ~ "Tier 3",
    TRUE ~ "Tier 4"
  )
}

candidate_reversal <- tbl_cgp %>%
  filter(reversal_class == "Candidate reversal") %>%
  mutate(
    source_collection = "MSigDB C2:CGP",
    direction_bonus = ifelse(perturbation_set_direction == "DOWN gene set", 0.10, ifelse(perturbation_set_direction == "UP gene set", 0.05, 0.02)),
    evidence_bonus = pmin(matched_drug_tokens_n, 3) * 0.08 +
      pmin(matched_gene_tokens_n, 3) * 0.03 +
      pmin(matched_pathway_tokens_n, 3) * 0.03,
    actionability_bonus = case_when(
      actionability_class == "Directly drug-like" ~ 0.18,
      actionability_class == "Mechanistically actionable" ~ 0.12,
      actionability_class == "Biology-supportive" ~ 0.07,
      actionability_class == "Phenotype-supportive" ~ 0.04,
      actionability_class == "Contextual only" ~ 0.02,
      TRUE ~ 0.01
    ),
    priority_score = safe_rescale(absNES) * 0.35 +
      safe_rescale(neglog10_fdr) * 0.30 +
      safe_rescale(setSize) * 0.05 +
      evidence_bonus + actionability_bonus + direction_bonus,
    priority_tier = priority_tier_fun(priority_score)
  ) %>%
  arrange(desc(priority_score), p.adjust, NES)

concordant_perturbagens <- tbl_cgp %>%
  filter(reversal_class == "Concordant / phenocopying") %>%
  mutate(source_collection = "MSigDB C2:CGP") %>%
  arrange(p.adjust, desc(absNES))

context_summary <- candidate_reversal %>%
  count(context_family, sort = TRUE, name = "n_signatures") %>%
  mutate(prop = n_signatures / sum(n_signatures))

actionability_summary <- candidate_reversal %>%
  count(actionability_class, sort = TRUE, name = "n_signatures") %>%
  mutate(prop = n_signatures / sum(n_signatures))

priority_tier_summary <- candidate_reversal %>%
  count(priority_tier, actionability_class, sort = TRUE, name = "n") %>%
  group_by(priority_tier) %>%
  mutate(prop_within_tier = n / sum(n)) %>%
  ungroup()

hallmark_reversal_context <- tbl_h %>% filter(p.adjust < 0.25) %>% arrange(p.adjust, desc(absNES))
canonical_context <- tbl_cp %>% filter(p.adjust < 0.25) %>% arrange(p.adjust, desc(absNES))
reactome_context <- tbl_reactome_gsea %>% filter(p.adjust < 0.25) %>% arrange(p.adjust, desc(absNES))

signature_overlap_summary <- bind_rows(tbl_ora_risk, tbl_ora_protective) %>%
  mutate(neglog10_fdr = -log10(pmax(p.adjust, 1e-300))) %>%
  arrange(signature_side, p.adjust)

signature_overlap_summary_cp <- bind_rows(tbl_ora_risk_cp, tbl_ora_protective_cp) %>%
  mutate(neglog10_fdr = -log10(pmax(p.adjust, 1e-300))) %>%
  arrange(signature_side, p.adjust)

signature_overlap_summary_reactome <- bind_rows(tbl_reactome_risk, tbl_reactome_protective) %>%
  arrange(signature_side, p.adjust)

extract_leading_edge <- function(gsea_tbl, top_n = 10, max_terms = 10) {
  if (is.null(gsea_tbl) || nrow(gsea_tbl) == 0 || !"core_enrichment" %in% names(gsea_tbl)) return(tibble())
  top_tbl <- gsea_tbl %>%
    filter(reversal_class == "Candidate reversal") %>%
    arrange(p.adjust, desc(absNES)) %>%
    utils::head(max_terms)
  if (nrow(top_tbl) == 0) return(tibble())

  purrr::map_dfr(seq_len(nrow(top_tbl)), function(i) {
    ce <- top_tbl$core_enrichment[i]
    genes <- unique(unlist(strsplit(as.character(ce), "/", fixed = TRUE)))
    genes <- genes[nzchar(genes)]
    if (length(genes) == 0) return(tibble())
    genes <- utils::head(genes, top_n)
    mapped <- anno_select(org.Hs.eg.db, keys = genes, columns = c("SYMBOL", "GENENAME"), keytype = "ENTREZID") %>%
      distinct(ENTREZID, .keep_all = TRUE)
    tibble(
      term_id = top_tbl$ID[i],
      perturbagen = top_tbl$perturbagen[i],
      priority_score = top_tbl$priority_score[i],
      ENTREZID = genes
    ) %>%
      left_join(mapped, by = "ENTREZID") %>%
      mutate(gene_symbol = coalesce(SYMBOL, ENTREZID)) %>%
      dplyr::select(term_id, perturbagen, priority_score, ENTREZID, gene_symbol, GENENAME)
  })
}

leading_edge_tbl <- extract_leading_edge(candidate_reversal, top_n = top_leading_edge, max_terms = min(10, shortlist_n))

family_direction_summary <- candidate_reversal %>%
  count(context_family, perturbation_set_direction, name = "n_signatures") %>%
  group_by(context_family) %>%
  mutate(prop = n_signatures / sum(n_signatures)) %>%
  ungroup()

candidate_shortlist_base <- candidate_reversal %>%
  group_by(perturbagen) %>%
  summarise(
    best_NES = min(NES, na.rm = TRUE),
    best_absNES = max(absNES, na.rm = TRUE),
    best_fdr = min(p.adjust, na.rm = TRUE),
    best_priority = max(priority_score, na.rm = TRUE),
    priority_tier = priority_tier_fun(best_priority),
    context_family = context_family[which.max(priority_score)][1],
    actionability_class = actionability_class[which.max(priority_score)][1],
    is_druglike = any(is_druglike, na.rm = TRUE),
    n_supporting_signatures = n(),
    contrast_support = paste(sort(unique(perturbation_set_direction)), collapse = "; "),
    matched_drug_tokens = matched_drug_tokens[which.max(priority_score)][1],
    matched_gene_tokens = matched_gene_tokens[which.max(priority_score)][1],
    matched_pathway_tokens = matched_pathway_tokens[which.max(priority_score)][1],
    .groups = "drop"
  ) %>%
  mutate(neglog10_fdr = -log10(pmax(best_fdr, 1e-300))) %>%
  arrange(desc(best_priority), best_fdr, best_NES)

candidate_shortlist <- candidate_shortlist_base
if (druglike_only_shortlist) {
  ds <- candidate_shortlist_base %>% filter(is_druglike | actionability_class %in% c("Directly drug-like", "Mechanistically actionable"))
  if (nrow(ds) > 0) candidate_shortlist <- ds
}
candidate_shortlist <- utils::head(candidate_shortlist, shortlist_n)

gene_direction_summary <- gene_tbl %>%
  count(direction, stability_tier, name = "n_genes") %>%
  group_by(direction) %>%
  mutate(prop = n_genes / sum(n_genes)) %>%
  ungroup()

annotation_evidence_summary <- candidate_reversal %>%
  summarise(
    n_reversal = n(),
    n_druglike = sum(is_druglike, na.rm = TRUE),
    n_with_drug_token = sum(matched_drug_tokens_n > 0, na.rm = TRUE),
    n_with_gene_token = sum(matched_gene_tokens_n > 0, na.rm = TRUE),
    n_with_pathway_token = sum(matched_pathway_tokens_n > 0, na.rm = TRUE)
  )

source_summary <- tibble(
  analysis_layer = c(
    "GSEA: perturbagen signatures (C2:CGP)",
    "GSEA: Hallmark",
    "GSEA: canonical pathways (C2:CP)",
    "GSEA: ReactomePA",
    "ORA: risk-promoting genes vs C2:CGP",
    "ORA: protective genes vs C2:CGP",
    "ORA: risk-promoting genes vs C2:CP",
    "ORA: protective genes vs C2:CP",
    "ORA: risk-promoting genes vs ReactomePA",
    "ORA: protective genes vs ReactomePA",
    "Candidate reversal perturbational signatures",
    "Concordant / phenocopying signatures",
    "Shortlisted perturbational nominations",
    "Leading-edge gene rows"
  ),
  n_terms = c(
    nrow(tbl_cgp), nrow(tbl_h), nrow(tbl_cp), nrow(tbl_reactome_gsea),
    nrow(tbl_ora_risk), nrow(tbl_ora_protective), nrow(tbl_ora_risk_cp), nrow(tbl_ora_protective_cp),
    nrow(tbl_reactome_risk), nrow(tbl_reactome_protective),
    nrow(candidate_reversal), nrow(concordant_perturbagens), nrow(candidate_shortlist), nrow(leading_edge_tbl)
  )
)

interpretation_notes <- tribble(
  ~component, ~note,
  "Entity dictionaries", "Drug tokens are derived from KEGG drug names when KEGGREST is available; gene tokens are derived from org.Hs.eg.db SYMBOL and ALIAS mappings; pathway tokens are derived from Reactome, KEGG, and MSigDB canonical collections.",
  "C2:CGP collection", "CGP remains a heterogeneous perturbational collection. Package-driven dictionaries improve evidence provenance but do not convert every CGP term into a bona fide small-molecule perturbation.",
  "Reversal candidates", "Negative NES indicates inverse connectivity to the LUAD ferroptosis-related signature; these terms support perturbational reversal hypotheses.",
  "Actionability class", "Directly drug-like and mechanistically actionable signatures should be prioritized over purely phenotypic or lineage-context signatures for medicinal-chemistry follow-up.",
  "Reactome layer", "ReactomePA is added as an orthogonal pathway-level validation layer to support mechanism-oriented interpretation in a journal-ready workflow."
)

# -----------------------------
# tables
# -----------------------------

safe_write(gene_tbl, file.path(tab_dir, "table_signature_genes_annotated.tsv"))
safe_write(tbl_cgp, file.path(tab_dir, "table_cgp_gsea_full.tsv"))
safe_write(tbl_cp, file.path(tab_dir, "table_canonical_gsea_full.tsv"))
safe_write(tbl_h, file.path(tab_dir, "table_hallmark_gsea_full.tsv"))
safe_write(tbl_reactome_gsea, file.path(tab_dir, "table_reactome_gsea_full.tsv"))
safe_write(candidate_reversal, file.path(tab_dir, "table_reversal_perturbational_candidates.tsv"))
safe_write(concordant_perturbagens, file.path(tab_dir, "table_concordant_phenocopying_signatures.tsv"))
safe_write(candidate_shortlist, file.path(tab_dir, "table_perturbational_nomination_shortlist.tsv"))
safe_write(context_summary, file.path(tab_dir, "table_context_family_summary.tsv"))
safe_write(actionability_summary, file.path(tab_dir, "table_actionability_class_summary.tsv"))
safe_write(priority_tier_summary, file.path(tab_dir, "table_priority_tier_summary.tsv"))
safe_write(family_direction_summary, file.path(tab_dir, "table_family_direction_summary.tsv"))
safe_write(tbl_ora_risk, file.path(tab_dir, "table_risk_promoting_c2_overlap.tsv"))
safe_write(tbl_ora_protective, file.path(tab_dir, "table_protective_c2_overlap.tsv"))
safe_write(tbl_ora_risk_cp, file.path(tab_dir, "table_risk_promoting_canonical_overlap.tsv"))
safe_write(tbl_ora_protective_cp, file.path(tab_dir, "table_protective_canonical_overlap.tsv"))
safe_write(tbl_reactome_risk, file.path(tab_dir, "table_risk_promoting_reactome_overlap.tsv"))
safe_write(tbl_reactome_protective, file.path(tab_dir, "table_protective_reactome_overlap.tsv"))
safe_write(signature_overlap_summary, file.path(tab_dir, "table_signature_overlap_summary.tsv"))
safe_write(signature_overlap_summary_cp, file.path(tab_dir, "table_signature_canonical_overlap_summary.tsv"))
safe_write(signature_overlap_summary_reactome, file.path(tab_dir, "table_signature_reactome_overlap_summary.tsv"))
safe_write(hallmark_reversal_context, file.path(tab_dir, "table_hallmark_context.tsv"))
safe_write(canonical_context, file.path(tab_dir, "table_canonical_context.tsv"))
safe_write(reactome_context, file.path(tab_dir, "table_reactome_context.tsv"))
safe_write(leading_edge_tbl, file.path(tab_dir, "table_leading_edge_genes_reversal.tsv"))
safe_write(gene_direction_summary, file.path(tab_dir, "table_gene_direction_summary.tsv"))
safe_write(annotation_evidence_summary, file.path(tab_dir, "table_annotation_evidence_summary.tsv"))
safe_write(source_summary, file.path(tab_dir, "table_source_summary.tsv"))
safe_write(interpretation_notes, file.path(tab_dir, "table_interpretation_notes.tsv"))

# -----------------------------
# figures
# -----------------------------

if (nrow(candidate_reversal) > 0) {
  df <- top_n_df(candidate_reversal, max_reverse_terms)
  df <- df %>% mutate(label = fct_reorder(wrap_label(perturbagen, 38), NES))
  p1 <- ggplot(df, aes(x = NES, y = label, fill = neglog10_fdr)) +
    geom_col(width = 0.72) +
    scale_fill_viridis_c(option = "plasma", direction = -1, name = expression(-log[10](adj.~italic(p)))) +
    labs(
      title = "Perturbational signatures predicted to reverse the LUAD ferroptosis program",
      subtitle = "More negative NES indicates stronger inverse connectivity relative to the disease-associated signature",
      x = "Normalized enrichment score (NES)",
      y = NULL
    ) +
    q1_theme(11) +
    theme(axis.text.y = element_text(size = 8))
  plot_save(p1, file.path(fig_dir, "Fig_R1_signature_reversal_barplot"), width = 11.5, height = 8.5)
}

if (nrow(candidate_reversal) > 0) {
  df <- top_n_df(candidate_reversal, max_reverse_terms)
  df <- df %>% mutate(label = fct_reorder(wrap_label(perturbagen, 38), priority_score))
  p2 <- ggplot(df, aes(x = priority_score, y = label, size = absNES, color = neglog10_fdr, shape = actionability_class)) +
    geom_point(alpha = 0.9) +
    scale_color_viridis_c(option = "magma", direction = -1, name = expression(-log[10](adj.~italic(p)))) +
    scale_size_continuous(name = "|NES|") +
    labs(
      title = "Priority-ranked reversal perturbational signatures",
      subtitle = "Priority score integrates effect magnitude, significance, package-driven annotation evidence, direction, and actionability class",
      x = "Priority score",
      y = NULL,
      shape = "Actionability class"
    ) +
    q1_theme(11)
  plot_save(p2, file.path(fig_dir, "Fig_R2_priority_bubble"), width = 12, height = 8.5)
}

if (nrow(context_summary) > 1) {
  df <- top_n_df(context_summary, 8, arrange_by = "n_signatures")
  df <- df %>% mutate(context_family = fct_reorder(context_family, n_signatures))
  p3 <- ggplot(df, aes(x = n_signatures, y = context_family, fill = n_signatures)) +
    geom_col(width = 0.72) +
    scale_fill_viridis_c(option = "cividis", direction = -1, guide = "none") +
    labs(
      title = "Perturbational context classes represented among reversal signatures",
      x = "Number of significant reversal signatures",
      y = NULL
    ) +
    q1_theme(11)
  plot_save(p3, file.path(fig_dir, "Fig_R3_context_family_summary"), width = 10.5, height = 6)
}

if (nrow(actionability_summary) > 0) {
  df <- actionability_summary %>% mutate(actionability_class = fct_reorder(actionability_class, n_signatures))
  p4 <- ggplot(df, aes(x = n_signatures, y = actionability_class, fill = prop)) +
    geom_col(width = 0.72) +
    scale_fill_viridis_c(option = "inferno", direction = -1, name = "Proportion") +
    labs(
      title = "Actionability classes among reversal signatures",
      subtitle = "Directly drug-like calls are package-backed rather than hard-coded",
      x = "Number of signatures",
      y = NULL
    ) +
    q1_theme(10)
  plot_save(p4, file.path(fig_dir, "Fig_R4_actionability_summary"), width = 10.5, height = 6)
}

if (nrow(signature_overlap_summary) > 0) {
  df <- signature_overlap_summary %>% filter(p.adjust < 0.25)
  df <- top_n_df(df, max_overlap_terms, arrange_by = "neglog10_fdr")
  if (nrow(df) > 0) {
    df <- df %>%
      mutate(Description = ifelse(is.na(Description) | Description == "", ID, Description),
             Description = fct_reorder(wrap_label(Description, 40), neglog10_fdr),
             signature_side = factor(signature_side, levels = c("Risk-promoting", "Protective")))
    p5 <- ggplot(df, aes(x = neglog10_fdr, y = Description, fill = signature_side)) +
      geom_col(width = 0.72) +
      facet_wrap(~ signature_side, scales = "free_y") +
      scale_fill_manual(values = c("Risk-promoting" = "#C0392B", "Protective" = "#2980B9")) +
      labs(
        title = "Directional overlap with perturbational C2 gene sets",
        subtitle = "ORA highlights differential alignment of risk-promoting and protective signature genes",
        x = expression(-log[10](adj.~italic(p))),
        y = NULL,
        fill = NULL
      ) +
      q1_theme(10) +
      theme(axis.text.y = element_text(size = 7.5))
    plot_save(p5, file.path(fig_dir, "Fig_R5_directional_overlap_facets"), width = 13, height = 8)
  }
}

if (nrow(canonical_context) > 0) {
  df <- top_n_df(canonical_context, max_canonical_terms)
  df <- df %>% mutate(Description = ifelse(is.na(Description) | Description == "", ID, Description),
                      Description = fct_reorder(wrap_label(Description, 40), NES))
  p6 <- ggplot(df, aes(x = NES, y = Description, fill = neglog10_fdr)) +
    geom_col(width = 0.72) +
    scale_fill_viridis_c(option = "viridis", direction = -1, name = expression(-log[10](adj.~italic(p)))) +
    labs(
      title = "Canonical pathway context of the LUAD ferroptosis signature",
      subtitle = "C2 canonical pathways contextualize the biology underlying perturbational reversal hypotheses",
      x = "Normalized enrichment score (NES)",
      y = NULL
    ) +
    q1_theme(11) +
    theme(axis.text.y = element_text(size = 8))
  plot_save(p6, file.path(fig_dir, "Fig_R6_canonical_context_barplot"), width = 12, height = 8)
}

if (nrow(candidate_reversal) > 0) {
  df <- top_n_df(candidate_reversal, max_reverse_terms)
  p7 <- ggplot(df, aes(x = NES, y = neglog10_fdr, color = context_family, size = priority_score)) +
    geom_point(alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    labs(
      title = "Effect-size versus significance landscape for reversal signatures",
      subtitle = "Upper-left region indicates stronger and more significant inverse connectivity",
      x = "Normalized enrichment score (NES)",
      y = expression(-log[10](adj.~italic(p))),
      color = "Context family",
      size = "Priority score"
    ) +
    q1_theme(10)
  plot_save(p7, file.path(fig_dir, "Fig_R7_reversal_landscape"), width = 11, height = 8)
}

if (nrow(priority_tier_summary) > 0 && length(unique(priority_tier_summary$actionability_class)) > 1) {
  p8 <- ggplot(priority_tier_summary, aes(x = priority_tier, y = n, fill = actionability_class)) +
    geom_col(width = 0.72) +
    labs(
      title = "Priority-tier composition of reversal perturbational signatures",
      subtitle = "Priority tiers summarize integrated evidence strength across actionability classes",
      x = "Priority tier",
      y = "Number of signatures",
      fill = "Actionability class"
    ) +
    q1_theme(10)
  plot_save(p8, file.path(fig_dir, "Fig_R8_priority_tier_stacked_bar"), width = 11, height = 7)
}

if (nrow(leading_edge_tbl) > 0) {
  led_counts <- leading_edge_tbl %>%
    count(perturbagen, gene_symbol, name = "n") %>%
    group_by(perturbagen) %>%
    mutate(rank_within = dplyr::row_number(dplyr::desc(n))) %>%
    filter(rank_within <= min(8, top_leading_edge)) %>%
    ungroup()
  led_counts$perturbagen <- fct_reorder(led_counts$perturbagen, led_counts$n, .fun = max)
  p9 <- ggplot(led_counts, aes(x = gene_symbol, y = perturbagen, fill = n)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "mako", direction = -1, name = "Count") +
    labs(
      title = "Leading-edge gene reuse across top reversal signatures",
      subtitle = "Genes repeatedly appearing in leading edges may indicate convergent reversal biology",
      x = "Leading-edge gene",
      y = NULL
    ) +
    q1_theme(10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
          axis.text.y = element_text(size = 8))
  plot_save(p9, file.path(fig_dir, "Fig_R9_leading_edge_heatmap"), width = 12, height = 7)
}

if (nrow(candidate_shortlist) > 0) {
  df <- candidate_shortlist %>% mutate(label = fct_reorder(wrap_label(perturbagen, 38), best_priority))
  p10 <- ggplot(df, aes(x = best_priority, y = label, color = actionability_class)) +
    geom_segment(aes(x = 0, xend = best_priority, y = label, yend = label), linewidth = 0.8, alpha = 0.7) +
    geom_point(aes(size = neglog10_fdr), alpha = 0.95) +
    labs(
      title = "Shortlisted perturbational nominations",
      subtitle = "Shortlist integrates inverse connectivity, significance, support, and package-derived annotation evidence",
      x = "Best priority score",
      y = NULL,
      color = "Actionability class",
      size = expression(-log[10](adj.~italic(p)))
    ) +
    q1_theme(10)
  plot_save(p10, file.path(fig_dir, "Fig_R10_shortlist_lollipop"), width = 11.5, height = 8.5)
}

if (nrow(family_direction_summary) > 0 && length(unique(family_direction_summary$context_family)) > 1) {
  p11 <- ggplot(family_direction_summary, aes(x = perturbation_set_direction, y = context_family, fill = n_signatures)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "rocket", direction = -1) +
    labs(
      title = "Directionality distribution across perturbational context families",
      x = "Perturbation-set direction",
      y = NULL,
      fill = "No. signatures"
    ) +
    q1_theme(10)
  plot_save(p11, file.path(fig_dir, "Fig_R11_family_direction_heatmap"), width = 9.5, height = 6)
}

if (nrow(candidate_reversal) > 0) {
  df <- candidate_reversal %>%
    count(actionability_class, is_druglike, name = "n") %>%
    mutate(actionability_class = fct_reorder(actionability_class, n))
  p12 <- ggplot(df, aes(x = n, y = actionability_class, fill = is_druglike)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c("TRUE" = "#2C7FB8", "FALSE" = "#BDBDBD"), name = "Drug-like flag") +
    labs(
      title = "Drug-like versus non-drug-like composition of reversal signatures",
      x = "Number of signatures",
      y = NULL
    ) +
    q1_theme(10)
  plot_save(p12, file.path(fig_dir, "Fig_R12_druglike_composition"), width = 10.5, height = 6)
}

if (nrow(reactome_context) > 0) {
  df <- top_n_df(reactome_context, max_canonical_terms)
  df <- df %>% mutate(Description = ifelse(is.na(Description) | Description == "", ID, Description),
                      Description = fct_reorder(wrap_label(Description, 40), NES))
  p13 <- ggplot(df, aes(x = NES, y = Description, fill = neglog10_fdr)) +
    geom_col(width = 0.72) +
    scale_fill_viridis_c(option = "turbo", direction = -1, name = expression(-log[10](adj.~italic(p)))) +
    labs(
      title = "Reactome pathway context of the LUAD ferroptosis signature",
      subtitle = "Orthogonal pathway validation layer derived from ReactomePA",
      x = "Normalized enrichment score (NES)",
      y = NULL
    ) +
    q1_theme(11) +
    theme(axis.text.y = element_text(size = 8))
  plot_save(p13, file.path(fig_dir, "Fig_R13_reactome_context_barplot"), width = 12, height = 8)
}

# -----------------------------
# metadata + workbook
# -----------------------------

run_metadata <- tibble(
  key = c(
    "project_root", "gene_file", "script_name", "generated_at", "n_signature_genes",
    "n_ranked_genes", "n_cgp_terms", "n_canonical_terms", "n_hallmark_terms", "n_reactome_terms",
    "n_candidate_reversal", "n_concordant", "n_shortlist", "n_leading_edge_rows",
    "has_openxlsx", "has_keggrest", "msigdbr_version", "druglike_only_shortlist",
    "gene_dictionary_n", "drug_dictionary_n", "pathway_dictionary_n"
  ),
  value = c(
    root,
    normalizePath(gene_file, winslash = "/", mustWork = FALSE),
    script_basename,
    as.character(Sys.time()),
    as.character(nrow(gene_tbl)),
    as.character(length(rank_vec)),
    as.character(nrow(tbl_cgp)),
    as.character(nrow(tbl_cp)),
    as.character(nrow(tbl_h)),
    as.character(nrow(tbl_reactome_gsea)),
    as.character(nrow(candidate_reversal)),
    as.character(nrow(concordant_perturbagens)),
    as.character(nrow(candidate_shortlist)),
    as.character(nrow(leading_edge_tbl)),
    as.character(has_openxlsx),
    as.character(has_keggrest),
    as.character(utils::packageVersion("msigdbr")),
    as.character(druglike_only_shortlist),
    as.character(nrow(gene_dict)),
    as.character(nrow(drug_dict)),
    as.character(nrow(pathway_dict))
  )
)
safe_write(run_metadata, file.path(log_dir, "run_metadata.tsv"))
manifest <- tibble(path = list.files(out_dir, recursive = TRUE, full.names = FALSE))
safe_write(manifest, file.path(log_dir, "output_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

if (has_openxlsx) {
  wb <- openxlsx::createWorkbook()
  add_sheet <- function(name, df) {
    if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, as.data.frame(df))
  }
  add_sheet("genes", gene_tbl)
  add_sheet("cgp_gsea", tbl_cgp)
  add_sheet("canonical_gsea", tbl_cp)
  add_sheet("hallmark_gsea", tbl_h)
  add_sheet("reactome_gsea", tbl_reactome_gsea)
  add_sheet("reversal_candidates", candidate_reversal)
  add_sheet("concordant_candidates", concordant_perturbagens)
  add_sheet("shortlist", candidate_shortlist)
  add_sheet("context_summary", context_summary)
  add_sheet("actionability_summary", actionability_summary)
  add_sheet("priority_tiers", priority_tier_summary)
  add_sheet("family_direction", family_direction_summary)
  add_sheet("risk_overlap_cgp", tbl_ora_risk)
  add_sheet("protective_overlap_cgp", tbl_ora_protective)
  add_sheet("risk_overlap_cp", tbl_ora_risk_cp)
  add_sheet("protective_overlap_cp", tbl_ora_protective_cp)
  add_sheet("risk_overlap_reactome", tbl_reactome_risk)
  add_sheet("protective_overlap_reactome", tbl_reactome_protective)
  add_sheet("hallmark_context", hallmark_reversal_context)
  add_sheet("canonical_context", canonical_context)
  add_sheet("reactome_context", reactome_context)
  add_sheet("leading_edge", leading_edge_tbl)
  add_sheet("annotation_evidence", annotation_evidence_summary)
  add_sheet("source_summary", source_summary)
  add_sheet("interpretation_notes", interpretation_notes)
  openxlsx::saveWorkbook(wb, file.path(tab_dir, "signature_reversal_results_package_driven.xlsx"), overwrite = TRUE)
}

cat("[DONE] Package-driven signature-reversal / perturbational-connectivity analysis completed.\n")
