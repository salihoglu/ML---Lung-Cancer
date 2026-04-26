################################################################################
##                                                                            ##
##   COMPREHENSIVE PATHWAY ENRICHMENT ANALYSIS —                              ##
##   Ferroptosis-Related Prognostic Genes in Lung Adenocarcinoma (LUAD)       ##
##                                                                            ##
##   ANALYSIS MODULES:                                                        ##
##    A.  Gene annotation & ID mapping                                        ##
##    B.  ORA  — GO (BP / MF / CC)                                            ##
##    C.  ORA  — KEGG pathways                                                ##
##    D.  ORA  — Reactome pathways                                            ##
##    E.  ORA  — MSigDB Hallmark                                              ##
##    F.  ORA  — MSigDB C2 (BIOCARTA + PID + WikiPathways + KEGG_LEGACY)     ##
##    G.  ORA  — MSigDB C7 Immunologic (ImmuneSigDB)                         ##
##    H.  ORA  — MSigDB C8 Cell-type signatures                              ##
##    I.  ORA  — Disease Ontology (DOSE)                                      ##
##    J.  ORA  — Network of Cancer Genes (NCG)                                ##
##    K.  ORA  — DisGeNET                                                     ##
##    L.  GSEA — Hallmark (meta-Z ranked)                                     ##
##    M.  GSEA — GO-BP (meta-Z ranked)                                        ##
##    N.  GSEA — KEGG (meta-Z ranked)                                         ##
##    O.  GSEA — C2 Canonical (meta-Z ranked)                                 ##
##    P.  Direction-stratified GO/KEGG (risk-promoting vs protective)         ##
##    Q.  compareCluster — multi-group pathway comparison                     ##
##    R.  Visualisations: dotplot, barplot, cnetplot, emapplot, treeplot,     ##
##        ridgeplot, upset, heatmap, lollipop, bubble, GSEA running-score     ##
##    S.  Excel export (all results, multi-sheet)                             ##
##    T.  RData checkpoint                                                    ##
##                                                                            ##
##   Author  : Rana Salihoglu                                                 ##
##   Date    : 2026-04-12                                                     ##
##                                                                            ##
################################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# 0.  PACKAGE INSTALLATION & LOADING
# ═══════════════════════════════════════════════════════════════════════════════
bioc_pkgs <- c(
  "clusterProfiler", "enrichplot", "DOSE",
  "org.Hs.eg.db", "AnnotationDbi",
  "ReactomePA",
  "pathview",
  "ComplexHeatmap", "circlize",
  "BiocParallel", "limma"
)
cran_pkgs <- c(
  "msigdbr",
  "ggplot2", "ggrepel", "patchwork", "ggforce", "ggridges",
  "dplyr", "tidyr", "stringr", "purrr", "tibble", "forcats",
  "openxlsx",
  "RColorBrewer", "viridis", "scales", "colorspace",
  "igraph", "ggraph",
  "pheatmap",
  "UpSetR",
  "cowplot", "gridExtra"
)

install_if_missing <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("  [INSTALL] ", p)
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else {
        install.packages(p, repos = "https://cloud.r-project.org")
      }
    }
  }
}

install_if_missing(bioc_pkgs, bioc = TRUE)
install_if_missing(cran_pkgs, bioc = FALSE)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ReactomePA)
  library(msigdbr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ggridges)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(forcats)
  library(openxlsx)
  library(RColorBrewer)
  library(viridis)
  library(scales)
  library(colorspace)
  library(igraph)
  library(ggraph)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(UpSetR)
  library(cowplot)
  library(BiocParallel)
})

register(MulticoreParam(workers = max(1L, parallel::detectCores() - 1L)))
set.seed(2025)

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  PATHS & OUTPUT DIRECTORIES
# ═══════════════════════════════════════════════════════════════════════════════
BASE_DIR   <- "/home/rana/Desktop/PROJECT_16/test"
SCRIPT_DIR <- file.path(BASE_DIR, "bioinf_scripts/scripts")
OUT_DIR    <- file.path(SCRIPT_DIR, "pathway_enrichment_Q1_v2")

subdirs <- c("figures/ORA", "figures/GSEA", "figures/comparative",
             "figures/network", "figures/overview",
             "tables", "rdata", "pathview")
for (d in subdirs) dir.create(file.path(OUT_DIR, d),
                               showWarnings = FALSE, recursive = TRUE)

fig_ora  <- file.path(OUT_DIR, "figures/ORA")
fig_gsea <- file.path(OUT_DIR, "figures/GSEA")
fig_comp <- file.path(OUT_DIR, "figures/comparative")
fig_net  <- file.path(OUT_DIR, "figures/network")
fig_ov   <- file.path(OUT_DIR, "figures/overview")

message("Output root: ", OUT_DIR)

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  GLOBAL THEME & SAVE HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

## ── Publication theme ────────────────────────────────────────────────────────
q1_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title        = element_text(face = "bold", size = base_size + 2,
                                       hjust = 0, color = "#1A1A1A"),
      plot.subtitle     = element_text(size = base_size - 0.5, color = "#4D4D4D",
                                       hjust = 0, margin = margin(b = 6)),
      plot.caption      = element_text(size = 7.5, color = "#888888",
                                       hjust = 1, face = "italic",
                                       margin = margin(t = 8)),
      axis.title        = element_text(face = "bold", size = base_size),
      axis.text         = element_text(size = base_size - 1, color = "#2B2B2B"),
      axis.line         = element_line(linewidth = 0.45, color = "black"),
      axis.ticks        = element_line(linewidth = 0.35, color = "black"),
      panel.grid.major.x = element_line(linetype = "dashed",
                                         color = "#E0E0E0", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      legend.title      = element_text(face = "bold", size = base_size - 1),
      legend.text       = element_text(size = base_size - 2),
      legend.key.size   = unit(0.38, "cm"),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.margin     = margin(4, 4, 4, 4),
      plot.margin       = margin(12, 18, 12, 12),
      strip.background  = element_rect(fill = "#EFEFEF", color = "grey60",
                                        linewidth = 0.4),
      strip.text        = element_text(face = "bold", size = base_size - 1,
                                        color = "#222222")
    )
}

## Colour constants
COL_RISK <- c("Risk-promoting" = "#C0392B", "Protective" = "#2980B9")
COL_ONT  <- c("Biological Process" = "#4472C4",
               "Molecular Function" = "#ED7D31",
               "Cellular Component" = "#70AD47")

## ── Save function (PDF + TIFF 300 dpi always) ────────────────────────────────
save_fig <- function(p, name, subdir = fig_ora,
                     width = 10, height = 7, dpi = 300) {
  if (is.null(p)) return(invisible(NULL))
  pdf_path  <- file.path(subdir, paste0(name, ".pdf"))
  tiff_path <- file.path(subdir, paste0(name, ".tiff"))
  suppressMessages({
    ggplot2::ggsave(pdf_path,  plot = p, width = width, height = height)
    ggplot2::ggsave(tiff_path, plot = p, width = width, height = height,
                    dpi = dpi, compression = "lzw")
  })
  message("  [FIG] ", basename(pdf_path))
  invisible(p)
}

## ── Safe dotplot wrapper ─────────────────────────────────────────────────────
safe_dotplot <- function(res, n = 20, title = "", subtitle = "",
                          caption = "BH-adjusted p < 0.05",
                          palette = "viridis", direction = -1) {
  if (is.null(res) || nrow(res@result[res@result$p.adjust < 0.05, ]) == 0)
    return(NULL)
  dotplot(res, showCategory = n, font.size = 9) +
    scale_color_viridis_c(option = palette, direction = direction,
                           name = "Adj. p-value") +
    labs(title = title, subtitle = subtitle, caption = caption) +
    q1_theme(10) +
    theme(axis.text.y = element_text(size = 8))
}

## ── Ranked barplot helper ────────────────────────────────────────────────────
make_barplot <- function(res, title = "", subtitle = "",
                          caption = "BH-adjusted p < 0.05",
                          n = 25, fill_pal = "rocket",
                          wrap_w = 55) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  df <- res@result %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = n) %>%
    mutate(
      log10p    = -log10(p.adjust),
      wrap_desc = str_wrap(Description, wrap_w),
      wrap_desc = fct_reorder(wrap_desc, log10p)
    )
  if (nrow(df) == 0) return(NULL)
  ggplot(df, aes(x = log10p, y = wrap_desc, fill = log10p)) +
    geom_col(width = 0.72, alpha = 0.93) +
    geom_text(aes(label = paste0("n=", Count)),
              hjust = -0.1, size = 2.7, color = "#444444") +
    scale_fill_viridis_c(option = fill_pal, direction = -1,
                          name = expression(-log[10](p[adj]))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
    labs(x = expression(-log[10](Adjusted~italic(p)~value)),
         y = NULL, title = title, subtitle = subtitle, caption = caption) +
    q1_theme(10) +
    theme(axis.text.y = element_text(size = 8.5))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  LOAD & ANNOTATE GENE SIGNATURE
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[A] Loading gene signature ──────────────────────────────────────────")

## Locate CSV (multiple fallback paths)
gene_csv_candidates <- c(
  file.path(BASE_DIR, "bioinf_scripts/scripts",
            "possible_prognostic_genes_FULLTRAIN.csv"),
  file.path(BASE_DIR, "possible_prognostic_genes_FULLTRAIN.csv"),
  list.files(BASE_DIR, "possible_prognostic_genes_FULLTRAIN.csv",
             recursive = TRUE, full.names = TRUE)
)
gene_csv <- gene_csv_candidates[file.exists(gene_csv_candidates)][1]
if (is.na(gene_csv)) stop("Gene CSV not found — set path manually.")

gene_df <- read.csv(gene_csv, stringsAsFactors = FALSE) %>%
  mutate(
    direction = ifelse(z_meta > 0, "Risk-promoting", "Protective"),
    bag_tier  = case_when(
      bag_frac == 1.00 ~ "Core (8/8 cohorts)",
      bag_frac >= 0.75 ~ "Stable (6/8 cohorts)",
      TRUE             ~ "Variable"
    ),
    abs_z = abs(z_meta)
  ) %>%
  arrange(desc(abs_z))

message("  Genes loaded: ", nrow(gene_df))
message("  Risk-promoting (z>0): ",  sum(gene_df$direction == "Risk-promoting"))
message("  Protective     (z<0): ",  sum(gene_df$direction == "Protective"))

## ── ID mapping ───────────────────────────────────────────────────────────────
id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(gene_df$gene),
  columns = c("ENTREZID","SYMBOL","GENENAME","ENSEMBL","UNIPROT"),
  keytype = "SYMBOL"
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_df <- gene_df %>%
  left_join(id_map, by = c("gene" = "SYMBOL"))

gene_entrez <- na.omit(gene_df$ENTREZID)
message("  Mapped to Entrez: ", length(gene_entrez), " / ", nrow(gene_df))

## ── GSEA ranked vector (meta-Z; named by Entrez) ─────────────────────────────
gsea_vec <- gene_df %>%
  filter(!is.na(ENTREZID)) %>%
  arrange(desc(z_meta)) %>%
  { setNames(.$z_meta, .$ENTREZID) }

## ── Subset vectors ───────────────────────────────────────────────────────────
entrez_pos <- gene_df %>% filter(direction == "Risk-promoting",
                                  !is.na(ENTREZID)) %>% pull(ENTREZID)
entrez_neg <- gene_df %>% filter(direction == "Protective",
                                  !is.na(ENTREZID)) %>% pull(ENTREZID)
entrez_core <- gene_df %>% filter(bag_frac == 1.0,
                                   !is.na(ENTREZID)) %>% pull(ENTREZID)

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  GENE UNIVERSE
# ═══════════════════════════════════════════════════════════════════════════════
universe_file <- file.path(SCRIPT_DIR,
  "luad_r_analysis_bundle/q1_bioinfo_results/tables/candidate_gene_universe.tsv")

if (file.exists(universe_file)) {
  uni_sym <- read.table(universe_file, header = TRUE,
                         sep = "\t")[[1]]
  universe_entrez <- AnnotationDbi::select(
    org.Hs.eg.db, keys = uni_sym,
    columns = "ENTREZID", keytype = "SYMBOL"
  ) %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID)
  message("  Custom universe: ", length(universe_entrez), " genes")
} else {
  universe_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")
  message("  Full-genome universe: ", length(universe_entrez), " genes")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5.  [B]  GO ORA — BP / MF / CC
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[B] GO Enrichment (ORA) ────────────────────────────────────────────")

run_go_ora <- function(genes, ont, universe = universe_entrez) {
  enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE,
    minGSSize     = 10,
    maxGSSize     = 500
  )
}

go_bp <- run_go_ora(gene_entrez, "BP")
go_mf <- run_go_ora(gene_entrez, "MF")
go_cc <- run_go_ora(gene_entrez, "CC")

for (nm in c("go_bp","go_mf","go_cc")) {
  res <- get(nm)
  n   <- nrow(res@result[res@result$p.adjust < 0.05, ])
  message("  ", toupper(nm), ": ", n, " significant terms")
}

## ── 5a. Individual dotplots ───────────────────────────────────────────────────
p_go_bp_dot <- safe_dotplot(go_bp, 25,
  "GO Biological Process Enrichment",
  "Ferroptosis-related prognostic gene signature | LUAD (n=8 cohorts)",
  palette = "plasma")

p_go_mf_dot <- safe_dotplot(go_mf, 20,
  "GO Molecular Function Enrichment",
  "Ferroptosis-related prognostic gene signature | LUAD",
  palette = "viridis")

p_go_cc_dot <- safe_dotplot(go_cc, 15,
  "GO Cellular Component Enrichment",
  "Ferroptosis-related prognostic gene signature | LUAD",
  palette = "cividis")

save_fig(p_go_bp_dot, "Fig_B1_GO_BP_dotplot",   fig_ora, 10, 11)
save_fig(p_go_mf_dot, "Fig_B2_GO_MF_dotplot",   fig_ora, 10,  9)
save_fig(p_go_cc_dot, "Fig_B3_GO_CC_dotplot",   fig_ora, 10,  8)

## ── 5b. Individual barplots ───────────────────────────────────────────────────
save_fig(make_barplot(go_bp, "GO-BP Enrichment", n = 25, fill_pal = "plasma"),
         "Fig_B4_GO_BP_barplot", fig_ora, 11, 11)
save_fig(make_barplot(go_mf, "GO-MF Enrichment", n = 20, fill_pal = "viridis"),
         "Fig_B5_GO_MF_barplot", fig_ora, 11,  9)
save_fig(make_barplot(go_cc, "GO-CC Enrichment", n = 15, fill_pal = "cividis"),
         "Fig_B6_GO_CC_barplot", fig_ora, 11,  8)

## ── 5c. Combined GO barplot (all three ontologies) ────────────────────────────
make_top_n <- function(res, ont_label, n = 10) {
  res@result %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = n) %>%
    mutate(ont = ont_label, log10p = -log10(p.adjust),
           wrap_desc = str_wrap(Description, 42))
}
go_combined_df <- bind_rows(
  make_top_n(go_bp, "Biological Process", 12),
  make_top_n(go_mf, "Molecular Function",  8),
  make_top_n(go_cc, "Cellular Component",  6)
)
if (nrow(go_combined_df) > 0) {
  go_combined_df <- go_combined_df %>%
    mutate(wrap_desc = fct_reorder(wrap_desc, log10p))

  p_go_combined <- ggplot(go_combined_df,
                           aes(x = log10p, y = wrap_desc, fill = ont)) +
    geom_col(width = 0.72, alpha = 0.92) +
    geom_text(aes(label = paste0("n=", Count)),
              hjust = -0.12, size = 2.6, color = "#444444") +
    facet_wrap(~ ont, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = COL_ONT, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
    labs(x = expression(-log[10](Adjusted~italic(p)~value)), y = NULL,
         title    = "Gene Ontology Enrichment — All Ontologies",
         subtitle = "Top significant GO terms across BP, MF and CC",
         caption  = "BH-adjusted p < 0.05; n = signature genes in term") +
    q1_theme(10) +
    theme(axis.text.y = element_text(size = 8))
  save_fig(p_go_combined, "Fig_B7_GO_combined_barplot",
           fig_ora, 11, 15)
}

## ── 5d. GO-BP treeplot (semantic clustering) ─────────────────────────────────
if (nrow(go_bp@result[go_bp@result$p.adjust < 0.05, ]) >= 5) {
  go_bp_sim <- pairwise_termsim(go_bp)
  p_tree <- tryCatch({
    treeplot(go_bp_sim, showCategory = 30, nWords = 3,
             cex_category = 0.8, hclust_method = "ward.D",
             label_format_cladelab = 18,
             offset = rel(2)) +
      labs(title    = "GO-BP Semantic Clustering (Ward's hierarchical)",
           subtitle = "Ferroptosis-related prognostic genes | LUAD") +
      theme(plot.title    = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 9, color = "#555555"))
  }, error = function(e) {
    message("  treeplot skipped (", conditionMessage(e), ")")
    NULL
  })
  save_fig(p_tree, "Fig_B8_GO_BP_treeplot", fig_ora, 15, 10)
}

## ── 5e. GO-BP emapplot ───────────────────────────────────────────────────────
if (exists("go_bp_sim") &&
    nrow(go_bp@result[go_bp@result$p.adjust < 0.05, ]) >= 5) {
  p_emap_go <- tryCatch(
    emapplot(go_bp_sim, showCategory = 35, cex_label_category = 0.65,
             layout = "nicely", repel = TRUE) +
      labs(title    = "GO-BP Enrichment Map",
           subtitle = "Node size = gene count; edge width = term overlap") +
      theme_graph(base_family = "sans") +
      theme(plot.title    = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 9, color = "#444444")),
    error = function(e) { message("  emapplot skipped"); NULL }
  )
  save_fig(p_emap_go, "Fig_B9_GO_BP_emapplot", fig_net, 14, 12)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 6.  [C]  KEGG ORA
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[C] KEGG Enrichment (ORA) ──────────────────────────────────────────")

kegg_res <- enrichKEGG(
  gene          = gene_entrez,
  universe      = universe_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  minGSSize     = 10,
  maxGSSize     = 500
)
kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
message("  KEGG hits: ",
        nrow(kegg_res@result[kegg_res@result$p.adjust < 0.05, ]))

save_fig(safe_dotplot(kegg_res, 20, "KEGG Pathway Enrichment",
                       "Ferroptosis-related prognostic gene signature | LUAD",
                       palette = "inferno"),
         "Fig_C1_KEGG_dotplot", fig_ora, 10, 9)
save_fig(make_barplot(kegg_res, "KEGG Pathway Enrichment",
                       fill_pal = "inferno"),
         "Fig_C2_KEGG_barplot", fig_ora, 11, 9)

## KEGG emapplot
if (nrow(kegg_res@result[kegg_res@result$p.adjust < 0.05, ]) >= 3) {
  kegg_sim <- pairwise_termsim(kegg_res)
  p_emap_kegg <- tryCatch(
    emapplot(kegg_sim, showCategory = 25, cex_label_category = 0.65,
             layout = "nicely") +
      labs(title = "KEGG Pathway Enrichment Map") +
      theme_graph(base_family = "sans") +
      theme(plot.title = element_text(face = "bold", size = 12)),
    error = function(e) NULL
  )
  save_fig(p_emap_kegg, "Fig_C3_KEGG_emapplot", fig_net, 13, 11)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  [D]  REACTOME ORA
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[D] Reactome Enrichment (ORA) ──────────────────────────────────────")

reactome_res <- enrichPathway(
  gene          = gene_entrez,
  universe      = universe_entrez,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE,
  minGSSize     = 10,
  maxGSSize     = 500
)
message("  Reactome hits: ",
        nrow(reactome_res@result[reactome_res@result$p.adjust < 0.05, ]))

save_fig(safe_dotplot(reactome_res, 20, "Reactome Pathway Enrichment",
                       "Ferroptosis-related prognostic gene signature | LUAD",
                       palette = "mako"),
         "Fig_D1_Reactome_dotplot", fig_ora, 11, 10)
save_fig(make_barplot(reactome_res, "Reactome Pathway Enrichment",
                       fill_pal = "mako", wrap_w = 60),
         "Fig_D2_Reactome_barplot", fig_ora, 12, 10)

# ═══════════════════════════════════════════════════════════════════════════════
# 8.  [E]  MSigDB HALLMARK ORA
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[E] MSigDB Hallmark ORA ────────────────────────────────────────────")

get_msig <- function(category, subcategory = NULL) {
  args <- list(species = "Homo sapiens", category = category)
  if (!is.null(subcategory)) args$subcategory <- subcategory
  do.call(msigdbr, args) %>%
    dplyr::select(gs_name, entrez_gene)
}

hallmark_t2g <- get_msig("H") %>%
  mutate(gs_name = str_remove(gs_name, "HALLMARK_"))

hallmark_res <- enricher(
  gene          = gene_entrez,
  universe      = universe_entrez,
  TERM2GENE     = hallmark_t2g,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  minGSSize     = 5,
  maxGSSize     = 500
)
if (!is.null(hallmark_res))
  hallmark_res <- setReadable(hallmark_res,
                               OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

message("  Hallmark hits: ",
        if (!is.null(hallmark_res))
          nrow(hallmark_res@result[hallmark_res@result$p.adjust < 0.05, ])
        else 0)

save_fig(safe_dotplot(hallmark_res, 20, "MSigDB Hallmark Enrichment",
                       "Ferroptosis-related prognostic gene signature | LUAD",
                       palette = "turbo"),
         "Fig_E1_Hallmark_dotplot", fig_ora, 11, 9)

## Hallmark bubble chart
if (!is.null(hallmark_res) && nrow(hallmark_res) > 0) {
  hall_bub <- hallmark_res@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(
      GRn = sapply(GeneRatio, function(x) {
        p <- strsplit(x, "/")[[1]]; as.numeric(p[1]) / as.numeric(p[2])
      }),
      log10p    = -log10(p.adjust),
      wrap_desc = str_wrap(str_replace_all(Description, "_", " "), 38),
      wrap_desc = fct_reorder(wrap_desc, GRn)
    )
  p_hall_bub <- ggplot(hall_bub,
                        aes(x = GRn, y = wrap_desc,
                            size = Count, color = log10p)) +
    geom_point(alpha = 0.88) +
    scale_color_viridis_c(option = "inferno", direction = -1,
                           name = expression(-log[10](p[adj]))) +
    scale_size_continuous(name = "Gene count", range = c(4, 11)) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "Gene Ratio", y = NULL,
         title    = "MSigDB Hallmark Enrichment Bubble Chart",
         subtitle = "Bubble size = gene count; colour = significance",
         caption  = "BH-adjusted p < 0.05") +
    q1_theme(10) +
    theme(axis.text.y = element_text(size = 9))
  save_fig(p_hall_bub, "Fig_E2_Hallmark_bubble", fig_ora, 11, 9)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 9.  [F]  MSigDB C2 CANONICAL PATHWAYS ORA
#          (BIOCARTA · PID · WikiPathways · KEGG_LEGACY · REACTOME)
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[F] MSigDB C2 Canonical ORA ────────────────────────────────────────")

c2_subcats <- c("CP:BIOCARTA", "CP:PID",
                 "CP:WIKIPATHWAYS", "CP:KEGG_LEGACY", "CP:REACTOME")

c2_results <- setNames(
  lapply(c2_subcats, function(sc) {
    t2g <- get_msig("C2", sc)
    res <- enricher(
      gene          = gene_entrez,
      universe      = universe_entrez,
      TERM2GENE     = t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20,
      minGSSize     = 5,
      maxGSSize     = 500
    )
    if (!is.null(res))
      res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    n <- if (!is.null(res))
           nrow(res@result[res@result$p.adjust < 0.05, ]) else 0
    message("  ", sc, ": ", n, " hits")
    res
  }),
  c2_subcats
)

## Combined C2 barplot
c2_combined_df <- bind_rows(
  lapply(names(c2_results), function(sc) {
    res <- c2_results[[sc]]
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res@result %>%
      filter(p.adjust < 0.05) %>%
      slice_min(p.adjust, n = 6) %>%
      mutate(subcat = sc, log10p = -log10(p.adjust),
             wrap_desc = str_wrap(Description, 50))
  })
)

if (!is.null(c2_combined_df) && nrow(c2_combined_df) > 0) {
  c2_combined_df <- c2_combined_df %>%
    mutate(wrap_desc = fct_reorder(wrap_desc, log10p))

  p_c2_combined <- ggplot(c2_combined_df,
                           aes(x = log10p, y = wrap_desc,
                               fill = subcat)) +
    geom_col(width = 0.72, alpha = 0.9) +
    geom_text(aes(label = paste0("n=", Count)),
              hjust = -0.12, size = 2.5, color = "#444444") +
    facet_wrap(~ subcat, scales = "free_y", ncol = 1) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = expression(-log[10](Adjusted~italic(p)~value)), y = NULL,
         title    = "MSigDB C2 Canonical Pathway Enrichment",
         subtitle = "BIOCARTA · PID · WikiPathways · KEGG_LEGACY · REACTOME",
         caption  = "BH-adjusted p < 0.05; top 6 terms per sub-collection") +
    q1_theme(9) +
    theme(axis.text.y = element_text(size = 7.5))
  save_fig(p_c2_combined, "Fig_F1_C2_canonical_combined",
           fig_ora, 12, 16)

  ## Per-subcategory figures
  for (sc in names(c2_results)) {
    res <- c2_results[[sc]]
    nm  <- str_replace_all(sc, ":", "_")
    save_fig(
      safe_dotplot(res, 20, paste("MSigDB", sc),
                    "Ferroptosis-related prognostic gene signature | LUAD"),
      paste0("Fig_F_", nm, "_dotplot"), fig_ora, 11, 9
    )
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# 10. [G]  MSigDB C7 IMMUNOLOGIC ORA
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[G] MSigDB C7 Immunologic ORA ──────────────────────────────────────")

c7_t2g <- get_msig("C7", "IMMUNESIGDB")
c7_res  <- enricher(
  gene          = gene_entrez,
  universe      = universe_entrez,
  TERM2GENE     = c7_t2g,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  minGSSize     = 5,
  maxGSSize     = 300
)
if (!is.null(c7_res))
  c7_res <- setReadable(c7_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
message("  C7-ImmuneSigDB hits: ",
        if (!is.null(c7_res)) nrow(c7_res@result[c7_res@result$p.adjust < 0.05, ]) else 0)

save_fig(make_barplot(c7_res, "MSigDB C7 Immunologic Signature Enrichment",
                       subtitle = "ImmuneSigDB gene sets | top 25 significant terms",
                       n = 25, fill_pal = "turbo", wrap_w = 60),
         "Fig_G1_C7_Immunologic_barplot", fig_ora, 13, 11)

# ═══════════════════════════════════════════════════════════════════════════════
# 11. [H]  MSigDB C8 CELL-TYPE SIGNATURES ORA
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[H] MSigDB C8 Cell-type ORA ────────────────────────────────────────")

c8_t2g <- get_msig("C8")
c8_res  <- enricher(
  gene          = gene_entrez,
  universe      = universe_entrez,
  TERM2GENE     = c8_t2g,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  minGSSize     = 5,
  maxGSSize     = 500
)
if (!is.null(c8_res))
  c8_res <- setReadable(c8_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
message("  C8 cell-type hits: ",
        if (!is.null(c8_res)) nrow(c8_res@result[c8_res@result$p.adjust < 0.05, ]) else 0)

save_fig(make_barplot(c8_res, "MSigDB C8 Cell-Type Signature Enrichment",
                       subtitle = "Single-cell-defined cell-type gene sets | LUAD signature",
                       n = 20, fill_pal = "plasma", wrap_w = 58),
         "Fig_H1_C8_celltype_barplot", fig_ora, 12, 10)

# ═══════════════════════════════════════════════════════════════════════════════
# 12. [I]  DISEASE ONTOLOGY ORA (DOSE)
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[I] Disease Ontology ORA ───────────────────────────────────────────")

do_res <- tryCatch(
  enrichDO(
    gene          = gene_entrez,
    universe      = universe_entrez,
    ont           = "DO",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE,
    minGSSize     = 5,
    maxGSSize     = 500
  ),
  error = function(e) { message("  enrichDO error: ", e$message); NULL }
)
message("  Disease Ontology hits: ",
        if (!is.null(do_res)) nrow(do_res@result[do_res@result$p.adjust < 0.05, ]) else 0)

save_fig(safe_dotplot(do_res, 20, "Disease Ontology (DO) Enrichment",
                       "LUAD ferroptosis-related gene signature",
                       palette = "rocket"),
         "Fig_I1_DO_dotplot", fig_ora, 10, 9)
save_fig(make_barplot(do_res, "Disease Ontology Enrichment",
                       fill_pal = "rocket"),
         "Fig_I2_DO_barplot", fig_ora, 11, 9)

# ═══════════════════════════════════════════════════════════════════════════════
# 13. [J]  NETWORK OF CANCER GENES (NCG) ORA
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[J] NCG (Network of Cancer Genes) ORA ─────────────────────────────")

ncg_res <- tryCatch(
  enrichNCG(
    gene          = gene_entrez,
    universe      = universe_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE,
    minGSSize     = 5
  ),
  error = function(e) { message("  enrichNCG error: ", e$message); NULL }
)
message("  NCG hits: ",
        if (!is.null(ncg_res)) nrow(ncg_res@result[ncg_res@result$p.adjust < 0.05, ]) else 0)

save_fig(make_barplot(ncg_res, "Network of Cancer Genes (NCG) Enrichment",
                       fill_pal = "magma"),
         "Fig_J1_NCG_barplot", fig_ora, 11, 8)

# ═══════════════════════════════════════════════════════════════════════════════
# 14. [K]  DisGeNET ORA
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[K] DisGeNET ORA ───────────────────────────────────────────────────")

dgn_res <- tryCatch(
  enrichDGN(
    gene          = gene_entrez,
    universe      = universe_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE,
    minGSSize     = 5
  ),
  error = function(e) { message("  enrichDGN error: ", e$message); NULL }
)
message("  DisGeNET hits: ",
        if (!is.null(dgn_res)) nrow(dgn_res@result[dgn_res@result$p.adjust < 0.05, ]) else 0)

save_fig(make_barplot(dgn_res, "DisGeNET Disease-Gene Association Enrichment",
                       fill_pal = "mako", wrap_w = 58),
         "Fig_K1_DisGeNET_barplot", fig_ora, 12, 9)

# ═══════════════════════════════════════════════════════════════════════════════
# 15. GENE-CONCEPT NETWORK (cnetplot) — ORA results
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[NET] Gene-concept network plots ───────────────────────────────────")

## Named z_meta vector (by gene symbol) for node colouring
z_by_symbol <- setNames(gene_df$z_meta, gene_df$gene)

make_cnet <- function(res, title, subtitle = "", n_cat = 8) {
  if (is.null(res) || nrow(res@result[res@result$p.adjust < 0.05, ]) < 2)
    return(NULL)
  tryCatch({
    cnetplot(
      res,
      showCategory   = n_cat,
      foldChange     = z_by_symbol,
      node_label     = "all",
      cex_label_gene = 0.60,
      cex_label_category = 0.72,
      circular       = FALSE
    ) +
      scale_color_gradient2(
        low = "#2980B9", mid = "grey92", high = "#C0392B",
        midpoint = 0, name = "Meta-Z\n(risk direction)"
      ) +
      labs(title = title, subtitle = subtitle) +
      theme_graph(base_family = "sans") +
      theme(plot.title    = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 9, color = "#444444"),
            legend.position = "right")
  }, error = function(e) { message("  cnetplot skipped: ", e$message); NULL })
}

save_fig(make_cnet(go_bp,        "GO-BP Gene-Concept Network",
                    "Nodes coloured by meta-Z score (risk direction)"),
         "Fig_NET1_cnet_GO_BP", fig_net, 14, 12)
save_fig(make_cnet(kegg_res,     "KEGG Gene-Concept Network",
                    "Nodes coloured by meta-Z score"),
         "Fig_NET2_cnet_KEGG", fig_net, 14, 12)
save_fig(make_cnet(reactome_res, "Reactome Gene-Concept Network",
                    "Nodes coloured by meta-Z score"),
         "Fig_NET3_cnet_Reactome", fig_net, 14, 12)
save_fig(make_cnet(do_res,       "Disease Ontology Gene-Concept Network"),
         "Fig_NET4_cnet_DO", fig_net, 13, 11)

# ═══════════════════════════════════════════════════════════════════════════════
# 16. [L]  GSEA — MSigDB HALLMARK
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[L] GSEA — MSigDB Hallmark ─────────────────────────────────────────")

gsea_hallmark <- GSEA(
  geneList      = gsea_vec,
  TERM2GENE     = hallmark_t2g,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  minGSSize     = 5,
  maxGSSize     = 500,
  seed          = 42,
  eps           = 0,
  nPermSimple   = 10000
)
message("  GSEA Hallmark hits: ",
        nrow(gsea_hallmark@result[gsea_hallmark@result$p.adjust < 0.05, ]))

if (nrow(gsea_hallmark) > 0) {
  ## Dotplot split by direction
  p_gsea_h_dot <- dotplot(gsea_hallmark, showCategory = 20,
                           split = ".sign", font.size = 9) +
    facet_grid(~ .sign) +
    scale_color_gradient2(low = "#2980B9", mid = "grey88",
                          high = "#C0392B", midpoint = 0,
                          name = "NES") +
    scale_y_discrete(labels = function(x) str_replace_all(x, "_", " ")) +
    labs(title    = "GSEA — MSigDB Hallmark Gene Sets",
         subtitle = "Ranked by meta-Z score | LUAD multi-cohort signature",
         caption  = "NES = Normalised Enrichment Score; BH-adjusted p < 0.05") +
    q1_theme(10) +
    theme(axis.text.y = element_text(size = 8))
  save_fig(p_gsea_h_dot, "Fig_L1_GSEA_Hallmark_dotplot", fig_gsea, 13, 9)

  ## Ridge plot
  p_gsea_h_ridge <- ridgeplot(gsea_hallmark, showCategory = 20,
                               fill = "p.adjust") +
    scale_fill_viridis_c(option = "plasma", direction = -1,
                          name = "Adj. p-value") +
    scale_y_discrete(labels = function(x) str_replace_all(x, "_", " ")) +
    labs(title    = "GSEA Hallmark — Enrichment Score Distributions",
         subtitle = "Gene-level contribution density to NES") +
    q1_theme(10) +
    theme(axis.text.y = element_text(size = 8.5))
  save_fig(p_gsea_h_ridge, "Fig_L2_GSEA_Hallmark_ridgeplot", fig_gsea, 12, 9)

  ## Running score plots for top 9 terms
  top_h_ids <- gsea_hallmark@result %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 9) %>% pull(ID)

  if (length(top_h_ids) > 0) {
    pdf(file.path(fig_gsea, "Fig_L3_GSEA_Hallmark_running_score_top9.pdf"),
        width = 11, height = 7)
    for (id in top_h_ids) {
      p <- tryCatch(
        gseaplot2(gsea_hallmark, geneSetID = id,
                  title = str_replace_all(id, "_", " "),
                  color = "#C0392B", pvalue_table = TRUE, base_size = 9),
        error = function(e) NULL
      )
      if (!is.null(p)) print(p)
    }
    dev.off()
    message("  GSEA Hallmark running-score plots saved")
  }

  ## NES bar chart (Q1 style — both directions)
  nes_df <- gsea_hallmark@result %>%
    filter(p.adjust < 0.05) %>%
    arrange(NES) %>%
    mutate(
      term  = str_replace_all(Description, "_", " "),
      term  = str_wrap(term, 38),
      term  = factor(term, levels = unique(term)),
      sign  = ifelse(NES > 0, "Activated", "Suppressed"),
      stars = case_when(p.adjust < 0.001 ~ "***",
                        p.adjust < 0.01  ~ "**",
                        p.adjust < 0.05  ~ "*",
                        TRUE             ~ "")
    )
  if (nrow(nes_df) > 0) {
    p_nes_bar <- ggplot(nes_df, aes(x = NES, y = term, fill = sign)) +
      geom_col(width = 0.72, alpha = 0.9) +
      geom_text(aes(label = stars,
                    x = NES + ifelse(NES > 0, 0.05, -0.05)),
                hjust = ifelse(nes_df$NES > 0, 0, 1),
                size = 3.5, color = "#222222") +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
      scale_fill_manual(
        values = c("Activated" = "#C0392B", "Suppressed" = "#2980B9"),
        name = "Enrichment"
      ) +
      labs(x = "Normalised Enrichment Score (NES)", y = NULL,
           title    = "GSEA Hallmark — NES Summary",
           subtitle = "Activated vs suppressed gene sets in high-risk vs low-risk",
           caption  = "* p.adj<0.05  ** p.adj<0.01  *** p.adj<0.001") +
      q1_theme(10) +
      theme(axis.text.y = element_text(size = 8.5))
    save_fig(p_nes_bar, "Fig_L4_GSEA_Hallmark_NES_bar", fig_gsea, 11, 9)
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# 17. [M]  GSEA — GO-BP
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[M] GSEA — GO-BP ───────────────────────────────────────────────────")

gsea_go_bp <- gseGO(
  geneList      = gsea_vec,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  seed          = 42,
  eps           = 0,
  nPermSimple   = 10000
)
message("  GSEA GO-BP hits: ",
        nrow(gsea_go_bp@result[gsea_go_bp@result$p.adjust < 0.05, ]))

if (nrow(gsea_go_bp) > 0) {
  p_gsea_go_dot <- dotplot(gsea_go_bp, showCategory = 20, split = ".sign",
                            font.size = 9) +
    facet_grid(~ .sign) +
    labs(title    = "GSEA — GO Biological Process",
         subtitle = "Meta-Z ranked | LUAD ferroptosis-related gene signature") +
    q1_theme(10) +
    theme(axis.text.y = element_text(size = 8))
  save_fig(p_gsea_go_dot, "Fig_M1_GSEA_GO_BP_dotplot", fig_gsea, 13, 10)

  ## GO-BP GSEA ridgeplot
  p_gsea_go_ridge <- tryCatch(
    ridgeplot(gsea_go_bp, showCategory = 20, fill = "p.adjust") +
      scale_fill_viridis_c(option = "mako", direction = -1,
                            name = "Adj. p-value") +
      labs(title = "GSEA GO-BP — Enrichment Score Distributions") +
      q1_theme(10) +
      theme(axis.text.y = element_text(size = 8)),
    error = function(e) NULL
  )
  save_fig(p_gsea_go_ridge, "Fig_M2_GSEA_GO_BP_ridgeplot", fig_gsea, 12, 10)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 18. [N]  GSEA — KEGG
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[N] GSEA — KEGG ────────────────────────────────────────────────────")

gsea_kegg <- gseKEGG(
  geneList      = gsea_vec,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  seed          = 42,
  eps           = 0,
  nPermSimple   = 10000
)
gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
message("  GSEA KEGG hits: ",
        nrow(gsea_kegg@result[gsea_kegg@result$p.adjust < 0.05, ]))

if (nrow(gsea_kegg) > 0) {
  save_fig(
    dotplot(gsea_kegg, showCategory = 20, split = ".sign", font.size = 9) +
      facet_grid(~ .sign) +
      labs(title    = "GSEA — KEGG Pathways",
           subtitle = "Meta-Z ranked | LUAD ferroptosis-related gene signature") +
      q1_theme(10) +
      theme(axis.text.y = element_text(size = 8)),
    "Fig_N1_GSEA_KEGG_dotplot", fig_gsea, 13, 9
  )

  ## Running scores for top KEGG
  top_kegg_ids <- gsea_kegg@result %>%
    filter(p.adjust < 0.05) %>% slice_min(p.adjust, n = 6) %>% pull(ID)
  if (length(top_kegg_ids) > 0) {
    pdf(file.path(fig_gsea, "Fig_N2_GSEA_KEGG_running_score_top6.pdf"),
        width = 11, height = 7)
    for (id in top_kegg_ids) {
      p <- tryCatch(
        gseaplot2(gsea_kegg, geneSetID = id,
                  title = gsea_kegg@result$Description[gsea_kegg@result$ID == id],
                  color = "#2980B9", pvalue_table = TRUE, base_size = 9),
        error = function(e) NULL
      )
      if (!is.null(p)) print(p)
    }
    dev.off()
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# 19. [O]  GSEA — C2 CANONICAL (combined)
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[O] GSEA — C2 Canonical ────────────────────────────────────────────")

c2_all_t2g <- get_msig("C2")

gsea_c2 <- GSEA(
  geneList      = gsea_vec,
  TERM2GENE     = c2_all_t2g,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  seed          = 42,
  eps           = 0,
  nPermSimple   = 10000
)
gsea_c2 <- setReadable(gsea_c2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
message("  GSEA C2 hits: ",
        nrow(gsea_c2@result[gsea_c2@result$p.adjust < 0.05, ]))

if (nrow(gsea_c2) > 0) {
  save_fig(
    dotplot(gsea_c2, showCategory = 20, split = ".sign", font.size = 9) +
      facet_grid(~ .sign) +
      labs(title    = "GSEA — MSigDB C2 Canonical Pathways",
           subtitle = "Meta-Z ranked | LUAD ferroptosis-related gene signature") +
      q1_theme(10) +
      theme(axis.text.y = element_text(size = 7.5)),
    "Fig_O1_GSEA_C2_dotplot", fig_gsea, 13, 10
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# 20. [P]  DIRECTION-STRATIFIED ORA
#          Risk-promoting (z>0) vs Protective (z<0) separately
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[P] Direction-stratified ORA ───────────────────────────────────────")

run_go_relaxed <- function(genes, ont = "BP") {
  if (length(genes) < 5) return(NULL)
  enrichGO(
    gene = genes, universe = universe_entrez,
    OrgDb = org.Hs.eg.db, ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.10, qvalueCutoff = 0.30,
    readable = TRUE, minGSSize = 5, maxGSSize = 500
  )
}
run_kegg_relaxed <- function(genes) {
  if (length(genes) < 5) return(NULL)
  res <- enrichKEGG(
    gene = genes, universe = universe_entrez,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.10, qvalueCutoff = 0.30,
    minGSSize = 5, maxGSSize = 500
  )
  if (!is.null(res)) res <- setReadable(res, org.Hs.eg.db, "ENTREZID")
  res
}

go_bp_pos  <- run_go_relaxed(entrez_pos,  "BP")
go_bp_neg  <- run_go_relaxed(entrez_neg,  "BP")
go_mf_pos  <- run_go_relaxed(entrez_pos,  "MF")
go_mf_neg  <- run_go_relaxed(entrez_neg,  "MF")
kegg_pos   <- run_kegg_relaxed(entrez_pos)
kegg_neg   <- run_kegg_relaxed(entrez_neg)

## ── Diverging barplot (GO-BP) ─────────────────────────────────────────────────
make_dir_df <- function(res_pos, res_neg, n = 15) {
  get_top <- function(res, dir) {
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res@result %>%
      filter(p.adjust < 0.10) %>%
      slice_min(p.adjust, n = n) %>%
      mutate(direction = dir,
             log10p    = -log10(p.adjust),
             signed_p  = ifelse(dir == "Risk-promoting", log10p, -log10p),
             wrap_desc = str_wrap(Description, 44))
  }
  bind_rows(get_top(res_pos, "Risk-promoting"),
            get_top(res_neg, "Protective"))
}

dir_go_df <- make_dir_df(go_bp_pos, go_bp_neg, 15)
if (!is.null(dir_go_df) && nrow(dir_go_df) > 0) {
  dir_go_df <- dir_go_df %>%
    arrange(signed_p) %>%
    mutate(wrap_desc = factor(wrap_desc, levels = unique(wrap_desc)))

  p_dir_go <- ggplot(dir_go_df,
                      aes(x = signed_p, y = wrap_desc, fill = direction)) +
    geom_col(width = 0.73, alpha = 0.90) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    scale_fill_manual(values = COL_RISK, name = "Gene direction") +
    scale_x_continuous(labels = function(x) abs(round(x, 1)),
                        name  = expression(-log[10](Adjusted~italic(p)~value))) +
    labs(y = NULL,
         title    = "Direction-Stratified GO-BP Enrichment",
         subtitle = paste0(
           "Risk-promoting (n=", length(entrez_pos), ") vs ",
           "Protective (n=", length(entrez_neg), ") genes"
         ),
         caption  = "BH-adjusted p < 0.10; left = protective; right = risk-promoting") +
    q1_theme(10) +
    theme(axis.text.y = element_text(size = 8))
  save_fig(p_dir_go, "Fig_P1_direction_stratified_GO_BP",
           fig_comp, 13, 11)
}

dir_kegg_df <- make_dir_df(kegg_pos, kegg_neg, 12)
if (!is.null(dir_kegg_df) && nrow(dir_kegg_df) > 0) {
  dir_kegg_df <- dir_kegg_df %>%
    arrange(signed_p) %>%
    mutate(wrap_desc = factor(wrap_desc, levels = unique(wrap_desc)))

  p_dir_kegg <- ggplot(dir_kegg_df,
                        aes(x = signed_p, y = wrap_desc, fill = direction)) +
    geom_col(width = 0.73, alpha = 0.90) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    scale_fill_manual(values = COL_RISK, name = "Gene direction") +
    scale_x_continuous(labels = function(x) abs(round(x, 1)),
                        name  = expression(-log[10](Adjusted~italic(p)~value))) +
    labs(y = NULL,
         title    = "Direction-Stratified KEGG Enrichment",
         subtitle = "Risk-promoting vs Protective genes") +
    q1_theme(10) +
    theme(axis.text.y = element_text(size = 8.5))
  save_fig(p_dir_kegg, "Fig_P2_direction_stratified_KEGG",
           fig_comp, 13, 10)
}

## ── Core-only ORA (bag_frac == 1) ────────────────────────────────────────────
go_bp_core <- run_go_relaxed(entrez_core, "BP")
message("  Core-gene GO-BP hits: ",
        if (!is.null(go_bp_core))
          nrow(go_bp_core@result[go_bp_core@result$p.adjust < 0.10, ]) else 0)
save_fig(make_barplot(go_bp_core,
                       "GO-BP Enrichment — Core Genes Only (8/8 cohorts)",
                       subtitle = paste0("n=", length(entrez_core),
                                         " genes present in all 8 LUAD cohorts"),
                       n = 20, fill_pal = "plasma"),
         "Fig_P3_core_genes_GO_BP", fig_comp, 11, 9)

# ═══════════════════════════════════════════════════════════════════════════════
# 21. [Q]  compareCluster — MULTI-GROUP COMPARISON
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[Q] compareCluster — multi-group comparison ────────────────────────")

## Gene list: risk-promoting / protective / core-only groups
cluster_list <- list(
  "Risk-promoting\n(z>0)"       = entrez_pos,
  "Protective\n(z<0)"           = entrez_neg,
  "Core (8/8 cohorts)"          = entrez_core
)

## compareCluster — GO-BP
cc_go <- compareCluster(
  geneClusters  = cluster_list,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

## compareCluster — KEGG
cc_kegg <- compareCluster(
  geneClusters  = cluster_list,
  fun           = "enrichKEGG",
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20
)
cc_kegg <- setReadable(cc_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

## compareCluster — Hallmark
cc_hall <- compareCluster(
  geneClusters  = cluster_list,
  fun           = "enricher",
  TERM2GENE     = hallmark_t2g,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20
)
if (!is.null(cc_hall))
  cc_hall <- setReadable(cc_hall, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

## Plot helper
cc_dotplot <- function(cc_res, title, subtitle = "") {
  if (is.null(cc_res)) return(NULL)
  dotplot(cc_res, showCategory = 10, font.size = 8.5) +
    scale_color_viridis_c(option = "plasma", direction = -1,
                           name = "Adj. p-value") +
    labs(title = title, subtitle = subtitle,
         caption = "BH-adjusted p < 0.05; dot size = gene ratio") +
    q1_theme(10) +
    theme(axis.text.x = element_text(size = 8.5),
          axis.text.y = element_text(size = 8))
}

save_fig(cc_dotplot(cc_go,   "Multi-group GO-BP Comparison",
                     "Risk-promoting vs Protective vs Core genes"),
         "Fig_Q1_compareCluster_GO_BP",   fig_comp, 13, 11)
save_fig(cc_dotplot(cc_kegg, "Multi-group KEGG Comparison",
                     "Risk-promoting vs Protective vs Core genes"),
         "Fig_Q2_compareCluster_KEGG",    fig_comp, 13, 10)
save_fig(cc_dotplot(cc_hall, "Multi-group Hallmark Comparison",
                     "Risk-promoting vs Protective vs Core genes"),
         "Fig_Q3_compareCluster_Hallmark",fig_comp, 13, 10)

# ═══════════════════════════════════════════════════════════════════════════════
# 22.  UPSET PLOT — gene overlap across databases
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[UPSET] Gene-database membership ──────────────────────────────────")

collect_genes <- function(res, db_label, n_terms = 12) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res@result %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = n_terms) %>%
    dplyr::select(Description, geneID) %>%
    separate_rows(geneID, sep = "/") %>%
    mutate(db = db_label)
}

upset_raw <- bind_rows(
  collect_genes(go_bp,        "GO-BP"),
  collect_genes(go_mf,        "GO-MF"),
  collect_genes(kegg_res,     "KEGG"),
  collect_genes(reactome_res, "Reactome"),
  collect_genes(hallmark_res, "Hallmark"),
  collect_genes(do_res,       "Disease\nOntology"),
  collect_genes(dgn_res,      "DisGeNET")
)

if (!is.null(upset_raw) && nrow(upset_raw) > 0) {
  upset_wide <- upset_raw %>%
    distinct(geneID, db) %>%
    mutate(val = 1) %>%
    pivot_wider(names_from = db, values_from = val, values_fill = 0) %>%
    as.data.frame()

  db_cols <- setdiff(names(upset_wide), "geneID")

  if (length(db_cols) >= 2) {
    db_colors <- brewer.pal(max(3, length(db_cols)), "Set2")[seq_along(db_cols)]
    pdf(file.path(fig_net, "Fig_UPSET_gene_db_membership.pdf"),
        width = 13, height = 7)
    print(upset(
      upset_wide[, db_cols],
      sets            = db_cols,
      order.by        = "freq",
      decreasing      = TRUE,
      mb.ratio        = c(0.6, 0.4),
      text.scale      = c(1.4, 1.2, 1.0, 1.0, 1.2, 0.9),
      mainbar.y.label = "Intersection Size (Genes)",
      sets.x.label    = "Genes per Database",
      point.size      = 3.2,
      line.size       = 0.9,
      sets.bar.color  = db_colors,
      main.bar.color  = "#2F4F8F",
      matrix.color    = "#2F4F8F"
    ))
    dev.off()
    message("  UpSet plot saved")
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# 23.  FERROPTOSIS PATHWAY MEMBERSHIP HEATMAP
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[HEAT] Ferroptosis pathway heatmap ────────────────────────────────")

ferr_pathways <- list(
  "GPX4/GSH axis"           = c("GPX4","GSS","GSR","GGT1","SLC7A11","SLC3A2",
                                  "GCLC","GCLM"),
  "Iron metabolism"          = c("TF","TFRC","FTH1","FTL","ACO1","STEAP3",
                                  "HMOX1","SLC40A1","HAMP","BMP6"),
  "Lipid peroxidation"       = c("ALOX5","ALOX5AP","ALOX12","ALOX15","ACSL4",
                                  "LPCAT3","PTGS2","PLA2G4A"),
  "Oxidative stress"         = c("NOX1","NOX4","CYBB","NCF1","SOD1","CAT",
                                  "PRDX1","TXNRD1","NQO1"),
  "Apoptosis/Caspase"        = c("CASP1","CASP8","FAS","FASLG","TNF","RIPK2",
                                  "CFLAR","CASP3","BAX","BCL2"),
  "Immune checkpoint"        = c("CD274","CTLA4","CD28","CD40LG","LAG3",
                                  "HAVCR2","PDCD1"),
  "T-cell signalling"        = c("IL7R","FYN","PRKCB","BTK","ITM2A","BCL11B",
                                  "IRF8","MEF2C","CD3E","LCK"),
  "Inflammation/Cytokines"   = c("IL1B","CXCL10","TLR2","AIM2","GBP5","IRF1",
                                  "NR3C1","IFNG","TNF","NFKB1"),
  "Autophagy/mTOR"           = c("BECN1","ATG5","ATG7","MTOR","SQSTM1","MAP1LC3B"),
  "p53/Cell cycle"           = c("TP53","CDKN1A","MDM2","RB1","E2F1","CCND1")
)

all_ferr <- unique(unlist(ferr_pathways))
sig_ferr  <- intersect(all_ferr, gene_df$gene)

if (length(sig_ferr) >= 5) {
  mat_ferr <- sapply(names(ferr_pathways), function(pw) {
    as.integer(sig_ferr %in% ferr_pathways[[pw]])
  })
  rownames(mat_ferr) <- sig_ferr

  z_vals    <- gene_df$z_meta[match(sig_ferr, gene_df$gene)]
  dir_vals  <- gene_df$direction[match(sig_ferr, gene_df$gene)]
  tier_vals <- gene_df$bag_tier[match(sig_ferr, gene_df$gene)]

  row_ha <- HeatmapAnnotation(
    which       = "row",
    Direction   = dir_vals,
    Bag_Tier    = tier_vals,
    Meta_Z      = anno_barplot(z_vals, gp = gpar(fill = ifelse(z_vals > 0,
                                                                "#C0392B",
                                                                "#2980B9"))),
    col = list(
      Direction = c("Risk-promoting" = "#C0392B", "Protective" = "#2980B9"),
      Bag_Tier  = c("Core (8/8 cohorts)"   = "#2C3E50",
                    "Stable (6/8 cohorts)"  = "#7F8C8D",
                    "Variable"             = "#BDC3C7")
    ),
    annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
  )

  col_fun <- colorRamp2(c(0, 1), c("white", "#2166AC"))

  pdf(file.path(fig_ov, "Fig_HEAT_ferroptosis_pathway_membership.pdf"),
      width = 14, height = max(6, length(sig_ferr) * 0.35 + 4))
  ht <- Heatmap(
    mat_ferr,
    name              = "Member",
    col               = col_fun,
    right_annotation  = row_ha,
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    row_names_gp      = gpar(fontsize = 9),
    column_names_gp   = gpar(fontsize = 9, fontface = "bold"),
    column_names_rot  = 35,
    rect_gp           = gpar(col = "grey85", lwd = 0.5),
    cell_fun          = function(j, i, x, y, width, height, fill) {
      if (mat_ferr[i, j] == 1)
        grid.text("●", x, y, gp = gpar(fontsize = 8, col = "white"))
    },
    column_title      = "Ferroptosis Pathway Gene Membership",
    column_title_gp   = gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(title = "Member",
                                labels = c("No","Yes"),
                                at = c(0, 1))
  )
  draw(ht)
  dev.off()
  message("  Ferroptosis heatmap saved")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 24.  GENE SIGNATURE OVERVIEW — LOLLIPOP + ANNOTATION STRIPS
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[OV] Gene signature overview ───────────────────────────────────────")

gene_ov <- gene_df %>%
  arrange(z_meta) %>%
  mutate(
    gene      = factor(gene, levels = gene),
    risk_cat  = case_when(
      z_meta >  3 ~ "Strong risk (z>3)",
      z_meta >  0 ~ "Risk (0<z≤3)",
      z_meta > -3 ~ "Protective (-3≤z<0)",
      TRUE         ~ "Strong protective (z<-3)"
    ),
    risk_cat  = factor(risk_cat,
                       levels = c("Strong risk (z>3)","Risk (0<z≤3)",
                                  "Protective (-3≤z<0)","Strong protective (z<-3)"))
  )

col_risk_cat <- c(
  "Strong risk (z>3)"        = "#922B21",
  "Risk (0<z≤3)"             = "#E59866",
  "Protective (-3≤z<0)"      = "#7FB3D3",
  "Strong protective (z<-3)" = "#1A5276"
)

p_lollipop <- ggplot(gene_ov,
                      aes(x = z_meta, y = gene, color = risk_cat)) +
  geom_segment(aes(x = 0, xend = z_meta, y = gene, yend = gene),
               linewidth = 0.65, alpha = 0.75) +
  geom_point(aes(size = bag_frac), alpha = 0.92, shape = 19) +
  geom_vline(xintercept = 0, linetype = "dashed",
             linewidth  = 0.45, color = "#555555") +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted",
             linewidth  = 0.35, color = "#AAAAAA") +
  scale_color_manual(values = col_risk_cat, name = "Effect category") +
  scale_size_continuous(
    name   = "Bag fraction\n(stability)",
    breaks = c(0.75, 1.0),
    labels = c("0.75 (6/8)", "1.00 (8/8)"),
    range  = c(2, 5)
  ) +
  labs(
    x        = "Meta-Z Score (prognostic direction across 8 LUAD cohorts)",
    y        = NULL,
    title    = "Ferroptosis-Related Prognostic Gene Signature",
    subtitle = paste0(nrow(gene_ov), " ML-selected genes | 8 GEO cohorts | ",
                      "SHAP + WGCNA validated"),
    caption  = "Point size = bootstrap bag fraction; dashed = 0; dotted = ±2"
  ) +
  q1_theme(10) +
  theme(axis.text.y    = element_text(size = 7.5),
        legend.position = "right")
save_fig(p_lollipop, "Fig_OV1_gene_signature_lollipop",
         fig_ov, 11, 15)

## Z-score volcano-style scatter
p_zscore_scatter <- ggplot(gene_df,
                            aes(x = abs_z, y = bag_frac,
                                color = direction, size = abs_z,
                                label = gene)) +
  geom_point(alpha = 0.80) +
  geom_text_repel(
    data   = gene_df %>% filter(abs_z > 2 | bag_frac == 1.0),
    size   = 2.8, max.overlaps = 20, seed = 42,
    segment.color = "grey60", segment.size = 0.3
  ) +
  geom_hline(yintercept = 0.875, linetype = "dashed",
             color = "#888888", linewidth = 0.4) +
  geom_vline(xintercept = 2, linetype = "dashed",
             color = "#888888", linewidth = 0.4) +
  scale_color_manual(values = c("Risk-promoting" = "#C0392B",
                                 "Protective" = "#2980B9"),
                     name = "Direction") +
  scale_size_continuous(range = c(2, 7), guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.70, 1.02)) +
  labs(
    x        = "|Meta-Z score| (effect magnitude)",
    y        = "Bag fraction (bootstrap stability)",
    title    = "Gene Stability vs Effect Magnitude",
    subtitle = "Gene selection across 8 LUAD cohorts",
    caption  = "Dashed lines: |z| > 2 and bag fraction > 0.875"
  ) +
  q1_theme(10)
save_fig(p_zscore_scatter, "Fig_OV2_stability_vs_effect",
         fig_ov, 9, 7)

# ═══════════════════════════════════════════════════════════════════════════════
# 25.  COMPREHENSIVE MULTI-DB DOTPLOT PANEL (patchwork)
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[PANEL] Multi-database summary panel ──────────────────────────────")

make_mini_bar <- function(res, db_label, n = 8, pal = "viridis") {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  df <- res@result %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = n) %>%
    mutate(log10p    = -log10(p.adjust),
           wrap_desc = str_wrap(Description, 35),
           wrap_desc = fct_reorder(wrap_desc, log10p))
  if (nrow(df) == 0) return(NULL)
  ggplot(df, aes(x = log10p, y = wrap_desc, fill = log10p)) +
    geom_col(width = 0.70, alpha = 0.90) +
    scale_fill_viridis_c(option = pal, direction = -1, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = expression(-log[10](p[adj])), y = NULL,
         title = db_label) +
    theme_classic(base_size = 8) +
    theme(
      plot.title   = element_text(face = "bold", size = 9),
      axis.text.y  = element_text(size = 7),
      axis.text.x  = element_text(size = 7),
      plot.margin  = margin(6, 8, 6, 6)
    )
}

panels <- list(
  make_mini_bar(go_bp,        "GO — Biological Process", pal = "plasma"),
  make_mini_bar(go_mf,        "GO — Molecular Function", pal = "viridis"),
  make_mini_bar(kegg_res,     "KEGG Pathways",           pal = "inferno"),
  make_mini_bar(reactome_res, "Reactome",                pal = "mako"),
  make_mini_bar(hallmark_res, "MSigDB Hallmark",         pal = "turbo"),
  make_mini_bar(do_res,       "Disease Ontology",        pal = "rocket"),
  make_mini_bar(dgn_res,      "DisGeNET",                pal = "magma"),
  make_mini_bar(c7_res,       "MSigDB C7 Immunologic",   pal = "cividis")
)
panels <- Filter(Negate(is.null), panels)

if (length(panels) >= 2) {
  ncols   <- min(4, length(panels))
  panel_p <- wrap_plots(panels, ncol = ncols) +
    plot_annotation(
      title    = "Pathway Enrichment Summary — All Databases",
      subtitle = paste0(
        "Ferroptosis-related prognostic gene signature | LUAD | ",
        nrow(gene_df), " genes across 8 cohorts"
      ),
      caption  = "BH-adjusted p < 0.05; top 8 terms per database",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "#444444"),
        plot.caption  = element_text(size = 8, color = "#888888",
                                      face = "italic")
      )
    )
  save_fig(panel_p, "Fig_PANEL_multidb_summary",
           fig_ov, width = 20, height = 16)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 26.  EXCEL EXPORT — ALL RESULTS
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[XL] Exporting Excel tables ────────────────────────────────────────")

## ── Unified clean_result: works for both ORA and GSEA result objects ─────────
## ORA  columns: GeneRatio, BgRatio, Count, geneID
## GSEA columns: setSize, enrichmentScore, NES, core_enrichment
## Both share  : ID, Description, pvalue, p.adjust, qvalue
clean_result <- function(res, db_label) {
  if (is.null(res) || nrow(res) == 0) return(NULL)

  df <- res@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(Database = db_label)

  if (nrow(df) == 0) return(NULL)

  ## Detect result type by presence of NES column
  is_gsea <- "NES" %in% names(df)

  base_cols  <- c("Database", "ID", "Description",
                  "pvalue", "p.adjust", "qvalue")
  ora_extras  <- c("GeneRatio", "BgRatio", "Count", "geneID")
  gsea_extras <- c("setSize", "enrichmentScore", "NES",
                   "rank", "leading_edge", "core_enrichment")

  extra_cols <- if (is_gsea) gsea_extras else ora_extras
  keep_cols  <- c(base_cols, intersect(extra_cols, names(df)))

  df %>%
    dplyr::select(all_of(keep_cols)) %>%
    arrange(p.adjust)
}

## ── Gene signature table: guard optional columns ──────────────────────────────
sig_cols_available <- intersect(
  c("gene", "z_meta", "z_std", "bag_freq", "bag_frac",
    "direction", "bag_tier", "n_cohorts_used",
    "ENTREZID", "ENSEMBL", "UNIPROT", "GENENAME"),
  names(gene_df)
)

## Build named list of all result frames
all_tables <- list(
  "Gene_Signature"     = gene_df %>%
    dplyr::select(all_of(sig_cols_available)) %>%
    arrange(desc(abs(z_meta))),
  "GO_BP"              = clean_result(go_bp,        "GO-BP"),
  "GO_MF"              = clean_result(go_mf,        "GO-MF"),
  "GO_CC"              = clean_result(go_cc,        "GO-CC"),
  "KEGG"               = clean_result(kegg_res,     "KEGG"),
  "Reactome"           = clean_result(reactome_res, "Reactome"),
  "Hallmark_ORA"       = clean_result(hallmark_res, "Hallmark"),
  "C7_Immunologic"     = clean_result(c7_res,       "C7-ImmuneSigDB"),
  "C8_CellType"        = clean_result(c8_res,       "C8-CellType"),
  "Disease_Ontology"   = clean_result(do_res,       "DisOntology"),
  "NCG"                = clean_result(ncg_res,      "NCG"),
  "DisGeNET"           = clean_result(dgn_res,      "DisGeNET"),
  "GSEA_Hallmark"      = clean_result(gsea_hallmark,"GSEA-Hallmark"),
  "GSEA_GO_BP"         = clean_result(gsea_go_bp,  "GSEA-GO-BP"),
  "GSEA_KEGG"          = clean_result(gsea_kegg,   "GSEA-KEGG"),
  "GSEA_C2"            = clean_result(gsea_c2,     "GSEA-C2"),
  "GO_BP_RiskGenes"    = clean_result(go_bp_pos,   "GO-BP (Risk-promoting)"),
  "GO_BP_Protective"   = clean_result(go_bp_neg,   "GO-BP (Protective)"),
  "KEGG_RiskGenes"     = clean_result(kegg_pos,    "KEGG (Risk-promoting)"),
  "KEGG_Protective"    = clean_result(kegg_neg,    "KEGG (Protective)"),
  "GO_BP_CoreGenes"    = clean_result(go_bp_core,  "GO-BP (Core 8/8)")
)

## Add C2 subcollection results
for (sc in names(c2_results)) {
  nm <- paste0("C2_", str_replace_all(sc, "[^A-Za-z0-9]", "_"))
  all_tables[[nm]] <- clean_result(c2_results[[sc]], sc)
}

all_tables <- Filter(Negate(is.null), all_tables)

## Build workbook
wb <- createWorkbook()

hdr_style <- createStyle(
  fontName = "Calibri", fontSize = 10,
  fontColour = "white", fgFill = "#1F3864",
  halign = "center", valign = "center",
  textDecoration = "bold",
  border = "Bottom", borderColour = "white",
  wrapText = TRUE
)
body_style <- createStyle(
  fontName = "Calibri", fontSize = 9,
  halign = "left", valign = "top",
  border = "Bottom", borderColour = "#D9D9D9"
)
alt_style <- createStyle(
  fontName = "Calibri", fontSize = 9,
  fgFill = "#F2F7FF", halign = "left"
)
pval_style <- createStyle(numFmt = "0.00E+00")
num_style  <- createStyle(numFmt = "0.000")

for (sn in names(all_tables)) {
  df <- all_tables[[sn]]
  if (is.null(df) || nrow(df) == 0) next

  addWorksheet(wb, sn)
  writeData(wb, sn, df, headerStyle = hdr_style)
  nr <- nrow(df)
  nc <- ncol(df)

  addStyle(wb, sn, body_style, rows = 2:(nr + 1), cols = 1:nc,
           gridExpand = TRUE)
  ## Alternating row shading
  if (nr >= 2) {
    alt_rows <- seq(3, nr + 1, by = 2)
    addStyle(wb, sn, alt_style, rows = alt_rows, cols = 1:nc,
             gridExpand = TRUE, stack = TRUE)
  }
  ## p-value formatting
  p_cols <- which(names(df) %in% c("pvalue","p.adjust","qvalue"))
  if (length(p_cols))
    addStyle(wb, sn, pval_style, rows = 2:(nr + 1), cols = p_cols,
             gridExpand = TRUE, stack = TRUE)

  setColWidths(wb, sn, cols = 1:nc, widths = "auto")
  setRowHeights(wb, sn, rows = 1, heights = 30)
  freezePane(wb, sn, firstRow = TRUE)
}

## Summary sheet
summary_tbl <- tibble(
  Sheet    = names(all_tables),
  N_terms  = sapply(all_tables, nrow),
  Database = str_replace_all(names(all_tables), "_", " ")
)
addWorksheet(wb, "SUMMARY", tabColour = "#C0392B")
writeData(wb, "SUMMARY", summary_tbl, headerStyle = hdr_style)
setColWidths(wb, "SUMMARY", cols = 1:3, widths = c(25, 12, 35))
freezePane(wb, "SUMMARY", firstRow = TRUE)

xl_path <- file.path(OUT_DIR, "tables",
                      "Pathway_Enrichment_ALL_Results_Q1_v2.xlsx")
saveWorkbook(wb, xl_path, overwrite = TRUE)
message("  Excel saved: ", basename(xl_path))

# ═══════════════════════════════════════════════════════════════════════════════
# 27.  SAVE R OBJECTS
# ═══════════════════════════════════════════════════════════════════════════════
message("\n[RD] Saving RData ──────────────────────────────────────────────────")

save(
  ## Gene info
  gene_df, gene_entrez, gsea_vec,
  entrez_pos, entrez_neg, entrez_core,
  ## ORA results
  go_bp, go_mf, go_cc,
  kegg_res, reactome_res,
  hallmark_res, c7_res, c8_res,
  do_res, ncg_res, dgn_res,
  c2_results,
  ## Direction-stratified
  go_bp_pos, go_bp_neg, go_mf_pos, go_mf_neg,
  kegg_pos, kegg_neg, go_bp_core,
  ## GSEA results
  gsea_hallmark, gsea_go_bp, gsea_kegg, gsea_c2,
  ## compareCluster
  cc_go, cc_kegg, cc_hall,
  file = file.path(OUT_DIR, "rdata",
                   "pathway_enrichment_ALL_results.RData")
)
message("  RData saved.")

# ═══════════════════════════════════════════════════════════════════════════════
# 28.  SESSION INFO
# ═══════════════════════════════════════════════════════════════════════════════
sink(file.path(OUT_DIR, "session_info.txt"))
cat("Pathway Enrichment Analysis — Q1 v2\n")
cat("Date:", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()

# ───────────────────────────────────────────────────────────────────────────────
message("\n", strrep("═", 65))
message("  PIPELINE COMPLETE")
message(strrep("═", 65))
message("  Root   : ", OUT_DIR)
message("  Figures: ", OUT_DIR, "/figures/")
message("    ├── ORA/        — dotplots, barplots, emaplots")
message("    ├── GSEA/       — dotplots, ridgeplots, running-score PDFs")
message("    ├── comparative/— direction-stratified, compareCluster")
message("    ├── network/    — cnetplots, emaplots, UpSet")
message("    └── overview/   — lollipop, heatmap, panel")
message("  Tables : ", xl_path)
message("  RData  : ", OUT_DIR, "/rdata/pathway_enrichment_ALL_results.RData")
message(strrep("═", 65))
