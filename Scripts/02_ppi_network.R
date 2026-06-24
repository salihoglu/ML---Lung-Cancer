################################################################################
##  PROTEIN-PROTEIN INTERACTION NETWORK ANALYSIS                              ##
##  Ferroptosis-Related Prognostic Genes | LUAD                               ##
##  PPI figure generation with tidygraph + ggraph                             ##
##  Rana Salihoglu                                                            ##
################################################################################

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
bioc_pkgs <- c(
  "org.Hs.eg.db", "AnnotationDbi", "STRINGdb",
  "ComplexHeatmap", "circlize"
)

cran_pkgs <- c(
  "igraph", "tidygraph", "ggraph", "ggplot2", "ggrepel", "patchwork",
  "dplyr", "tidyr", "tibble", "stringr", "purrr", "forcats",
  "openxlsx", "jsonlite", "httr",
  "RColorBrewer", "viridis", "scales", "pheatmap", "corrplot"
)

install_if_missing <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (bioc) {
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else {
        install.packages(p, repos = "https://cloud.r-project.org")
      }
    }
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

install_if_missing(bioc_pkgs, bioc = TRUE)
install_if_missing(cran_pkgs, bioc = FALSE)

suppressPackageStartupMessages({
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(forcats)
  library(openxlsx)
  library(jsonlite)
  library(httr)
  library(RColorBrewer)
  library(viridis)
  library(scales)
  library(pheatmap)
  library(corrplot)
  library(ComplexHeatmap)
  library(circlize)
  library(STRINGdb)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

set.seed(2025)

# explicit igraph bindings
.igraph_env <- list(
  simplify              = igraph::simplify,
  components            = igraph::components,
  edge_density          = igraph::edge_density,
  is_connected          = igraph::is_connected,
  mean_distance         = igraph::mean_distance,
  diameter              = igraph::diameter,
  transitivity          = igraph::transitivity,
  assortativity_degree  = igraph::assortativity_degree,
  degree                = igraph::degree,
  betweenness           = igraph::betweenness,
  closeness             = igraph::closeness,
  eigen_centrality      = igraph::eigen_centrality,
  page_rank             = igraph::page_rank,
  strength              = igraph::strength,
  neighbors             = igraph::neighbors,
  induced_subgraph      = igraph::induced_subgraph,
  delete_vertices       = igraph::delete_vertices,
  max_cliques           = igraph::max_cliques,
  fit_power_law         = igraph::fit_power_law,
  graph_from_data_frame = igraph::graph_from_data_frame,
  as_data_frame         = igraph::as_data_frame,
  vcount                = igraph::vcount,
  ecount                = igraph::ecount,
  write_graph           = igraph::write_graph
)
list2env(.igraph_env, envir = globalenv())

# ── 1. PATHS & DIRECTORIES ────────────────────────────────────────────────────
BASE_DIR   <- "/test"
SCRIPT_DIR <- file.path(BASE_DIR, "bioinf_scripts/scripts")
OUT_DIR    <- file.path(SCRIPT_DIR, "ppi_network")

for (d in c(
  "figures/network",
  "figures/topology",
  "figures/hubs",
  "figures/subnetwork",
  "figures/robustness",
  "tables",
  "graph",
  "rdata"
)) {
  dir.create(file.path(OUT_DIR, d), showWarnings = FALSE, recursive = TRUE)
}

fig_net  <- file.path(OUT_DIR, "figures/network")
fig_topo <- file.path(OUT_DIR, "figures/topology")
fig_hub  <- file.path(OUT_DIR, "figures/hubs")
fig_sub  <- file.path(OUT_DIR, "figures/subnetwork")
fig_rob  <- file.path(OUT_DIR, "figures/robustness")

STRING_VERSION      <- "12.0"
STRING_SCORE_CUTOFF <- 400L
SPECIES_ID          <- 9606L

message("Output root : ", OUT_DIR)

# ── 2. THEME & SAVE HELPERS ───────────────────────────────────────────────────
q1_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2,
                                hjust = 0, color = "#1A1A1A"),
      plot.subtitle = element_text(size = base_size - 0.5, color = "#4D4D4D",
                                   hjust = 0, margin = margin(b = 6)),
      plot.caption = element_text(size = 7.5, color = "#888888",
                                  hjust = 1, face = "italic"),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 1, color = "#2B2B2B"),
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.35),
      panel.grid.major = element_line(linetype = "dashed",
                                      color = "#E5E5E5", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      legend.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(12, 18, 12, 12),
      strip.background = element_rect(fill = "#EFEFEF", color = "grey60"),
      strip.text = element_text(face = "bold", size = base_size - 1)
    )
}

save_fig <- function(p, name, subdir = fig_net, width = 10, height = 8, dpi = 300) {
  if (is.null(p)) return(invisible(NULL))
  tryCatch({
    suppressMessages({
      ggplot2::ggsave(
        file.path(subdir, paste0(name, ".pdf")),
        plot = p, width = width, height = height, device = cairo_pdf, bg = "white"
      )
      ggplot2::ggsave(
        file.path(subdir, paste0(name, ".png")),
        plot = p, width = width, height = height, dpi = dpi, bg = "white"
      )
      ggplot2::ggsave(
        file.path(subdir, paste0(name, ".tiff")),
        plot = p, width = width, height = height,
        dpi = dpi, compression = "lzw", bg = "white"
      )
    })
    message("  [FIG] ", name)
  }, error = function(e) {
    message("  [FIG SKIP] ", name, " — ", e$message)
  })
  invisible(p)
}

rescale_safe <- function(x, to = c(20, 80)) {
  if (length(x) == 0) return(numeric(0))
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(mean(to), length(x)))
  }
  scales::rescale(x, to = to, from = rng)
}

make_comm_palette <- function(n) {
  if (n <= 0) return(character(0))
  base_cols <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2",
    "#B07AA1", "#EDC948", "#9C755F", "#BAB0AC", "#2F4B7C",
    "#D45087", "#00A6ED", "#7F3C8D", "#11A579", "#3969AC"
  )
  if (n <= length(base_cols)) {
    return(setNames(base_cols[seq_len(n)], as.character(seq_len(n))))
  }
  setNames(colorRampPalette(base_cols)(n), as.character(seq_len(n)))
}

# ── 2b. PPI PLOTTING HELPERS ───────────────────────────────────────────
make_tbl_graph_from_lcc <- function(g_lcc, node_tbl) {
  edge_df <- igraph::as_data_frame(g_lcc, what = "edges") %>%
    transmute(from = from, to = to, weight = weight)
  
  node_df <- data.frame(name = igraph::V(g_lcc)$name, stringsAsFactors = FALSE) %>%
    left_join(
      node_tbl %>%
        select(
          gene, degree, weighted_degree, betweenness, closeness,
          eigenvector, pagerank, mcc, hub_score, is_hub,
          z_meta, abs_z, direction, bag_tier, bag_frac, community
        ),
      by = c("name" = "gene")
    ) %>%
    mutate(
      degree      = replace_na(degree, 0),
      hub_score   = replace_na(hub_score, 0),
      z_meta      = replace_na(z_meta, 0),
      abs_z       = replace_na(abs_z, 0),
      direction   = replace_na(direction, "Unknown"),
      bag_tier    = replace_na(bag_tier, "Unknown"),
      community   = replace_na(community, 0L),
      is_hub      = replace_na(is_hub, FALSE),
      label_main  = ifelse(is_hub | degree >= quantile(degree, 0.90, na.rm = TRUE), name, ""),
      label_hub   = ifelse(is_hub, name, "")
    )
  
  tidygraph::tbl_graph(
    nodes = node_df,
    edges = edge_df,
    directed = FALSE
  )
}

layout_tbl_graph <- function(tbl_g, layout = "stress") {
  ggraph::create_layout(tbl_g, layout = layout)
}

zoom_layout <- function(layout_obj, zoom = 1.20) {

  x0 <- mean(range(layout_obj$x, na.rm = TRUE))
  y0 <- mean(range(layout_obj$y, na.rm = TRUE))
  layout_obj$x <- x0 + (layout_obj$x - x0) / zoom
  layout_obj$y <- y0 + (layout_obj$y - y0) / zoom
  layout_obj
}

shape_legend_guide <- guide_legend(
  order = 1,
  override.aes = list(
    size = 6,
    fill = "grey75",
    colour = "#333333",
    stroke = 0.7,
    alpha = 1
  )
)

size_legend_guide <- guide_legend(
  order = 2,
  override.aes = list(
    shape = 21,
    fill = "grey80",
    colour = "#333333",
    alpha = 1,
    stroke = 0.5
  )
)

fill_colorbar_guide <- guide_colorbar(
  order = 3,
  barheight = unit(45, "mm"),
  frame.colour = "grey40",
  ticks.colour = "grey30"
)

fill_discrete_guide <- guide_legend(
  order = 3,
  override.aes = list(
    shape = 21,
    size = 6,
    colour = "#333333",
    stroke = 0.4,
    alpha = 1
  )
)

plot_ppi_main <- function(layout_obj, score_cutoff, string_version) {
  ggraph(layout_obj) +
    geom_edge_link(
      aes(width = weight, alpha = weight),
      colour = "#B8C0CC", show.legend = FALSE, lineend = "round"
    ) +
    scale_edge_width(range = c(0.15, 1.2)) +
    scale_edge_alpha(range = c(0.05, 0.45)) +
    geom_node_point(
      aes(size = degree, fill = z_meta, shape = is_hub),
      colour = "white", stroke = 0.55, alpha = 0.97
    ) +
    scale_shape_manual(
      values = c("FALSE" = 21, "TRUE" = 23),
      labels = c("Other", "Hub"),
      name = "Node type",
      guide = shape_legend_guide
    ) +
    scale_size_continuous(
      range = c(3.4, 12.5),
      name = "Degree",
      guide = size_legend_guide
    ) +
    scale_fill_gradient2(
      low = "#2C7FB8",
      mid = "grey95",
      high = "#D7301F",
      midpoint = 0,
      name = "Meta-Z",
      guide = fill_colorbar_guide
    ) +
    geom_node_text(
      aes(label = label_hub),
      repel = TRUE,
      size = 3.2,
      family = "sans",
      fontface = "bold",
      colour = "#1A1A1A",
      bg.colour = "white",
      bg.r = 0.09,
      max.overlaps = Inf,
      seed = 2025
    ) +
    labs(
      title = "Ferroptosis-Related Protein-Protein Interaction Network",
      subtitle = paste0(
        "STRING v", string_version,
        " | score ≥ ", score_cutoff,
        " | publication-grade layout (tidygraph + ggraph)"
      ),
      caption = "Node size = degree; node fill = meta-Z; diamond = hub"
    ) +
    theme_graph(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, colour = "#555555"),
      legend.position = "right",
      legend.box = "vertical"
    )
}

plot_ppi_modules <- function(layout_obj, modularity_q) {
  n_comm <- max(layout_obj$community, na.rm = TRUE)
  pal_comm <- make_comm_palette(n_comm)
  
  layout_obj$community_f <- factor(layout_obj$community)
  
  ggraph(layout_obj) +
    geom_edge_link(
      aes(alpha = weight),
      colour = "#CDD3DB", show.legend = FALSE
    ) +
    scale_edge_alpha(range = c(0.04, 0.35)) +
    geom_node_point(
      aes(size = degree, fill = community_f, shape = is_hub),
      colour = "white", stroke = 0.5, alpha = 0.96
    ) +
    scale_shape_manual(
      values = c("FALSE" = 21, "TRUE" = 23),
      labels = c("Other", "Hub"),
      name = "Node type",
      guide = shape_legend_guide
    ) +
    scale_fill_manual(
      values = pal_comm,
      name = "Module",
      guide = fill_discrete_guide
    ) +
    scale_size_continuous(
      range = c(3.0, 10.5),
      name = "Degree",
      guide = size_legend_guide
    ) +
    geom_node_text(
      aes(label = label_hub),
      repel = TRUE,
      size = 3.0,
      family = "sans",
      fontface = "bold",
      colour = "#1A1A1A",
      bg.colour = "white",
      bg.r = 0.09,
      max.overlaps = Inf,
      seed = 2025
    ) +
    labs(
      title = "PPI Network — Community Structure",
      subtitle = paste0("Louvain modules | modularity Q = ", round(modularity_q, 3)),
      caption = "Node colour = module; diamond = hub"
    ) +
    theme_graph(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "#555555"),
      legend.position = "right",
      legend.box = "vertical"
    )
}

plot_ppi_hub_context <- function(layout_obj) {
  ggraph(layout_obj) +
    geom_edge_link(
      aes(width = weight, alpha = weight),
      colour = "#9EA7B3", show.legend = FALSE, lineend = "round"
    ) +
    scale_edge_width(range = c(0.2, 1.8)) +
    scale_edge_alpha(range = c(0.08, 0.55)) +
    geom_node_point(
      aes(size = hub_score, fill = z_meta, shape = is_hub),
      colour = "white", stroke = 0.6, alpha = 0.97
    ) +
    scale_shape_manual(
      values = c("FALSE" = 21, "TRUE" = 23),
      labels = c("Neighbour", "Hub"),
      name = "Node type",
      guide = shape_legend_guide
    ) +
    scale_size_continuous(
      range = c(4.5, 16.5),
      name = "Hub score",
      guide = size_legend_guide
    ) +
    scale_fill_gradient2(
      low = "#2C7FB8",
      mid = "grey95",
      high = "#D7301F",
      midpoint = 0,
      name = "Meta-Z",
      guide = fill_colorbar_guide
    ) +
    geom_node_text(
      aes(label = label_main),
      repel = TRUE,
      size = 3.2,
      family = "sans",
      fontface = "bold",
      colour = "#1A1A1A",
      bg.colour = "white",
      bg.r = 0.09,
      max.overlaps = Inf,
      seed = 2025
    ) +
    labs(
      title = "Hub-Centric PPI Subnetwork",
      subtitle = "Hub genes and first-degree neighbors",
      caption = "Node size = composite hub score; fill = meta-Z"
    ) +
    theme_graph(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "#555555"),
      legend.position = "right",
      legend.box = "vertical"
    )
}

# ── 3. LOAD & ANNOTATE GENE SIGNATURE ────────────────────────────────────────
message("\n[A] Loading gene signature")

gene_csv <- Filter(file.exists, c(
  file.path(SCRIPT_DIR, "possible_prognostic_genes_FULLTRAIN.csv"),
  file.path(BASE_DIR,   "possible_prognostic_genes_FULLTRAIN.csv"),
  list.files(BASE_DIR, "possible_prognostic_genes_FULLTRAIN.csv",
             recursive = TRUE, full.names = TRUE)
))[1]

if (is.na(gene_csv)) stop("Gene CSV not found.")

gene_df <- read.csv(gene_csv, stringsAsFactors = FALSE) %>%
  mutate(
    direction = ifelse(z_meta > 0, "Risk-promoting", "Protective"),
    bag_tier = case_when(
      bag_frac == 1.00 ~ "Core (8/8)",
      bag_frac >= 0.75 ~ "Stable (6/8)",
      TRUE             ~ "Variable"
    ),
    abs_z = abs(z_meta)
  ) %>%
  arrange(desc(abs_z))

if (!"gene" %in% names(gene_df) && "feature" %in% names(gene_df)) {
  gene_df$gene <- sub("^RNA__", "", gene_df$feature)
}

stopifnot(
  "'gene' column required" = "gene" %in% names(gene_df),
  "no NA in gene" = !any(is.na(gene_df$gene))
)

gene_symbols <- unique(gene_df$gene)
message(
  "  Genes : ", length(gene_symbols),
  "  (risk=", sum(gene_df$direction == "Risk-promoting"),
  ", protective=", sum(gene_df$direction == "Protective"), ")"
)

id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_symbols,
  columns = c("ENTREZID", "GENENAME", "UNIPROT"),
  keytype = "SYMBOL"
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_df <- merge(
  gene_df, id_map,
  by.x = "gene", by.y = "SYMBOL",
  all.x = TRUE, sort = FALSE
)

stopifnot("'gene' must survive merge" = "gene" %in% names(gene_df))

message("  Mapped to Entrez: ", sum(!is.na(gene_df$ENTREZID)), " / ", nrow(gene_df))

# ── 4. STRING INTERACTIONS ───────────────────────────────────────────────────
message("\n[B] STRING interactions")

string_rest_api <- function(genes, species = 9606L, score_cut = 400L) {
  url <- "https://string-db.org/api/tsv/network"
  body <- list(
    identifiers = paste(genes, collapse = "%0d"),
    species = as.character(species),
    caller_identity = "PROJECT16_LUAD_PPI"
  )
  
  resp <- tryCatch(
    httr::POST(url, body = body, encode = "form", httr::timeout(60)),
    error = function(e) NULL
  )
  if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)
  
  txt <- httr::content(resp, "text", encoding = "UTF-8")
  df <- tryCatch(
    read.table(text = txt, sep = "\t", header = TRUE,
               stringsAsFactors = FALSE, quote = ""),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  if ("preferredName_A" %in% names(df)) {
    df <- df %>%
      rename(
        gene_A = preferredName_A,
        gene_B = preferredName_B,
        combined_score = score
      ) %>%
      mutate(combined_score = as.numeric(combined_score) * 1000)
  } else if (!"gene_A" %in% names(df)) {
    df <- df %>%
      rename_with(~"gene_A", matches("^stringId_A|^protein1")) %>%
      rename_with(~"gene_B", matches("^stringId_B|^protein2"))
  }
  
  df %>%
    filter(
      combined_score >= score_cut,
      gene_A != gene_B,
      gene_A %in% genes,
      gene_B %in% genes
    ) %>%
    mutate(weight = combined_score / 1000) %>%
    select(gene_A, gene_B, combined_score, weight)
}

interactions <- NULL
string_db <- NULL

string_db <- tryCatch({
  db <- STRINGdb$new(
    version = STRING_VERSION,
    species = SPECIES_ID,
    score_threshold = STRING_SCORE_CUTOFF,
    network_type = "full",
    input_directory = file.path(OUT_DIR, "rdata")
  )
  message("  STRINGdb initialised (v", STRING_VERSION, ")")
  db
}, error = function(e) {
  message("  v", STRING_VERSION, " failed: ", e$message)
  tryCatch({
    db <- STRINGdb$new(
      version = "11.5",
      species = SPECIES_ID,
      score_threshold = STRING_SCORE_CUTOFF,
      network_type = "full",
      input_directory = file.path(OUT_DIR, "rdata")
    )
    message("  STRINGdb initialised ")
    db
  }, error = function(e2) {
    message("  STRINGdb unavailable: ", e2$message)
    NULL
  })
})

if (!is.null(string_db)) {
  mp <- tryCatch(
    string_db$map(gene_df, "gene", removeUnmappedRows = TRUE),
    error = function(e) {
      message("  map() failed: ", e$message)
      NULL
    }
  )
  
  if (!is.null(mp) && nrow(mp) > 0) {
    id2sym <- setNames(mp$gene, mp$STRING_id)
    raw <- tryCatch(
      string_db$get_interactions(mp$STRING_id),
      error = function(e) {
        message("  get_interactions() failed: ", e$message)
        NULL
      }
    )
    
    if (!is.null(raw) && nrow(raw) > 0) {
      interactions <- raw %>%
        filter(combined_score >= STRING_SCORE_CUTOFF) %>%
        mutate(
          gene_A = id2sym[from],
          gene_B = id2sym[to],
          weight = combined_score / 1000
        ) %>%
        filter(!is.na(gene_A), !is.na(gene_B), gene_A != gene_B) %>%
        select(gene_A, gene_B, combined_score, weight)
      message("  Package interactions: ", nrow(interactions))
    }
  }
}

if (is.null(interactions) || nrow(interactions) == 0) {
  message("  Trying REST API fallback …")
  interactions <- string_rest_api(gene_symbols, SPECIES_ID, STRING_SCORE_CUTOFF)
  if (!is.null(interactions) && nrow(interactions) > 0) {
    message("  REST API interactions: ", nrow(interactions))
  }
}

if (is.null(interactions) || nrow(interactions) == 0) {
  stop("No interactions retrieved. Check internet / lower STRING_SCORE_CUTOFF.")
}

write.csv(
  interactions,
  file.path(OUT_DIR, "tables", "STRING_interactions.csv"),
  row.names = FALSE
)

message("  Unique genes in edges: ",
        length(unique(c(interactions$gene_A, interactions$gene_B))))

# ── 5. BUILD IGRAPH ──────────────────────────────────────────────────────────
message("\n[C] Building igraph")

net_genes <- unique(c(interactions$gene_A, interactions$gene_B))

vtx_df <- data.frame(
  name      = net_genes,
  z_meta    = gene_df$z_meta[match(net_genes, gene_df$gene)],
  abs_z     = gene_df$abs_z[match(net_genes, gene_df$gene)],
  direction = gene_df$direction[match(net_genes, gene_df$gene)],
  bag_tier  = gene_df$bag_tier[match(net_genes, gene_df$gene)],
  bag_frac  = gene_df$bag_frac[match(net_genes, gene_df$gene)],
  stringsAsFactors = FALSE
)

vtx_df$z_meta    <- ifelse(is.na(vtx_df$z_meta), 0, vtx_df$z_meta)
vtx_df$abs_z     <- ifelse(is.na(vtx_df$abs_z), 0, vtx_df$abs_z)
vtx_df$direction <- ifelse(is.na(vtx_df$direction), "Unknown", vtx_df$direction)
vtx_df$bag_tier  <- ifelse(is.na(vtx_df$bag_tier), "Unknown", vtx_df$bag_tier)
vtx_df           <- vtx_df[!duplicated(vtx_df$name), ]

g <- igraph::graph_from_data_frame(
  d = as.data.frame(interactions[, c("gene_A", "gene_B", "weight")]),
  directed = FALSE,
  vertices = vtx_df
)

g <- igraph::simplify(
  g,
  remove.multiple = TRUE,
  remove.loops = TRUE,
  edge.attr.comb = list(weight = "max")
)

message("  Full network — nodes: ", igraph::vcount(g),
        " | edges: ", igraph::ecount(g))

comp <- igraph::components(g)
lcc_idx <- which(comp$membership == which.max(comp$csize))
g_lcc <- igraph::induced_subgraph(g, lcc_idx)

message("  LCC — nodes: ", igraph::vcount(g_lcc),
        " | edges: ", igraph::ecount(g_lcc))

global_metrics <- data.frame(
  Metric = c(
    "Nodes (LCC)", "Edges (LCC)", "Density",
    "Avg degree", "Avg path length", "Diameter",
    "Clustering coeff", "N components", "Largest comp size",
    "Assortativity"
  ),
  Value = c(
    igraph::vcount(g_lcc),
    igraph::ecount(g_lcc),
    round(igraph::edge_density(g_lcc), 5),
    round(mean(igraph::degree(g_lcc)), 3),
    round(igraph::mean_distance(g_lcc, unconnected = TRUE), 3),
    igraph::diameter(g_lcc, unconnected = TRUE),
    round(igraph::transitivity(g_lcc, type = "global"), 4),
    igraph::components(g)$no,
    max(igraph::components(g)$csize),
    round(igraph::assortativity_degree(g_lcc), 4)
  ),
  stringsAsFactors = FALSE
)

print(global_metrics)
write.csv(
  global_metrics,
  file.path(OUT_DIR, "tables", "network_global_metrics.csv"),
  row.names = FALSE
)

# ── 6. CENTRALITY & HUB IDENTIFICATION ──────────────────────────────────────
message("\n[D] Centrality metrics")

deg    <- igraph::degree(g_lcc, normalized = FALSE)
str_wt <- igraph::strength(g_lcc, weights = igraph::E(g_lcc)$weight)
bet    <- igraph::betweenness(g_lcc, normalized = TRUE, weights = igraph::E(g_lcc)$weight)
clos   <- igraph::closeness(g_lcc, normalized = TRUE, weights = igraph::E(g_lcc)$weight)
eig    <- igraph::eigen_centrality(g_lcc, weights = igraph::E(g_lcc)$weight)$vector
pr     <- igraph::page_rank(g_lcc, weights = igraph::E(g_lcc)$weight)$vector

node_names <- igraph::V(g_lcc)$name
mcc_score <- setNames(
  sapply(node_names, function(v) {
    nb <- igraph::neighbors(g_lcc, v)
    if (length(nb) == 0) return(0L)
    sg <- igraph::induced_subgraph(g_lcc, c(v, nb$name))
    cl <- igraph::max_cliques(sg, min = 2)
    sum(sapply(cl, length))
  }),
  node_names
)

node_tbl <- data.frame(
  gene            = node_names,
  degree          = as.integer(deg[node_names]),
  weighted_degree = as.numeric(str_wt[node_names]),
  betweenness     = as.numeric(bet[node_names]),
  closeness       = as.numeric(clos[node_names]),
  eigenvector     = as.numeric(eig[node_names]),
  pagerank        = as.numeric(pr[node_names]),
  mcc             = as.integer(mcc_score[node_names]),
  stringsAsFactors = FALSE
)

node_tbl <- merge(
  node_tbl,
  gene_df[, c("gene", "z_meta", "abs_z", "direction", "bag_tier", "bag_frac")],
  by = "gene", all.x = TRUE, sort = FALSE
)

node_tbl$z_meta    <- ifelse(is.na(node_tbl$z_meta), 0, node_tbl$z_meta)
node_tbl$direction <- ifelse(is.na(node_tbl$direction), "Unknown", node_tbl$direction)
node_tbl$bag_tier  <- ifelse(is.na(node_tbl$bag_tier), "Unknown", node_tbl$bag_tier)

n <- nrow(node_tbl)

node_tbl <- node_tbl %>%
  mutate(
    r_deg   = rank(degree) / n,
    r_bet   = rank(betweenness) / n,
    r_clos  = rank(closeness) / n,
    r_eig   = rank(eigenvector) / n,
    r_pr    = rank(pagerank) / n,
    r_mcc   = rank(mcc) / n,
    hub_score = (r_deg + r_bet + r_clos + r_eig + r_pr + r_mcc) / 6,
    top10_deg  = degree      >= quantile(degree, 0.90),
    top20_bet  = betweenness >= quantile(betweenness, 0.80),
    top20_clos = closeness   >= quantile(closeness, 0.80),
    top20_eig  = eigenvector >= quantile(eigenvector, 0.80),
    top20_pr   = pagerank    >= quantile(pagerank, 0.80),
    top20_mcc  = mcc         >= quantile(mcc, 0.80),
    n_criteria = top10_deg + top20_bet + top20_clos + top20_eig + top20_pr + top20_mcc,
    is_hub = as.logical(n_criteria >= 3)
  ) %>%
  arrange(desc(hub_score))

stopifnot(
  "is_hub must be logical" = is.logical(node_tbl$is_hub),
  "no NA in is_hub" = !any(is.na(node_tbl$is_hub)),
  "gene col in node_tbl" = "gene" %in% names(node_tbl)
)

hub_genes <- node_tbl$gene[node_tbl$is_hub]

message("  Hub genes (≥3 criteria): ", length(hub_genes))
message("  Hubs: ", paste(sort(hub_genes), collapse = ", "))

write.csv(
  node_tbl,
  file.path(OUT_DIR, "tables", "node_centrality.csv"),
  row.names = FALSE
)

# ── 7. COMMUNITY DETECTION ───────────────────────────────────────────────────
message("\n[E] Community detection")

comm_louvain <- igraph::cluster_louvain(g_lcc, weights = igraph::E(g_lcc)$weight)
comm_walk    <- igraph::cluster_walktrap(g_lcc, weights = igraph::E(g_lcc)$weight)

message(
  "  Louvain  — Q=", round(igraph::modularity(comm_louvain), 4),
  " | modules=", length(unique(igraph::membership(comm_louvain)))
)
message(
  "  Walktrap — Q=", round(igraph::modularity(comm_walk), 4),
  " | modules=", length(unique(igraph::membership(comm_walk)))
)

node_tbl$community <- as.integer(igraph::membership(comm_louvain)[node_tbl$gene])

igraph::V(g_lcc)$community <- igraph::membership(comm_louvain)[igraph::V(g_lcc)$name]
igraph::V(g_lcc)$is_hub    <- igraph::V(g_lcc)$name %in% hub_genes
igraph::V(g_lcc)$degree_v  <- igraph::degree(g_lcc)
igraph::V(g_lcc)$z_meta    <- node_tbl$z_meta[match(igraph::V(g_lcc)$name, node_tbl$gene)]
igraph::V(g_lcc)$hub_score <- node_tbl$hub_score[match(igraph::V(g_lcc)$name, node_tbl$gene)]

module_tbl <- node_tbl %>%
  group_by(community) %>%
  summarise(
    n_genes = n(),
    genes = paste(sort(gene), collapse = ", "),
    n_hubs = sum(is_hub),
    hub_genes_m = paste(sort(gene[is_hub]), collapse = ", "),
    mean_degree = round(mean(degree), 2),
    mean_z = round(mean(z_meta, na.rm = TRUE), 3),
    pct_risk = round(mean(direction == "Risk-promoting", na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genes))

write.csv(
  module_tbl,
  file.path(OUT_DIR, "tables", "community_modules.csv"),
  row.names = FALSE
)

# ── 8. NETWORK VISUALISATIONS ────────────────────────────────────────────────
message("\n[F] Network visualisations")

tbl_g_lcc <- make_tbl_graph_from_lcc(g_lcc, node_tbl)

layout_main   <- layout_tbl_graph(tbl_g_lcc, layout = "stress")
layout_fr     <- layout_tbl_graph(tbl_g_lcc, layout = "fr")
layout_circle <- layout_tbl_graph(tbl_g_lcc, layout = "kk")

# zoomed layouts
layout_main_zoom   <- zoom_layout(layout_main, zoom = 1.28)
layout_fr_zoom     <- zoom_layout(layout_fr, zoom = 1.24)
layout_circle_zoom <- zoom_layout(layout_circle, zoom = 1.18)

# F1. Main publication PPI figure
p_full <- plot_ppi_main(
  layout_obj = layout_main_zoom,
  score_cutoff = STRING_SCORE_CUTOFF,
  string_version = STRING_VERSION
)
save_fig(p_full, "Fig_F1_main_PPI_network", fig_net, 14, 12)

# F2. Community network
p_comm <- plot_ppi_modules(
  layout_obj = layout_main_zoom,
  modularity_q = igraph::modularity(comm_louvain)
)
save_fig(p_comm, "Fig_F2_module_colored_PPI", fig_net, 14, 12)

# F3. Alternative clean figure for supplement
p_full_fr <- ggraph(layout_fr_zoom) +
  geom_edge_link(
    aes(width = weight, alpha = weight),
    colour = "#C6CDD6", show.legend = FALSE
  ) +
  scale_edge_width(range = c(0.15, 1.0)) +
  scale_edge_alpha(range = c(0.04, 0.35)) +
  geom_node_point(
    aes(size = degree, fill = bag_tier, shape = is_hub),
    colour = "white", stroke = 0.5, alpha = 0.97
  ) +
  scale_shape_manual(
    values = c("FALSE" = 21, "TRUE" = 23),
    labels = c("Other", "Hub"),
    name = "Node type",
    guide = shape_legend_guide
  ) +
  scale_size_continuous(
    range = c(3.0, 11.0),
    name = "Degree",
    guide = size_legend_guide
  ) +
  scale_fill_manual(
    values = c(
      "Core (8/8)" = "#1A1A2E",
      "Stable (6/8)" = "#5C6BC0",
      "Variable" = "#B0BEC5",
      "Unknown" = "#D0D5DD"
    ),
    name = "Bagging tier",
    guide = fill_discrete_guide
  ) +
  geom_node_text(
    aes(label = label_hub),
    repel = TRUE,
    size = 3.0,
    family = "sans",
    fontface = "bold",
    colour = "#1A1A1A",
    bg.colour = "white",
    bg.r = 0.09,
    max.overlaps = Inf,
    seed = 2025
  ) +
  labs(
    title = "PPI Network — Stability-Oriented View",
    subtitle = "Node fill = bagging stability tier",
    caption = "Supplementary figure"
  ) +
  theme_graph(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "#555555"),
    legend.position = "right",
    legend.box = "vertical"
  )
save_fig(p_full_fr, "Fig_F3_PPI_stability_view", fig_net, 14, 12)

# F4. Hub-context subnetwork
hub_nbr <- unique(c(
  hub_genes,
  unlist(lapply(
    hub_genes[hub_genes %in% igraph::V(g_lcc)$name],
    function(h) igraph::neighbors(g_lcc, h)$name
  ))
))

g_hub <- igraph::induced_subgraph(g_lcc, hub_nbr)
tbl_g_hub <- make_tbl_graph_from_lcc(g_hub, node_tbl)
layout_hub <- layout_tbl_graph(tbl_g_hub, layout = "fr")
layout_hub_zoom <- zoom_layout(layout_hub, zoom = 1.32)

p_hub <- plot_ppi_hub_context(layout_hub_zoom)
save_fig(p_hub, "Fig_F4_hub_context_subnetwork", fig_hub, 13, 11)

# F5. Chord diagram — inter-module interactions
message("  Building chord diagram …")

edge_df <- igraph::as_data_frame(g_lcc, what = "edges")
comm_map <- setNames(paste0("M", node_tbl$community), node_tbl$gene)

chord_raw <- edge_df %>%
  mutate(modA = comm_map[from], modB = comm_map[to]) %>%
  filter(!is.na(modA), !is.na(modB)) %>%
  group_by(modA, modB) %>%
  summarise(n = n(), .groups = "drop")

if (nrow(chord_raw) > 0) {
  mods <- sort(unique(c(chord_raw$modA, chord_raw$modB)))
  cmat <- matrix(0, length(mods), length(mods), dimnames = list(mods, mods))
  
  for (i in seq_len(nrow(chord_raw))) {
    a <- chord_raw$modA[i]
    b <- chord_raw$modB[i]
    v <- chord_raw$n[i]
    cmat[a, b] <- cmat[a, b] + v
    cmat[b, a] <- cmat[b, a] + v
  }
  
  diag(cmat) <- diag(cmat) / 2
  chord_col <- setNames(
    colorRampPalette(brewer.pal(min(length(mods), 8), "Set2"))(length(mods)),
    mods
  )
  
  pdf(file.path(fig_net, "Fig_F5_module_chord_diagram.pdf"), width = 10, height = 10)
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4, clock.wise = FALSE)
  chordDiagram(
    cmat,
    grid.col = chord_col,
    transparency = 0.35,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.08)
  )
  circos.trackPlotRegion(
    track.index = 1, bg.border = NA,
    panel.fun = function(x, y) {
      sn <- get.cell.meta.data("sector.index")
      xl <- get.cell.meta.data("xlim")
      circos.text(
        mean(xl), 1.0, sn,
        facing = "clockwise", niceFacing = TRUE,
        adj = c(0, 0.5), cex = 1, font = 2
      )
    }
  )
  title("Inter-module Interactions", sub = "Chord width = edge count between modules")
  circos.clear()
  dev.off()
  message("  [FIG] Fig_F5_module_chord_diagram.pdf")
}

# F6. Centrality heatmap
message("  Building centrality heatmap …")

top_n <- min(40, nrow(node_tbl))
heat_df <- node_tbl[order(-node_tbl$hub_score), ][seq_len(top_n), ]

heat_mat <- scale(as.matrix(
  heat_df[, c("degree", "betweenness", "closeness", "eigenvector", "pagerank", "mcc")]
))
rownames(heat_mat) <- heat_df$gene
colnames(heat_mat) <- c("Degree", "Betweenness", "Closeness", "Eigenvector", "PageRank", "MCC")

row_ha <- HeatmapAnnotation(
  which = "row",
  Direction = heat_df$direction,
  Hub = ifelse(heat_df$is_hub, "Hub", "Other"),
  MetaZ = anno_barplot(
    heat_df$z_meta,
    gp = gpar(
      fill = ifelse(heat_df$z_meta > 0, "#C0392B", "#2980B9"),
      col = NA
    ),
    border = FALSE
  ),
  col = list(
    Direction = c(
      "Risk-promoting" = "#C0392B",
      "Protective" = "#2980B9",
      "Unknown" = "#95A5A6"
    ),
    Hub = c("Hub" = "#F39C12", "Other" = "#D5DBDB")
  ),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
)

col_fun_h <- colorRamp2(c(-3, 0, 3), c("#2980B9", "white", "#C0392B"))

pdf(file.path(fig_hub, "Fig_F6_centrality_heatmap.pdf"),
    width = 11, height = max(7, top_n * 0.28 + 4))
draw(Heatmap(
  heat_mat,
  name = "Z-score",
  col = col_fun_h,
  right_annotation = row_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8.5, fontface = ifelse(heat_df$is_hub, "bold", "plain")),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  column_title = paste0("Centrality Metrics — Top ", top_n, " Genes"),
  column_title_gp = gpar(fontsize = 11, fontface = "bold")
))
dev.off()
message("  [FIG] Fig_F6_centrality_heatmap.pdf")

# F7. Hub bubble ranking
hub_bub <- node_tbl %>%
  slice_max(hub_score, n = 25) %>%
  mutate(gene = fct_reorder(gene, hub_score))

hub_min <- if (any(hub_bub$is_hub)) min(hub_bub$hub_score[hub_bub$is_hub]) else NA

p_hub_bub <- ggplot(
  hub_bub,
  aes(x = hub_score, y = gene, size = degree, color = z_meta, shape = direction)
) +
  { if (!is.na(hub_min)) geom_vline(xintercept = hub_min, linetype = "dashed",
                                    color = "#E74C3C", linewidth = 0.5) } +
  geom_point(alpha = 0.88) +
  scale_color_gradient2(low = "#2980B9", mid = "grey80", high = "#C0392B",
                        midpoint = 0, name = "Meta-Z") +
  scale_size_continuous(range = c(3, 12), name = "Degree") +
  scale_shape_manual(values = c("Risk-promoting" = 17, "Protective" = 16, "Unknown" = 15),
                     name = "Direction") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(
    x = "Composite Hub Score",
    y = NULL,
    title = "Hub Gene Ranking",
    subtitle = paste0("Top 25 by composite centrality | ", length(hub_genes), " hubs identified"),
    caption = "Dashed = hub threshold (≥3 criteria)"
  ) +
  q1_theme(10) +
  theme(
    axis.text.y = element_text(
      size = 8.5,
      face = ifelse(levels(hub_bub$gene) %in% hub_genes, "bold", "plain")
    )
  )

save_fig(p_hub_bub, "Fig_F7_hub_bubble", fig_hub, 11, 10)

# F8. Degree distribution
deg_df <- data.frame(k = igraph::degree(g_lcc)) %>%
  count(k) %>%
  mutate(pk = n / sum(n), log_k = log10(k + 1), log_pk = log10(pk))

pl <- tryCatch(igraph::fit_power_law(igraph::degree(g_lcc) + 1), error = function(e) NULL)

pl_lbl <- ""
if (!is.null(pl)) {
  alpha_val <- tryCatch(round(as.numeric(pl$alpha), 2), error = function(e) NA)
  pval <- NA
  for (pfield in c("KS.p", "p.value", "pvalue", "ks.p")) {
    v <- pl[[pfield]]
    if (!is.null(v) && is.numeric(v)) {
      pval <- round(v, 3)
      break
    }
  }
  pl_lbl <- paste0(
    if (!is.na(alpha_val)) paste0("α = ", alpha_val) else "",
    if (!is.na(pval)) paste0("\np = ", pval) else ""
  )
}

p_deg <- ggplot(deg_df, aes(x = log_k, y = log_pk)) +
  geom_point(color = "#2980B9", size = 2.5, alpha = 0.80) +
  geom_smooth(method = "lm", se = TRUE, color = "#C0392B",
              linewidth = 0.7, linetype = "dashed", alpha = 0.12) +
  annotate(
    "text",
    x = max(deg_df$log_k) * 0.6,
    y = max(deg_df$log_pk) * 0.9,
    label = pl_lbl,
    hjust = 0,
    size = 3.5
  ) +
  labs(
    x = "Degree k (log10 scale)",
    y = "P(k) (log10)",
    title = "Degree Distribution — Log-Log",
    subtitle = paste0(igraph::vcount(g_lcc), " nodes | ", igraph::ecount(g_lcc), " edges"),
    caption = "Linear fit tests scale-free topology"
  ) +
  q1_theme(11)

save_fig(p_deg, "Fig_F8_degree_distribution", fig_topo, 9, 7)

# F9. Centrality correlation
cent_cols <- c("degree", "betweenness", "closeness", "eigenvector", "pagerank", "mcc", "hub_score")
corr_mat <- cor(node_tbl[, cent_cols], method = "spearman", use = "complete.obs")
colnames(corr_mat) <- rownames(corr_mat) <- c(
  "Degree", "Betweenness", "Closeness", "Eigenvector", "PageRank", "MCC", "HubScore"
)

pdf(file.path(fig_topo, "Fig_F9_centrality_correlation.pdf"), width = 8, height = 7)
corrplot(
  corr_mat,
  method = "color",
  type = "upper",
  order = "hclust",
  col = colorRampPalette(c("#2980B9", "white", "#C0392B"))(200),
  tl.col = "black",
  tl.srt = 40,
  tl.cex = 0.95,
  addCoef.col = "black",
  number.cex = 0.80,
  diag = FALSE,
  title = "Centrality — Spearman Correlation",
  mar = c(0, 0, 2, 0)
)
dev.off()
message("  [FIG] Fig_F9_centrality_correlation.pdf")

# F10. Radar chart
top_hubs <- node_tbl$gene[node_tbl$is_hub][order(-node_tbl$hub_score[node_tbl$is_hub])]
top_hubs <- top_hubs[seq_len(min(8, length(top_hubs)))]

if (length(top_hubs) >= 3) {
  radar_df <- node_tbl[node_tbl$gene %in% top_hubs,
                       c("gene", "r_deg", "r_bet", "r_clos", "r_eig", "r_pr", "r_mcc")] %>%
    pivot_longer(-gene, names_to = "metric", values_to = "val") %>%
    mutate(
      metric = recode(
        metric,
        r_deg = "Degree",
        r_bet = "Betweenness",
        r_clos = "Closeness",
        r_eig = "Eigenvector",
        r_pr = "PageRank",
        r_mcc = "MCC"
      ),
      metric = factor(metric, levels = c("Degree", "Betweenness", "Closeness", "Eigenvector", "PageRank", "MCC"))
    )
  
  p_radar <- ggplot(radar_df, aes(x = metric, y = val, group = gene, color = gene, fill = gene)) +
    geom_polygon(alpha = 0.10) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    coord_polar() +
    scale_y_continuous(limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1)) +
    scale_color_brewer(palette = "Dark2", name = "Gene") +
    scale_fill_brewer(palette = "Dark2", name = "Gene") +
    facet_wrap(~gene, ncol = 4) +
    labs(
      title = "Hub Centrality Profiles",
      subtitle = "Rank-normalised scores (0 = lowest, 1 = highest)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 8, face = "bold"),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "none"
    )
  
  save_fig(p_radar, "Fig_F10_hub_radar", fig_hub, 14, 9)
}

# ── 9. DIRECTION-STRATIFIED SUBNETWORKS ──────────────────────────────────────
message("\n[H] Direction-stratified subnetworks")

make_sub_plot <- function(keep_genes, g_full, node_tbl, title, subtitle, zoom = 1.28) {
  vkeep <- intersect(keep_genes, igraph::V(g_full)$name)
  if (length(vkeep) < 3) return(NULL)
  
  g_sub <- igraph::induced_subgraph(g_full, vkeep)
  tbl_sub <- make_tbl_graph_from_lcc(g_sub, node_tbl)
  lay_sub <- layout_tbl_graph(tbl_sub, layout = "fr")
  lay_sub <- zoom_layout(lay_sub, zoom = zoom)
  
  ggraph(lay_sub) +
    geom_edge_link(
      aes(alpha = weight, width = weight),
      colour = "#AAB3BF", show.legend = FALSE
    ) +
    scale_edge_alpha(range = c(0.08, 0.55)) +
    scale_edge_width(range = c(0.2, 1.4)) +
    geom_node_point(
      aes(fill = z_meta, size = degree, shape = is_hub),
      colour = "white", stroke = 0.5, alpha = 0.96
    ) +
    scale_shape_manual(
      values = c("FALSE" = 21, "TRUE" = 23),
      labels = c("Other", "Hub"),
      name = "Node type",
      guide = shape_legend_guide
    ) +
    scale_fill_gradient2(
      low = "#2980B9", mid = "grey92", high = "#C0392B",
      midpoint = 0, name = "Meta-Z",
      guide = fill_colorbar_guide
    ) +
    scale_size_continuous(
      range = c(3.0, 10.0),
      name = "Degree",
      guide = size_legend_guide
    ) +
    geom_node_text(
      aes(label = ifelse(is_hub | degree >= quantile(degree, 0.80), name, "")),
      repel = TRUE,
      size = 2.9,
      family = "sans",
      fontface = "bold",
      colour = "#1A1A1A",
      bg.colour = "white",
      bg.r = 0.09,
      max.overlaps = Inf,
      seed = 2025
    ) +
    labs(title = title, subtitle = subtitle) +
    theme_graph(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.box = "vertical"
    )
}

p_risk <- make_sub_plot(
  gene_df$gene[gene_df$direction == "Risk-promoting"],
  g_lcc, node_tbl,
  "Risk-Promoting Gene Subnetwork",
  paste0("z_meta > 0 | n=", sum(gene_df$direction == "Risk-promoting")),
  zoom = 1.34
)

p_prot <- make_sub_plot(
  gene_df$gene[gene_df$direction == "Protective"],
  g_lcc, node_tbl,
  "Protective Gene Subnetwork",
  paste0("z_meta < 0 | n=", sum(gene_df$direction == "Protective")),
  zoom = 1.34
)

save_fig(p_risk, "Fig_H1_risk_subnetwork", fig_sub, 12, 10)
save_fig(p_prot, "Fig_H2_protective_subnetwork", fig_sub, 12, 10)

if (!is.null(p_risk) && !is.null(p_prot)) {
  save_fig(
    p_risk + p_prot +
      plot_annotation(
        title = "Direction-Stratified PPI Subnetworks",
        theme = theme(plot.title = element_text(face = "bold", size = 13))
      ),
    "Fig_H3_dual_panel", fig_sub, 22, 10
  )
}

# ── 10. NETWORK ROBUSTNESS ────────────────────────────────────────────────────
message("\n[K] Robustness simulation")

attack_sim <- function(g, strategy = c("degree", "betweenness", "random"), steps = NULL) {
  strategy <- match.arg(strategy)
  n_orig <- igraph::vcount(g)
  if (is.null(steps)) steps <- min(n_orig, 40L)
  
  get_order <- function(gt) {
    vn <- igraph::V(gt)$name
    switch(
      strategy,
      degree      = vn[order(igraph::degree(gt), decreasing = TRUE)],
      betweenness = vn[order(igraph::betweenness(gt), decreasing = TRUE)],
      random      = sample(vn)
    )
  }
  
  removal_queue <- get_order(g)
  lcc_vec <- numeric(steps)
  gt <- g
  
  for (i in seq_len(steps)) {
    target <- NA_character_
    while (length(removal_queue) > 0) {
      cand <- removal_queue[1]
      removal_queue <- removal_queue[-1]
      if (cand %in% igraph::V(gt)$name) {
        target <- cand
        break
      }
    }
    if (is.na(target)) break
    
    gt <- igraph::delete_vertices(gt, target)
    lcc_vec[i] <- if (igraph::vcount(gt) == 0) 0 else max(igraph::components(gt)$csize)
    
    if (strategy != "random" && igraph::vcount(gt) > 0) {
      removal_queue <- get_order(gt)
    }
  }
  
  data.frame(
    step = seq_len(steps),
    frac_rem = seq_len(steps) / n_orig,
    lcc_frac = lcc_vec / n_orig,
    strategy = strategy,
    stringsAsFactors = FALSE
  )
}

rob_deg <- attack_sim(g_lcc, "degree")
rob_bet <- attack_sim(g_lcc, "betweenness")

rob_rnd_list <- lapply(seq_len(5L), function(s) {
  set.seed(s)
  attack_sim(g_lcc, "random")
})

rob_rnd_raw <- do.call(rbind, rob_rnd_list)
rob_rnd <- rob_rnd_raw %>%
  group_by(step, frac_rem) %>%
  summarise(lcc_frac = mean(lcc_frac), strategy = "random", .groups = "drop") %>%
  as.data.frame()

rob_all <- rbind(rob_deg, rob_bet, rob_rnd) %>%
  mutate(
    strategy = factor(
      strategy,
      levels = c("degree", "betweenness", "random"),
      labels = c("Targeted (degree)", "Targeted (betweenness)", "Random")
    )
  )

p_rob <- ggplot(rob_all, aes(x = frac_rem * 100, y = lcc_frac * 100,
                             color = strategy, linetype = strategy)) +
  geom_line(linewidth = 0.9, alpha = 0.90) +
  scale_color_manual(
    values = c(
      "Targeted (degree)" = "#C0392B",
      "Targeted (betweenness)" = "#E67E22",
      "Random" = "#2980B9"
    ),
    name = "Attack strategy"
  ) +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotted"),
    name = "Attack strategy"
  ) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), name = "Nodes removed (%)") +
  scale_y_continuous(labels = function(y) paste0(y, "%"), name = "LCC size (%)") +
  labs(
    title = "Network Robustness — Targeted vs Random Attack",
    caption = "Random = mean of 5 simulations"
  ) +
  q1_theme(11)

save_fig(p_rob, "Fig_K1_robustness", fig_rob, 10, 7)
write.csv(rob_all, file.path(OUT_DIR, "tables", "robustness.csv"), row.names = FALSE)

# ── 11. GRAPH EXPORT ──────────────────────────────────────────────────────────
message("\n[L] Graph export")

write.csv(node_tbl, file.path(OUT_DIR, "graph", "nodes.csv"), row.names = FALSE)
write.csv(igraph::as_data_frame(g_lcc, "edges"),
          file.path(OUT_DIR, "graph", "edges.csv"), row.names = FALSE)

igraph::write_graph(
  g_lcc,
  file.path(OUT_DIR, "graph", "ppi_network.graphml"),
  format = "graphml"
)

saveRDS(tbl_g_lcc, file.path(OUT_DIR, "graph", "ppi_tbl_graph.rds"))
message("  GraphML + CSV + tbl_graph exported")

# ── 12. EXCEL EXPORT ─────────────────────────────────────────────────────────
message("\n[M] Excel export")

wb <- createWorkbook()

hdr <- createStyle(
  fontName = "Calibri", fontSize = 10, fontColour = "white",
  fgFill = "#1F3864", halign = "center", textDecoration = "bold",
  border = "Bottom", wrapText = TRUE
)

alt <- createStyle(fgFill = "#F2F7FF")
hub_style <- createStyle(fgFill = "#FFF3CD", fontColour = "#856404",
                         textDecoration = "bold")

add_ws <- function(wb, name, df) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr)
  nr <- nrow(df)
  nc <- ncol(df)
  
  addStyle(
    wb, name,
    createStyle(fontName = "Calibri", fontSize = 9, border = "Bottom",
                borderColour = "#D9D9D9"),
    rows = 2:(nr + 1), cols = 1:nc, gridExpand = TRUE
  )
  
  if (nr >= 2) {
    addStyle(
      wb, name, alt,
      rows = seq(3, nr + 1, 2), cols = 1:nc,
      gridExpand = TRUE, stack = TRUE
    )
  }
  
  setColWidths(wb, name, 1:nc, "auto")
  freezePane(wb, name, firstRow = TRUE)
}

add_ws(wb, "Global_Metrics", global_metrics)

add_ws(
  wb, "Node_Centrality",
  node_tbl[order(-node_tbl$hub_score), ] %>%
    mutate(across(where(is.numeric), ~round(., 5)))
)

add_ws(
  wb, "Hub_Genes",
  node_tbl[node_tbl$is_hub,
           c("gene", "degree", "betweenness", "closeness", "eigenvector",
             "pagerank", "mcc", "hub_score", "n_criteria",
             "z_meta", "direction", "bag_tier")] %>%
    arrange(desc(hub_score))
)

add_ws(wb, "Community_Modules", module_tbl)

add_ws(
  wb, "STRING_Interactions",
  interactions[order(-interactions$combined_score), ][1:min(5000, nrow(interactions)), ]
)

add_ws(
  wb, "Robustness",
  rob_all %>% mutate(across(where(is.numeric), ~round(., 4)))
)

hub_rows <- which(node_tbl[order(-node_tbl$hub_score), "is_hub"]) + 1
if (length(hub_rows) > 0) {
  addStyle(
    wb, "Node_Centrality", hub_style,
    rows = hub_rows, cols = 1:ncol(node_tbl),
    gridExpand = TRUE, stack = TRUE
  )
}

xl <- file.path(OUT_DIR, "tables", "PPI_Network_Results_Q1.xlsx")
saveWorkbook(wb, xl, overwrite = TRUE)
message("  Excel: ", basename(xl))

# ── 13. SAVE RDATA ────────────────────────────────────────────────────────────
message("\n[N] Saving RData")

save(
  gene_df, gene_symbols, interactions,
  g, g_lcc, tbl_g_lcc, node_tbl, hub_genes,
  comm_louvain, comm_walk, module_tbl,
  global_metrics, rob_all,
  file = file.path(OUT_DIR, "rdata", "ppi_results.RData")
)

sink(file.path(OUT_DIR, "session_info.txt"))
cat("PPI Network Q1 v3.1 —", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()

message("\n", strrep("═", 66))
message("  PPI NETWORK COMPLETE")
message("  Root   : ", OUT_DIR)
message("  Hubs   : ", paste(sort(hub_genes), collapse = ", "))
message("  Figures: ", file.path(OUT_DIR, "figures"))
message("  Graph  : ", file.path(OUT_DIR, "graph"))
message("  Excel  : ", xl)
message(strrep("═", 66))
