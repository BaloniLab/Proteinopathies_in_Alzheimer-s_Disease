# ══════════════════════════════════════════════════════════════════════════
# Protein x metabolite Spearman correlation networks (dual-omics matched donors, AsymAD + AD)

#
#   Donors with BOTH proteomics and metabolomics (matched through the metadata crosswalk; the two
#   omics use different donor ids). Spearman correlation between differentially abundant proteins
#   (adj.P < 0.05, script 05) and metabolites (script 08; adj.P < 0.05 in the main network, < 0.1 in the
#   relaxed network); pairs with n >= 10, BH across pairs; edges = |rho| > 0.50 (main) / 0.40 (relaxed);
#   Cytoscape edge / node tables and GraphML files. Implemented as one function, run_network(),
#   called for the main and for the relaxed definition.
#
# INPUTS
#   R objects: meta_cc_8reg_n323_2026_09_02 and the session objects created by scripts 05 and 08:
#   protein_expr, res (proteins), expr_metab_filt, res_metab_named (metabolites)
# OUTPUTS
#   results/cor_prot_metab_long_DE_dualomics_matched[_0.1_p].csv
#   results/cytoscape_edges_/nodes_DE_dualomics_matched*.csv
#   results/network_AsymAD_vs_AD_DE_dualomics_matched[_relaxed].graphml
# RUN AFTER
#   05_proteomics_DEG_full_cohort.R and 08_metabolomics_DEG_full_cohort.R (same R session)
#
# Run from the project root (see README.md). The analysis scripts keep their
# intermediate objects in the R session; run_all.R sources them in order.
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")
require_packages(c("dplyr", "tibble", "reshape2", "igraph"))
load_inputs(c("meta_cc_8reg_n323_2026_09_02"))

library(dplyr)
library(tibble)
library(reshape2)
library(igraph)

for (obj in c("protein_expr", "res", "expr_metab_filt", "res_metab_named"))
  if (!exists(obj)) stop("object `", obj, "` not found: run scripts 05 and 08 first (same R session)")

# -----------------------------------------------------------------------------
# 0. crosswalk: which donors have BOTH a proteomics column and a metabolomics
#    column, and what their id is under each scheme
# -----------------------------------------------------------------------------
id_crosswalk <- meta_cc_8reg_n323_2026_09_02 %>%
  filter(EmoryStrictDx.2019 %in% c("AD", "AsymAD")) %>%
  transmute(
    prot_id   = as.character(...1),     # matches expr_filt's colnames
    metab_id  = as.character(projid),   # matches expr_metab_filt's colnames
    diagnosis = factor(EmoryStrictDx.2019, levels = c("AD", "AsymAD"))
  )

dual_ids <- id_crosswalk %>%
  filter(prot_id %in% colnames(protein_expr),
         metab_id %in% colnames(expr_metab_filt))
expr_filt <- protein_expr
cat("Donors with BOTH proteomics and metabolomics:", nrow(dual_ids), "\n")
print(table(dual_ids$diagnosis))

# -----------------------------------------------------------------------------
# 1. build the two matrices with IDENTICAL sample order (by donor, via the
#    crosswalk) — then relabel metabolomics' columns to the matching
#    proteomics id, purely so the two matrices carry one shared, unambiguous
#    sample label for everything downstream (the two ID schemes never get
#    compared as strings anywhere else in this script)
# -----------------------------------------------------------------------------
expr_prot_matched  <- expr_filt[, dual_ids$prot_id, drop = FALSE]
expr_metab_matched <- expr_metab_filt[, dual_ids$metab_id, drop = FALSE]

stopifnot(ncol(expr_prot_matched) == ncol(expr_metab_matched))

colnames(expr_metab_matched) <- colnames(expr_prot_matched)  # now both = prot_id, order-matched
stopifnot(identical(colnames(expr_prot_matched), colnames(expr_metab_matched)))

# -----------------------------------------------------------------------------
# 2-7. correlation network, run once per definition of "differentially abundant"
#   metab_adj_cut : adj.P.Val cut-off for the metabolites entering the network
#   rho_cut       : |rho| above which a protein-metabolite pair becomes an edge
#   prot_adj_cut  : adj.P.Val cut-off for the proteins (0.05 in both networks)
# -----------------------------------------------------------------------------
run_network <- function(metab_adj_cut, rho_cut, prot_adj_cut = 0.05,
                        f_long, f_edges, f_nodes, f_graph) {

  sig_proteins    <- res             %>% filter(adj.P.Val < prot_adj_cut)  %>% pull(protein)
  sig_metabolites <- res_metab_named %>% filter(adj.P.Val < metab_adj_cut) %>% pull(metabolite)
  cat("DE proteins (adj.P.Val <", prot_adj_cut, "):", length(sig_proteins), "\n")
  cat("DE metabolites (adj.P.Val <", metab_adj_cut, "):", length(sig_metabolites), "\n")

  # subset the dual-omics-matched matrices to the DE features
  prot_mat_aligned  <- expr_prot_matched[rownames(expr_prot_matched) %in% sig_proteins, , drop = FALSE]
  metab_mat_aligned <- expr_metab_matched[rownames(expr_metab_matched) %in% sig_metabolites, , drop = FALSE]
  cat("Matched samples:", ncol(prot_mat_aligned), "\n")
  cat("Proteins in network:", nrow(prot_mat_aligned), " Metabolites in network:", nrow(metab_mat_aligned), "\n")
  stopifnot(ncol(prot_mat_aligned) == ncol(metab_mat_aligned))

  # Spearman correlation matrix + n / p-value per pair
  cor_matrix <- cor(t(prot_mat_aligned), t(metab_mat_aligned),
                    method = "spearman", use = "pairwise.complete.obs")

  ind_prot  <- !is.na(prot_mat_aligned)
  ind_metab <- !is.na(metab_mat_aligned)
  n_mat     <- ind_prot %*% t(ind_metab)

  t_stat <- cor_matrix * sqrt((n_mat - 2) / (1 - cor_matrix^2))
  p_mat  <- 2 * pt(-abs(t_stat), df = n_mat - 2)
  p_mat[n_mat < 3] <- NA

  cor_long <- melt(cor_matrix, varnames = c("protein", "metabolite"), value.name = "rho") %>%
    left_join(melt(n_mat, varnames = c("protein", "metabolite"), value.name = "n"),
              by = c("protein", "metabolite")) %>%
    left_join(melt(p_mat, varnames = c("protein", "metabolite"), value.name = "p_value"),
              by = c("protein", "metabolite")) %>%
    filter(!is.na(rho), n >= 10) %>%
    mutate(adj_p     = p.adjust(p_value, method = "BH"),
           direction = ifelse(rho > 0, "Positive", "Negative"),
           abs_rho   = abs(rho)) %>%
    arrange(adj_p)

  cat("Pairs after n>=10 filter:", nrow(cor_long), "\n")
  write.csv(cor_long, file.path(RES_DIR, f_long), row.names = FALSE)

  # network edges for Cytoscape
  edge_list <- cor_long %>%
    filter(abs_rho > rho_cut) %>%
    select(source = protein, target = metabolite, rho, abs_rho, direction, n, p_value, adj_p)
  write.csv(edge_list, file.path(RES_DIR, f_edges), row.names = FALSE)
  cat("Edges (|rho| >", rho_cut, "):", nrow(edge_list), "\n")

  # node table: DE direction / magnitude from the full-cohort DEG results (res, res_metab_named)
  node_proteins <- res %>%
    filter(protein %in% unique(edge_list$source)) %>%
    transmute(id = protein, type = "Protein", logFC, adj.P.Val,
              sig = case_when(adj.P.Val < prot_adj_cut & logFC > 0 ~ "Up",
                              adj.P.Val < prot_adj_cut & logFC < 0 ~ "Down",
                              TRUE ~ "NS"))

  node_metabolites <- res_metab_named %>%
    filter(metabolite %in% unique(edge_list$target)) %>%
    transmute(id = metabolite, type = "Metabolite", logFC, adj.P.Val,
              sig = case_when(adj.P.Val < metab_adj_cut & logFC > 0 ~ "Up",
                              adj.P.Val < metab_adj_cut & logFC < 0 ~ "Down",
                              TRUE ~ "NS"))

  node_list <- bind_rows(node_proteins, node_metabolites)
  write.csv(node_list, file.path(RES_DIR, f_nodes), row.names = FALSE)
  cat("Nodes:", nrow(node_list), " (", nrow(node_proteins), "proteins,", nrow(node_metabolites), "metabolites)\n")

  # GraphML export
  g_export <- graph_from_data_frame(
    d        = edge_list %>% select(from = source, to = target, rho, abs_rho, direction, n, p_value, adj_p),
    vertices = node_list %>% rename(name = id),
    directed = FALSE
  )
  V(g_export)$degree <- degree(g_export)
  write_graph(g_export, file.path(RES_DIR, f_graph), format = "graphml")
  cat("Nodes:", vcount(g_export), " Edges:", ecount(g_export), "\nSaved:", f_graph, "\n\n")

  invisible(list(cor_long = cor_long, edge_list = edge_list, node_list = node_list))
}

# main network: metabolites adj.P < 0.05, edges |rho| > 0.50
net_main <- run_network(metab_adj_cut = 0.05, rho_cut = 0.50,
                        f_long  = "cor_prot_metab_long_DE_dualomics_matched.csv",
                        f_edges = "cytoscape_edges_DE_dualomics_matched.csv",
                        f_nodes = "cytoscape_nodes_DE_dualomics_matched.csv",
                        f_graph = "network_AsymAD_vs_AD_DE_dualomics_matched.graphml")

# relaxed network: metabolites adj.P < 0.1, edges |rho| > 0.40
net_relaxed <- run_network(metab_adj_cut = 0.10, rho_cut = 0.40,
                           f_long  = "cor_prot_metab_long_DE_dualomics_matched_0.1_p.csv",
                           f_edges = "cytoscape_edges_DE_dualomics_matched_0.1_p.csv",
                           f_nodes = "cytoscape_nodes_DE_dualomics_matched_p_0.1.csv",
                           f_graph = "network_AsymAD_vs_AD_DE_dualomics_matched_relaxed.graphml")
