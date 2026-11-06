# ══════════════════════════════════════════════════════════════════════════
# One heatmap per Reactome pathway for the JT-INCREASING proteins (table-19) and the JT-DECREASING
# proteins (table-20): the TOP_N pathways of each table (default 50), with significance stars versus
# Control inside the cells; plus all six pairwise group comparisons for all proteins.

#
# INPUTS    session objects expr_sub, pheno (and optionally res_jt) created by the proteomics JT script
#           (jonckheere_test_control_res_vol_ad.R -- NOT part of this release, see README);
#           data/reactome/Reactome_Pathways_2024_table-19.txt and table-20.txt
# OUTPUTS   figures/pathway_heatmaps/<increasing|decreasing>/*.png; results/pathway_heatmap_source_data.csv,
#           results/stats_vs_control_all_proteins.csv, results/pairwise_group_stats_adjP_lt_0.05.csv
#
#
#   Top N  = the TOP_N smallest adjusted P-values in each Reactome table.
#   Rows   = the proteins listed for that pathway in the Reactome table
#            (JT-significant proteins that belong to the pathway).
#   Columns= Control | Resilient AsymAD | Vulnerable AsymAD | AD.
#   Colour = mean of the per-protein z-score over the donors of that group.
#            Each protein is z-scored across ALL donors first, so the cell is
#            "how many SD above/below the cohort average this group sits".
#            One colour scale for all 20 heatmaps. Rows are clustered.
#   Stars  = each group vs Control, per protein, from limma:
#                expr ~ group + age + sex + APOE + PMI   (Control = reference)
#            fitted on ALL proteins in expr_sub, BH-adjusted across all of them
#            separately for each comparison.
#            * adj.P < 0.05   ** < 0.01   *** < 0.001   (no star = not significant)
#            The Control column is the reference, so it carries no stars.
#            NOTE: the colour is the raw group mean; the star is covariate-adjusted.
#   Names  : every protein is labelled up to MAX_LABELS rows; larger pathways
#            label only the N_MARK most JT-significant proteins.
#
# Requires in the session: expr_sub (proteins x samples, columns = sample_id),
#   pheno (sample_id, group_label, and age_death_z / msex / apoe_risk / pmi_z
#   for the adjusted test). res_jt is used only to pick which names to label.
# ══════════════════════════════════════════════════════════════════════════

source("R/00_setup.R")
require_packages(c("dplyr", "tidyr", "tibble", "readr", "limma", "ComplexHeatmap", "circlize"))

library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(limma)
library(ComplexHeatmap)
library(circlize)
library(grid)

REACTOME_FILE <- c(increasing = "data/reactome/Reactome_Pathways_2024_table-19.txt",
                   decreasing = "data/reactome/Reactome_Pathways_2024_table-20.txt")
TOP_N        <- 50
OUT_DIR      <- "figures/pathway_heatmaps"
GROUP_LEVELS <- c("Control", "Resilient AsymAD", "Vulnerable AsymAD", "AD")
MAX_LABELS   <- 70    # label every protein up to this many rows
N_MARK       <- 70    # bigger pathways: label only this many (most significant)
COVARIATES   <- c("age_death_z", "msex", "apoe_risk", "pmi_z")

stopifnot(exists("expr_sub"), exists("pheno"))

# ── 0. donors and groups ────────────────────────────────────────────────────
stopifnot(all(c("sample_id", "group_label") %in% names(pheno)))
bad <- setdiff(unique(as.character(pheno$group_label)), GROUP_LEVELS)
if (length(bad) > 0)
  stop("pheno$group_label has values not in GROUP_LEVELS: ", paste(bad, collapse = ", "))

grp <- pheno %>%
  transmute(sample_id, group = factor(as.character(group_label), levels = GROUP_LEVELS)) %>%
  filter(!is.na(group), sample_id %in% colnames(expr_sub))

n_per_group <- table(grp$group)
stopifnot(all(n_per_group[GROUP_LEVELS] > 0))
col_labels <- paste0(GROUP_LEVELS, "\n(n = ", as.integer(n_per_group[GROUP_LEVELS]), ")")
cat("Donors per group:\n"); print(n_per_group)

# ── 1. top-10 pathways per direction ────────────────────────────────────────
read_top <- function(direction) {
  read_tsv(REACTOME_FILE[[direction]], show_col_types = FALSE) %>%
    rename(adj_p = `Adjusted P-value`) %>%
    arrange(adj_p) %>%
    slice_head(n = TOP_N) %>%
    mutate(rank = row_number(), direction = direction)
}
pw <- bind_rows(lapply(names(REACTOME_FILE), read_top))

cat("\nPathways used:\n")
print(as.data.frame(pw %>% select(direction, rank, Term, Overlap, adj_p)))

# ── 2. per-protein z-score (all donors) -> mean per group ──────────────────
row_key <- setNames(rownames(expr_sub), toupper(rownames(expr_sub)))
match_rows <- function(g) {
  hit <- unname(row_key[toupper(trimws(g))])
  unique(hit[!is.na(hit)])
}

rows_needed <- unique(unlist(lapply(strsplit(pw$Genes, ";"), match_rows)))
cat("\nProteins across all 20 pathways found in expr_sub:", length(rows_needed), "\n")

X <- expr_sub[rows_needed, grp$sample_id, drop = FALSE]
Z <- t(scale(t(X)))                                   # z across ALL donors, per protein

group_means <- sapply(GROUP_LEVELS, function(g)
  rowMeans(Z[, grp$sample_id[grp$group == g], drop = FALSE], na.rm = TRUE))
good <- rownames(group_means)[rowSums(is.finite(group_means)) == length(GROUP_LEVELS)]
cat("Proteins with a finite mean in all 4 groups:", length(good), "of", nrow(group_means), "\n")

# one colour scale for every heatmap
LIM <- max(0.3, ceiling(quantile(abs(group_means[good, ]), 0.98) * 10) / 10)
col_fun <- colorRamp2(c(-LIM, 0, LIM), c("#2166AC", "#F7F7F7", "#B2182B"))
cat("Colour scale: +/-", LIM, "SD\n")

# ── 2b. each group vs Control: limma with covariates, ALL proteins, BH ─────
covars <- intersect(COVARIATES, names(pheno))
if (length(setdiff(COVARIATES, covars)) > 0)
  message("Covariates not found in pheno, model runs WITHOUT: ",
          paste(setdiff(COVARIATES, covars), collapse = ", "))

d <- grp %>% left_join(pheno %>% select(sample_id, all_of(covars)), by = "sample_id")
for (v in intersect(c("msex", "apoe_risk"), covars)) d[[v]] <- factor(d[[v]])
d <- d[complete.cases(d), ]                           # model.matrix would silently drop NA rows
stopifnot(all(GROUP_LEVELS %in% as.character(d$group)))
d$group <- factor(as.character(d$group), levels = GROUP_LEVELS)   # Control = reference

design <- model.matrix(as.formula(paste("~", paste(c("group", covars), collapse = " + "))), data = d)
coef_names <- paste0("group", GROUP_LEVELS[-1])
stopifnot(all(coef_names %in% colnames(design)))
cat("Stars model: ~", paste(c("group", covars), collapse = " + "),
    "| donors:", nrow(d), "| proteins:", nrow(expr_sub), "\n")

fit <- eBayes(lmFit(expr_sub[, d$sample_id, drop = FALSE], design))

res_vs_ctrl <- bind_rows(lapply(GROUP_LEVELS[-1], function(g) {
  topTable(fit, coef = paste0("group", g), number = Inf,
           adjust.method = "BH", sort.by = "none") %>%
    rownames_to_column("protein") %>%
    transmute(protein, group = g, logFC, P.Value, adj.P.Val)
}))

cat("Significant vs Control (BH < 0.05, all proteins):\n")
print(res_vs_ctrl %>% group_by(group) %>% summarise(n_sig = sum(adj.P.Val < 0.05, na.rm = TRUE)))

adj_wide <- res_vs_ctrl %>% select(protein, group, adj.P.Val) %>%
  pivot_wider(names_from = group, values_from = adj.P.Val)

stats_wide <- res_vs_ctrl %>%
  mutate(g = gsub(" ", "_", group)) %>%
  pivot_wider(id_cols = protein, names_from = g,
              values_from = c(logFC, adj.P.Val), names_sep = "_vs_Control_")

star_fun <- function(p) ifelse(is.na(p), "",
                               ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

# draws the stars, always black; text size follows the cell height
make_cell_fun <- function(S) {
  force(S)
  function(j, i, x, y, width, height, fill) {
    s <- S[i, j]
    if (!is.na(s) && nzchar(s)) {
      fs <- min(12, convertHeight(height, "pt", valueOnly = TRUE) * 0.9)
      grid.text(s, x, y - unit(fs * 0.25, "pt"),
                gp = gpar(fontsize = fs, fontface = "bold", col = "black"))
    }
  }
}

jt_p <- if (exists("res_jt")) setNames(res_jt$adj.P.Val, res_jt$protein) else NULL

# ── 3. one heatmap per pathway ─────────────────────────────────────────────
for (dd in names(REACTOME_FILE)) dir.create(file.path(OUT_DIR, dd), recursive = TRUE, showWarnings = FALSE)
source_rows <- list()

for (i in seq_len(nrow(pw))) {
  term      <- pw$Term[i]
  direction <- pw$direction[i]
  listed    <- strsplit(pw$Genes[i], ";")[[1]]
  genes     <- base::intersect(match_rows(listed), good)
  n         <- length(genes)
  cat(sprintf("%-10s %2d  %-50s listed %3d | in data %3d | plotted %3d\n",
              direction, pw$rank[i], term, length(listed), length(match_rows(listed)), n))
  if (n < 2) { cat("   SKIPPED (fewer than 2 proteins)\n"); next }
  
  M <- group_means[genes, , drop = FALSE]
  colnames(M) <- col_labels
  
  # stars matrix, same row/column order as M (Control column stays empty)
  S <- matrix("", nrow = n, ncol = length(GROUP_LEVELS), dimnames = list(genes, GROUP_LEVELS))
  A <- as.matrix(adj_wide[match(genes, adj_wide$protein), GROUP_LEVELS[-1]])
  S[, GROUP_LEVELS[-1]] <- star_fun(A)
  
  show_names <- n <= MAX_LABELS
  right_anno <- NULL
  if (!show_names) {
    top_lab <- if (!is.null(jt_p)) names(sort(jt_p[genes]))[seq_len(min(N_MARK, n))] else genes[seq_len(min(N_MARK, n))]
    right_anno <- rowAnnotation(
      labels = anno_mark(at = match(top_lab, rownames(M)), labels = top_lab,
                         labels_gp = gpar(fontsize = 10, fontface = "bold")))
  }
  
  body_h <- min(10, max(2.2, n * 0.17))
  
  ht <- Heatmap(M, name = "z", col = col_fun,
                cell_fun = make_cell_fun(S),
                cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
                show_row_names = show_names,
                row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                column_names_gp = gpar(fontsize = 11, fontface = "bold"),
                column_names_rot = 45,
                column_title = paste0(paste(strwrap(term, 40), collapse = "\n"),
                                      sprintf("\n%s | %d proteins | adj.P = %.1e",
                                              direction, n, pw$adj_p[i])),
                column_title_gp = gpar(fontsize = 13, fontface = "bold"),
                rect_gp = gpar(col = "white", lwd = 0.5),
                right_annotation = right_anno,
                width = unit(3.6, "in"), height = unit(body_h, "in"),
                heatmap_legend_param = list(title = "Mean z-score\n(per protein,\nall donors)",
                                            title_gp = gpar(fontsize = 10, fontface = "bold"),
                                            labels_gp = gpar(fontsize = 10, fontface = "bold")))
  
  safe  <- gsub("[^A-Za-z0-9]+", "_", term)
  fname <- file.path(OUT_DIR, direction, sprintf("%02d_%s.png", pw$rank[i], safe))
  png(fname, width = 8.5, height = body_h + 3.7, units = "in", res = 300)
  draw(ht, padding = unit(c(12, 5, 5, 5), "mm"))
  grid.text("Stars: group vs Control, limma adjusted for age, sex, APOE, PMI; BH across all proteins.  * <0.05  ** <0.01  *** <0.001",
            x = unit(5, "mm"), y = unit(3, "mm"), just = c("left", "bottom"),
            gp = gpar(fontsize = 8, col = "grey30"))
  dev.off()
  
  source_rows[[length(source_rows) + 1]] <- as.data.frame(M) %>%
    setNames(GROUP_LEVELS) %>%
    rownames_to_column("protein") %>%
    left_join(stats_wide, by = "protein") %>%
    mutate(direction = direction, rank = pw$rank[i], pathway = term,
           JT_adj_p = if (!is.null(jt_p)) unname(jt_p[protein]) else NA_real_,
           .before = 1)
}

write.csv(bind_rows(source_rows), file.path(RES_DIR, "pathway_heatmap_source_data.csv"), row.names = FALSE)
write.csv(res_vs_ctrl, file.path(RES_DIR, "stats_vs_control_all_proteins.csv"), row.names = FALSE)
cat("\nHeatmaps saved under", OUT_DIR, "(increasing/ and decreasing/)\n")



dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PAIRS <- list(c("Resilient AsymAD",  "Control"),
              c("Vulnerable AsymAD", "Control"),
              c("AD",                "Control"),
              c("Vulnerable AsymAD", "Resilient AsymAD"),
              c("AD",                "Resilient AsymAD"),
              c("AD",                "Vulnerable AsymAD"))

# contrast = effect(X) - effect(Y); Control is the reference, so its effect is 0
effect_vec <- function(g) {
  v <- setNames(numeric(ncol(design)), colnames(design))
  if (g != GROUP_LEVELS[1]) v[paste0("group", g)] <- 1
  v
}
cm <- sapply(PAIRS, function(p) effect_vec(p[1]) - effect_vec(p[2]))
colnames(cm) <- sapply(PAIRS, function(p) paste(p[1], "vs", p[2]))

fit_lm <- lmFit(expr_sub[, d$sample_id, drop = FALSE], design)
fit_pw <- eBayes(contrasts.fit(fit_lm, cm))

pairwise <- bind_rows(lapply(colnames(cm), function(cn) {
  topTable(fit_pw, coef = cn, number = Inf, adjust.method = "BH", sort.by = "none") %>%
    rownames_to_column("protein") %>%
    transmute(comparison = cn, protein, logFC, AveExpr, t, P.Value, adj.P.Val)
})) %>%
  mutate(comparison = factor(comparison, levels = colnames(cm)))

# sanity check: the vs-Control rows should reproduce the stars' statistics
chk <- pairwise %>% filter(comparison == "AD vs Control") %>%
  inner_join(res_vs_ctrl %>% filter(group == "AD"), by = "protein", suffix = c("", "_stars"))
cat("Max |adj.P difference| vs the star model (AD vs Control):",
    signif(max(abs(chk$adj.P.Val - chk$adj.P.Val_stars), na.rm = TRUE), 3), "\n")

# which of the 20 heatmap pathways each protein belongs to
pw_membership <- bind_rows(lapply(seq_len(nrow(pw)), function(i)
  tibble(protein = match_rows(strsplit(pw$Genes[i], ";")[[1]]),
         pathway = paste0(pw$direction[i], ": ", pw$Term[i])))) %>%
  group_by(protein) %>%
  summarise(heatmap_pathways = paste(pathway, collapse = " | "), .groups = "drop")

sig_pairwise <- pairwise %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05) %>%
  mutate(higher_in = ifelse(logFC > 0,
                            sub(" vs .*", "", comparison),
                            sub(".* vs ", "", comparison))) %>%
  left_join(pw_membership, by = "protein") %>%
  arrange(comparison, adj.P.Val)

write.csv(sig_pairwise, file.path(RES_DIR, "pairwise_group_stats_adjP_lt_0.05.csv"), row.names = FALSE)

cat("\nSignificant proteins (BH < 0.05, tested across all proteins) per comparison:\n")
print(pairwise %>% group_by(comparison) %>%
        summarise(n_sig = sum(adj.P.Val < 0.05, na.rm = TRUE), .groups = "drop"))
cat("Saved", nrow(sig_pairwise), "rows to",
    file.path(RES_DIR, "pairwise_group_stats_adjP_lt_0.05.csv"), "\n")
