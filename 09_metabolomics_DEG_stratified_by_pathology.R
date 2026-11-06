# ══════════════════════════════════════════════════════════════════════════
# Metabolomics differential abundance, AsymAD vs AD, within pathology strata (Johnson rubric)

#
#   Six independent limma models (diagnosis + age, sex, APOE genotype, PMI; AD = reference) within
#   CERAD stage 2 / 3, Braak III+IV / V and TDP-43 stage 1+2 / 3, with an 80% presence filter
#   computed within each stratum; stacked count plot; Venn / UpSet overlap across pathologies;
#   logFC heatmaps (rows = metabolites, then the flipped layout) for metabolites significant
#   (raw p < 0.05, same direction) in at least MIN_STRATA strata.
#
# INPUTS
#   R objects: Expr_metabolomics, dataset_1744_cross_sectional_08_30_2026 (see data/README.md)
#   results/metadata_metabolomics_johnson_with_pathology.csv (written by script 08)
# OUTPUTS
#   results/limma_metab_AsymAD_vs_AD_<stratum>.csv, results/limma_metab_AsymAD_vs_AD_all_strata_combined.csv
#   results/metab_johnson_*.txt
#   figures/deg_metab_*, venn_metab_*, upset_metab_*, heatmap_metab_multistratum_*
# RUN AFTER
#   08_metabolomics_DEG_full_cohort.R (writes the metadata CSV read here)
#
# Run from the project root (see README.md). The analysis scripts keep their
# intermediate objects in the R session; run_all.R sources them in order.
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")
require_packages(c("dplyr", "tidyr", "purrr", "limma", "tibble", "ggplot2", "ggvenn", "patchwork", "ComplexHeatmap"))
load_inputs(c("Expr_metabolomics", "dataset_1744_cross_sectional_08_30_2026"))

library(dplyr)
library(tidyr)
library(purrr)
library(limma)
library(tibble)
library(ggplot2)
library(ggvenn)
library(patchwork)
library(ComplexHeatmap)

DX_LEVELS <- c("AD", "AsymAD")

if ("...1" %in% names(Expr_metabolomics)) {
  Expr_metabolomics <- Expr_metabolomics %>% column_to_rownames(var = "...1")
}

# Unnamed "X - #####" features are removed at the source (as in scripts 07 and 08); a no-op when those
# scripts have just run in the same session.
x_cols <- grep("^X\\s*-\\s*\\d+", colnames(Expr_metabolomics), value = TRUE)
Expr_metabolomics <- Expr_metabolomics[, !colnames(Expr_metabolomics) %in% x_cols, drop = FALSE]

# ── 0. NEW metadata source: Johnson rubric diagnosis + dataset_1744 staging/
#      covariates, restricted to the matched complete-case cohort ─────────
meta_strata_source <- dataset_1744_cross_sectional_08_30_2026 %>%
  mutate(projid = as.character(projid)) %>%
  transmute(projid, age_death, msex, apoe_genotype, pmi,
            CERAD_std = 4 - ceradsc,   # reverse-coded raw ceradsc -> ascending severity
            braaksc, tdp_st4)

# Cohort metadata written by script 08 (read from disk if it is not already in the session)
if (!exists("metadata_metabolomics_johnson_with_pathology"))
  metadata_metabolomics_johnson_with_pathology <- read.csv(file.path(RES_DIR, "metadata_metabolomics_johnson_with_pathology.csv"))

meta_johnson <- metadata_metabolomics_johnson_with_pathology %>%
  filter(diagnosis %in% c("AD", "AsymAD")) %>%
  mutate(projid = as.character(projid)) %>%
  transmute(projid, diagnosis = factor(diagnosis, levels = DX_LEVELS)) %>%
  left_join(meta_strata_source, by = "projid") %>%
  filter(!is.na(age_death))   # restrict to the matched complete-case cohort (~326 donors)

cat("Johnson-rubric metadata built: n =", nrow(meta_johnson), "\n")
print(table(meta_johnson$diagnosis))

# ── 1. CHECK COUNTS FIRST — do this before trusting any of the 6 models ────
cat("\n=== Sample counts per stratum (metabolomics-matched donors only) ===\n\n")
check_counts <- function(stage_col, stage_values, label) {
  n <- meta_johnson %>%
    filter(.data[[stage_col]] %in% stage_values,
           projid %in% rownames(Expr_metabolomics)) %>%
    count(diagnosis)
  cat(label, ":\n"); print(n); cat("\n")
}
check_counts("CERAD_std", 2,      "CERAD stage 2")
check_counts("CERAD_std", 3,      "CERAD stage 3")
check_counts("braaksc",   c(3,4), "Braak Limbic (III+IV)")
check_counts("braaksc",   5,      "Braak Neocortical (V)")
check_counts("tdp_st4",   c(1,2), "TDP-43 Early (1+2)")
check_counts("tdp_st4",   3,      "TDP-43 Late (3)")

cat("=== Review the counts above before proceeding ===\n",
    "Rule of thumb: a 5-covariate model (diagnosis + age + sex + apoe + pmi)\n",
    "needs comfortably more than ~15-20 total samples per stratum to be stable.\n\n")

# ── 2. one function: build pheno for a stratum, 80% filter (computed FRESH,
#      per stratum, from raw Expr_metabolomics), fit limma ─────────────────
run_deg_stratum_metab <- function(stratum_name, stage_col, stage_values, min_pct = 0.80) {
  pheno <- meta_johnson %>%
    filter(.data[[stage_col]] %in% stage_values,
           projid %in% rownames(Expr_metabolomics)) %>%
    mutate(age_death_z   = as.numeric(scale(age_death)),
           pmi_z          = as.numeric(scale(pmi)),
           msex           = factor(msex),
           apoe_genotype  = factor(apoe_genotype))
  
  cat("\n==", stratum_name, "==\n")
  cat("n =", nrow(pheno), "\n")
  print(table(pheno$diagnosis))
  
  if (nrow(pheno) < 15 || n_distinct(pheno$diagnosis) < 2) {
    cat("SKIPPED:", stratum_name, "— too few samples or only one diagnosis present.\n")
    return(NULL)
  }
  
  expr <- t(Expr_metabolomics[pheno$projid, , drop = FALSE])
  stopifnot(identical(colnames(expr), pheno$projid))
  
  asymad_ids <- pheno$projid[pheno$diagnosis == "AsymAD"]
  ad_ids     <- pheno$projid[pheno$diagnosis == "AD"]
  pct_present <- function(x) mean(!is.na(x))
  pass_asymad <- apply(expr[, asymad_ids, drop = FALSE], 1, pct_present) >= min_pct
  pass_ad     <- apply(expr[, ad_ids,     drop = FALSE], 1, pct_present) >= min_pct
  keep_metabolites <- rownames(expr)[pass_asymad | pass_ad]
  expr_filt_s <- expr[keep_metabolites, , drop = FALSE]
  cat("Metabolites kept (80% presence, union, computed WITHIN this stratum):",
      length(keep_metabolites), "/", nrow(expr), "\n")
  
  if (length(keep_metabolites) < 5) {
    cat("SKIPPED:", stratum_name, "— fewer than 5 metabolites survive the presence filter.\n")
    return(NULL)
  }
  
  design <- model.matrix(~ diagnosis + age_death_z + msex + apoe_genotype + pmi_z, data = pheno)
  cat("Design matrix condition number (kappa):", round(kappa(design), 1), "\n")
  
  fit <- lmFit(expr_filt_s, design)
  fit <- eBayes(fit)
  
  res <- topTable(fit, coef = "diagnosisAsymAD", number = Inf, adjust.method = "BH") %>%
    tibble::rownames_to_column("metabolite") %>%
    arrange(adj.P.Val) %>%
    mutate(stratum = stratum_name)
  
  fname <- file.path(RES_DIR, paste0("limma_metab_AsymAD_vs_AD_", stratum_name, ".csv"))
  write.csv(res %>% select(-stratum), fname, row.names = FALSE)
  cat("Saved:", fname, "\n")
  cat("Significant (raw p<0.05):", sum(res$P.Value < 0.05, na.rm = TRUE),
      " | (BH adj.P<0.05):", sum(res$adj.P.Val < 0.05, na.rm = TRUE), "\n")
  
  res
}

# ── 3. the six strata ────────────────────────────────────────────────────────
strata_specs <- list(
  list(name = "CERAD_stage2",         stage_col = "CERAD_std", stage_values = 2),
  list(name = "CERAD_stage3",         stage_col = "CERAD_std", stage_values = 3),
  list(name = "Braak_Limbic_III_IV",  stage_col = "braaksc",   stage_values = c(3, 4)),
  list(name = "Braak_Neocortical_V",  stage_col = "braaksc",   stage_values = 5),
  list(name = "TDP43_Early_1_2",      stage_col = "tdp_st4",   stage_values = c(1, 2)),
  list(name = "TDP43_Late_3",         stage_col = "tdp_st4",   stage_values = 3)
)

all_res_metab <- map_dfr(strata_specs, function(s)
  run_deg_stratum_metab(s$name, s$stage_col, s$stage_values))

write.csv(all_res_metab, file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_all_strata_combined.csv"), row.names = FALSE)

# ── 4. stacked barplot ───────────────────────────
if (nrow(all_res_metab) > 0) {
  summary_counts <- all_res_metab %>%
    group_by(stratum) %>%
    summarise(
      Up_raw   = sum(logFC > 0 & P.Value  < 0.05, na.rm = TRUE),
      Down_raw = sum(logFC < 0 & P.Value  < 0.05, na.rm = TRUE),
      Up_adj   = sum(logFC > 0 & adj.P.Val < 0.05, na.rm = TRUE),
      Down_adj = sum(logFC < 0 & adj.P.Val < 0.05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-stratum, names_to = c("Direction", "Threshold"), names_sep = "_",
                 values_to = "n") %>%
    mutate(Threshold = recode(Threshold, raw = "p<0.05 (raw)", adj = "adj.P<0.05 (BH)"),
           Threshold = factor(Threshold, levels = c("p<0.05 (raw)", "adj.P<0.05 (BH)")),
           Direction = factor(Direction, levels = c("Up", "Down")),
           stratum   = factor(stratum, levels = sapply(strata_specs, `[[`, "name")))
  
  p <- ggplot(summary_counts, aes(x = stratum, y = n, fill = Direction)) +
    geom_col(position = "stack", width = 0.65) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_stack(vjust = 0.5),
              fontface = "bold", size = 3.5, color = "white") +
    facet_wrap(~ Threshold, ncol = 1) +
    scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue"), name = "Direction") +
    labs(x = NULL, y = "Number of DE metabolites",
         title = "DE metabolite counts across pathology strata (AsymAD vs AD, Johnson rubric)") +
    theme_bw(base_size = 13) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
          axis.text.y      = element_text(face = "bold", color = "black"),
          axis.title       = element_text(face = "bold", color = "black"),
          strip.text       = element_text(face = "bold", color = "black"),
          strip.background = element_rect(fill = "grey92", color = NA),
          plot.title       = element_text(face = "bold", hjust = 0.5),
          legend.title     = element_text(face = "bold"),
          legend.text      = element_text(face = "bold"))
  
  p
  dir.create("figures", showWarnings = FALSE)
  ggsave("figures/deg_metab_counts_stacked_barplot_johnson.png", p, width = 10, height = 9, dpi = 600, bg = "white")
} else {
  cat("\nNo strata produced results — check the counts printed above.\n")
}

# ── 5. 4 Venn diagrams + UpSet across severity/direction ───────
files <- list(
  CERAD_stage2        = file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_CERAD_stage2.csv"),
  CERAD_stage3        = file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_CERAD_stage3.csv"),
  Braak_Limbic_III_IV = file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_Braak_Limbic_III_IV.csv"),
  Braak_Neocortical_V = file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_Braak_Neocortical_V.csv"),
  TDP43_Early_1_2     = file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_TDP43_Early_1_2.csv"),
  TDP43_Late_3        = file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_TDP43_Late_3.csv")
)
res_list <- lapply(files, read.csv)

LOW_STRATA  <- c(Amyloid = "CERAD_stage2", Tau = "Braak_Limbic_III_IV", `TDP-43` = "TDP43_Early_1_2")
HIGH_STRATA <- c(Amyloid = "CERAD_stage3", Tau = "Braak_Neocortical_V", `TDP-43` = "TDP43_Late_3")

get_sig_set <- function(stratum_name, direction) {
  df <- res_list[[stratum_name]] %>%
    filter(!grepl("^X\\s*-\\s*\\d+", metabolite))
  if (direction == "Up") df %>% filter(logFC > 0, P.Value < 0.05) %>% pull(metabolite)
  else                   df %>% filter(logFC < 0, P.Value < 0.05) %>% pull(metabolite)
}

build_venn_list <- function(strata_map, direction) {
  setNames(lapply(strata_map, get_sig_set, direction = direction), names(strata_map))
}

low_up    <- build_venn_list(LOW_STRATA,  "Up")
low_down  <- build_venn_list(LOW_STRATA,  "Down")
high_up   <- build_venn_list(HIGH_STRATA, "Up")
high_down <- build_venn_list(HIGH_STRATA, "Down")

cat("Low  / Up   set sizes:  "); print(sapply(low_up, length))
cat("Low  / Down set sizes:  "); print(sapply(low_down, length))
cat("High / Up   set sizes:  "); print(sapply(high_up, length))
cat("High / Down set sizes:  "); print(sapply(high_down, length))

path_colors <- c("#56B4E9", "#E69F00", "#009E73")

make_venn <- function(venn_list, title_text) {
  ggvenn(
    venn_list,
    fill_color      = path_colors,
    fill_alpha      = 0.40,
    stroke_color    = "white",
    stroke_size     = 0.8,
    set_name_size   = 6,
    text_size       = 6,
    show_percentage = TRUE
  ) +
    labs(title = title_text) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

p_low_up    <- make_venn(low_up,    "Upregulated (higher in AsymAD) \u2014 Low pathology stage")
p_low_down  <- make_venn(low_down,  "Downregulated (higher in AD) \u2014 Low pathology stage")
p_high_up   <- make_venn(high_up,   "Upregulated (higher in AsymAD) \u2014 High pathology stage")
p_high_down <- make_venn(high_down, "Downregulated (higher in AD) \u2014 High pathology stage")

dir.create("figures", showWarnings = FALSE)
ggsave("figures/venn_metab_low_upregulated_johnson.png",    p_low_up,    width = 7, height = 7, dpi = 600, bg = "white")
ggsave("figures/venn_metab_low_downregulated_johnson.png",  p_low_down,  width = 7, height = 7, dpi = 600, bg = "white")
ggsave("figures/venn_metab_high_upregulated_johnson.png",   p_high_up,   width = 7, height = 7, dpi = 600, bg = "white")
ggsave("figures/venn_metab_high_downregulated_johnson.png", p_high_down, width = 7, height = 7, dpi = 600, bg = "white")

combined <- (p_low_up | p_low_down) / (p_high_up | p_high_down) +
  plot_annotation(
    title = "DE metabolite overlap across Amyloid / Tau / TDP-43 (raw p<0.05, Johnson rubric)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )
combined
ggsave("figures/venn_metab_all_four_combined_johnson.png", combined, width = 14, height = 14, dpi = 600, bg = "white")

triple_hits <- list(
  Low_Up    = Reduce(intersect, low_up),
  Low_Down  = Reduce(intersect, low_down),
  High_Up   = Reduce(intersect, high_up),
  High_Down = Reduce(intersect, high_down)
)
cat("\nTriple-hit metabolite set sizes:\n"); print(sapply(triple_hits, length)); print(triple_hits)
for (nm in names(triple_hits)) writeLines(triple_hits[[nm]], file.path(RES_DIR, paste0("metab_johnson_", nm, "_triple_pathology.txt")))

m <- make_comb_mat(triple_hits)
comb_order <- order(-comb_size(m))
cat("\nIntersection sizes (largest to smallest):\n")
print(data.frame(combination = comb_name(m)[comb_order], size = comb_size(m)[comb_order], row.names = NULL))

ht_upset <- UpSet(m,
                  comb_order = comb_order,
                  top_annotation = upset_top_annotation(
                    m, comb_order = comb_order, add_numbers = TRUE,
                    numbers_gp = gpar(fontsize = 10, fontface = "bold"), numbers_rot = 0,
                    gp = gpar(fill = "steelblue"), height = unit(4, "cm"),
                    annotation_name_gp = gpar(fontface = "bold")),
                  right_annotation = upset_right_annotation(
                    m, add_numbers = TRUE,
                    numbers_gp = gpar(fontsize = 9, fontface = "bold"),
                    gp = gpar(fill = "#E5789E"), width = unit(3, "cm"),
                    annotation_name_gp = gpar(fontface = "bold")),
                  row_names_gp = gpar(fontface = "bold", fontsize = 11),
                  pt_size = unit(4, "mm"), lwd = 2)

png("figures/upset_metab_triple_hits_johnson.png", width = 12, height = 6.5, units = "in", res = 300)
draw(ht_upset)
dev.off()
draw(ht_upset)

robust_up       <- intersect(triple_hits$Low_Up,   triple_hits$High_Up)
robust_down     <- intersect(triple_hits$Low_Down, triple_hits$High_Down)
flip_up_to_down <- intersect(triple_hits$Low_Up,   triple_hits$High_Down)
flip_down_to_up <- intersect(triple_hits$Low_Down, triple_hits$High_Up)

cat("\nRobust UP:", length(robust_up), "\n");     print(robust_up)
cat("Robust DOWN:", length(robust_down), "\n");   print(robust_down)
cat("Flip Up->Down:", length(flip_up_to_down), "\n"); print(flip_up_to_down)
cat("Flip Down->Up:", length(flip_down_to_up), "\n"); print(flip_down_to_up)

writeLines(robust_up,       file.path(RES_DIR, "metab_johnson_robust_up_both_severities.txt"))
writeLines(robust_down,     file.path(RES_DIR, "metab_johnson_robust_down_both_severities.txt"))
writeLines(flip_up_to_down, file.path(RES_DIR, "metab_johnson_flip_up_at_low_down_at_high.txt"))
writeLines(flip_down_to_up, file.path(RES_DIR, "metab_johnson_flip_down_at_low_up_at_high.txt"))

# ── 6. UpSet across all 6 strata, by direction ─────────────────
strata_names <- names(files)
up_sets   <- setNames(lapply(strata_names, get_sig_set, direction = "Up"),   strata_names)
down_sets <- setNames(lapply(strata_names, get_sig_set, direction = "Down"), strata_names)

cat("UP set sizes:\n");   print(sapply(up_sets, length))
cat("\nDOWN set sizes:\n"); print(sapply(down_sets, length))

build_upset <- function(set_list, bar_color) {
  m <- make_comb_mat(set_list)
  comb_order <- order(-comb_size(m))
  UpSet(m,
        comb_order = comb_order,
        top_annotation = upset_top_annotation(
          m, comb_order = comb_order, add_numbers = TRUE,
          numbers_gp = gpar(fontsize = 9, fontface = "bold"), numbers_rot = 0,
          gp = gpar(fill = bar_color), height = unit(4, "cm"),
          annotation_name_gp = gpar(fontface = "bold")),
        right_annotation = upset_right_annotation(
          m, add_numbers = TRUE,
          numbers_gp = gpar(fontsize = 8, fontface = "bold"),
          gp = gpar(fill = "#E5789E"), width = unit(3, "cm"),
          annotation_name_gp = gpar(fontface = "bold")),
        row_names_gp = gpar(fontface = "bold", fontsize = 10),
        pt_size = unit(3.5, "mm"), lwd = 1.8)
}

ht_up <- build_upset(up_sets, "firebrick")
png("figures/upset_metab_up_all_strata_johnson.png", width = 16, height = 7, units = "in", res = 300)
draw(ht_up, column_title = "DE metabolites UP (higher in AsymAD), raw p<0.05, across all 6 strata (Johnson)",
     column_title_gp = gpar(fontface = "bold", fontsize = 13))
dev.off()
draw(ht_up, column_title = "DE metabolites UP (higher in AsymAD), raw p<0.05, across all 6 strata (Johnson)",
     column_title_gp = gpar(fontface = "bold", fontsize = 13))

ht_down <- build_upset(down_sets, "steelblue")
png("figures/upset_metab_down_all_strata_johnson.png", width = 16, height = 7, units = "in", res = 300)
draw(ht_down, column_title = "DE metabolites DOWN (higher in AD), raw p<0.05, across all 6 strata (Johnson)",
     column_title_gp = gpar(fontface = "bold", fontsize = 13))
dev.off()
draw(ht_down, column_title = "DE metabolites DOWN (higher in AD), raw p<0.05, across all 6 strata (Johnson)",
     column_title_gp = gpar(fontface = "bold", fontsize = 13))

# ── 7. logFC heatmap for metabolites significant in MIN_STRATA+ strata ──────────────
MIN_STRATA <- 3

stratum_info <- tibble::tribble(
  ~file,                                                ~stratum,              ~pathology, ~severity,
  file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_CERAD_stage2.csv"),          "CERAD_stage2",        "Amyloid",  "Low",
  file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_CERAD_stage3.csv"),          "CERAD_stage3",        "Amyloid",  "High",
  file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_Braak_Limbic_III_IV.csv"),   "Braak_Limbic_III_IV", "Tau",      "Low",
  file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_Braak_Neocortical_V.csv"),   "Braak_Neocortical_V", "Tau",      "High",
  file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_TDP43_Early_1_2.csv"),       "TDP43_Early_1_2",     "TDP-43",   "Low",
  file.path(RES_DIR, "limma_metab_AsymAD_vs_AD_TDP43_Late_3.csv"),          "TDP43_Late_3",        "TDP-43",   "High"
)
res_list2 <- setNames(lapply(stratum_info$file, read.csv), stratum_info$stratum)

get_sig_set2 <- function(stratum_name, direction) {
  df <- res_list2[[stratum_name]] %>% filter(!grepl("^X\\s*-\\s*\\d+", metabolite))
  if (direction == "Up") df %>% filter(logFC > 0, P.Value < 0.05) %>% pull(metabolite)
  else                   df %>% filter(logFC < 0, P.Value < 0.05) %>% pull(metabolite)
}

up_sets2   <- setNames(lapply(stratum_info$stratum, get_sig_set2, direction = "Up"),   stratum_info$stratum)
down_sets2 <- setNames(lapply(stratum_info$stratum, get_sig_set2, direction = "Down"), stratum_info$stratum)

membership_counts <- function(set_list) {
  m <- make_comb_mat(set_list)
  bind_rows(lapply(comb_name(m), function(cn) {
    members <- extract_comb(m, cn)
    if (length(members) == 0) return(NULL)
    n_strata <- sum(as.integer(strsplit(cn, "")[[1]]))
    data.frame(metabolite = members, n_strata = n_strata)
  }))
}

up_multi   <- membership_counts(up_sets2)   %>% filter(n_strata >= MIN_STRATA) %>% pull(metabolite) %>% unique()
down_multi <- membership_counts(down_sets2) %>% filter(n_strata >= MIN_STRATA) %>% pull(metabolite) %>% unique()

cat("Up metabolites (", MIN_STRATA, "+ strata):", length(up_multi), "\n")
cat("Down metabolites (", MIN_STRATA, "+ strata):", length(down_multi), "\n")

if (length(up_multi) == 0 && length(down_multi) == 0) {
  stop("No metabolites meet MIN_STRATA = ", MIN_STRATA, " in either direction — lower MIN_STRATA and retry.")
}

metab_order <- c(up_multi, down_multi)

heat_data <- stratum_info %>%
  mutate(data = lapply(stratum, function(s) {
    res_list2[[s]] %>% filter(metabolite %in% metab_order) %>%
      select(metabolite, logFC, P.Value, adj.P.Val)
  })) %>%
  select(-file) %>%
  tidyr::unnest(data)

heat_data <- heat_data %>%
  mutate(
    metabolite = factor(metabolite, levels = rev(metab_order)),
    pathology  = factor(pathology, levels = c("Amyloid", "Tau", "TDP-43")),
    severity   = factor(severity, levels = c("Low", "High")),
    stars = case_when(P.Value < 0.001 ~ "***", P.Value < 0.01 ~ "**", P.Value < 0.05 ~ "*", TRUE ~ ""),
    label = stars,
    x_label = factor(paste0(pathology, "\n", severity),
                     levels = paste0(rep(c("Amyloid", "Tau", "TDP-43"), each = 2),
                                     "\n", rep(c("Low", "High"), 3)))
  )

divider_y <- length(down_multi) + 0.5

p_heat <- ggplot(heat_data, aes(x = x_label, y = metabolite, fill = logFC)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 6, fontface = "bold", lineheight = 0.85) +
  geom_hline(yintercept = divider_y, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "logFC\n(+ = higher in AsymAD)") +
  labs(x = NULL, y = NULL,
       title = paste0("Metabolites significant in ", MIN_STRATA, "+ strata (same direction, raw p<0.05, Johnson rubric)"),
       subtitle = "Top: consistently UP in AsymAD  |  Bottom: consistently UP in AD",
       caption = "Cell label: logFC (this stratum's own limma fit) and significance stars, raw p-value.\n* p<0.05  ** p<0.01  *** p<0.001") +
  theme_minimal(base_size = 12) +
  theme(axis.text     = element_text(face = "bold", color = "black"),
        axis.text.y   = element_text(size = 12, face = "bold"),
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 10),
        plot.caption  = element_text(hjust = 0, color = "black", size = 8),
        legend.title  = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid    = element_blank())

p_heat
ggsave("figures/heatmap_metab_multistratum_logFC_johnson_1.png", p_heat,
       width = 9, height = 12, dpi = 600, bg = "white")




# ── 8. flipped layout of the same heatmap ────────────────────────────────────
# Metabolites on the x-axis (UP in AsymAD on the left, UP in AD on the right), the six strata on the
# y-axis (Amyloid Low at the top). FILL = adjusted AsymAD - AD mean difference in SD units: the
# metabolite data are z-scored per metabolite, so limma's "logFC" is an SD difference, not a log2
# fold change.
# ──────────────────────────────────────────────────────────────────────────────
y_levels <- paste0(rep(c("Amyloid", "Tau", "TDP-43"), each = 2),
                   "\n", rep(c("Low", "High"), 3))

heat_data <- heat_data %>%
  mutate(
    metabolite  = factor(metabolite, levels = metab_order),   # UP set first (left), DOWN set after
    pathology   = factor(pathology, levels = c("Amyloid", "Tau", "TDP-43")),
    severity    = factor(severity,  levels = c("Low", "High")),
    stars       = case_when(P.Value < 0.001 ~ "***", P.Value < 0.01 ~ "**",
                            P.Value < 0.05  ~ "*",   TRUE ~ ""),
    stratum_lab = factor(paste0(pathology, "\n", severity),
                         levels = rev(y_levels))              # rev: first level is drawn at the bottom
  )

divider_x <- length(up_multi) + 0.5      # dashed line between the UP block and the DOWN block

p_heat <- ggplot(heat_data, aes(x = metabolite, y = stratum_lab, fill = logFC)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = stars), size = 6, fontface = "bold", lineheight = 0.85) +
  geom_vline(xintercept = divider_x, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Difference in metabolite level (SD)\n(+ = higher in AsymAD)") +
  scale_x_discrete(labels = function(x)          # ONE line break per name: before "(" if it has one, else after "/"
    ifelse(grepl(" (", x, fixed = TRUE),
           sub(" (", "\n(", x, fixed = TRUE),
           gsub("/", "/\n", x, fixed = TRUE))) +
  labs(x = NULL, y = NULL,
       title = paste0("Metabolites significant in ", MIN_STRATA,
                      "+ strata (same direction, raw p<0.05, Johnson rubric)"),
       subtitle = "Left: consistently UP in AsymAD  |  Right: consistently UP in AD",
       caption = paste0("Fill: AsymAD - AD mean difference in metabolite level within the stratum (SD units; ",
                        "metabolite data are z-scored per metabolite), from each stratum's own limma fit.\n",
                        "Stars: raw p-value. * p<0.05  ** p<0.01  *** p<0.001")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x     = element_text(face = "bold", color = "black", size = 12.5,
                                       angle = 60, hjust = 1, vjust = 1, lineheight = 0.85),
        axis.text.y     = element_text(face = "bold", color = "black", size = 12),
        plot.title      = element_text(face = "bold", hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey30", size = 10),
        plot.caption    = element_text(hjust = 0, color = "black", size = 8),
        legend.title    = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid      = element_blank(),
        plot.margin     = margin(t = 5.5, r = 5.5, b = 5.5, l = 45))   # room for the first rotated label

p_heat

n_m <- length(metab_order)
ggsave("figures/heatmap_metab_multistratum_SD_johnson_flipped.png", p_heat,
       width = max(8, 0.33 * n_m + 3), height = 7.5, dpi = 600, bg = "white")
