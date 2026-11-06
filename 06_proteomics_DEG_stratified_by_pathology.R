# ══════════════════════════════════════════════════════════════════════════
# Proteomics differential abundance, AsymAD vs AD, within pathology strata

#
#   Six independent limma models (diagnosis + age, sex, APOE risk, PMI) within CERAD stage 2 / 3,
#   Braak III+IV / V and TDP-43 stage 1+2 / 3; stacked count plot; Venn diagrams and robust (same-
#   direction) protein lists; logFC heatmap of the robust proteins; proteins unique to one stratum;
#   Reactome dot plots of the unique-protein lists.
#
# INPUTS
#   R objects: meta_cc_8reg_n323_2026_09_02, expr_filtered (see data/README.md)
#   data/reactome/Reactome_Pathways_2024_table-7.txt ... table-18.txt (Enrichr exports, see data/reactome/README.md)
# OUTPUTS
#   results/limma_AsymAD_vs_AD_<stratum>.csv, results/limma_AsymAD_vs_AD_all_strata_combined.csv
#   results/*_triple_pathology_proteins.txt, results/unique_protein/*
#   figures/deg_*, venn_*, heatmap_robust_*, unique_pathway_dotplots/
# RUN AFTER
#   none
#
# Run from the project root (see README.md). The analysis scripts keep their
# intermediate objects in the R session; run_all.R sources them in order.
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")
require_packages(c("dplyr", "tidyr", "purrr", "limma", "tibble", "ggplot2", "ggvenn", "patchwork", "readr", "stringr"))
load_inputs(c("meta_cc_8reg_n323_2026_09_02", "expr_filtered"))

library(dplyr)
library(tidyr)
library(purrr)
library(limma)
library(tibble)
library(ggplot2)

# ── 0. inputs: protein expression matrix + metadata ─────────────────────────
protein_expr <- expr_filtered %>% column_to_rownames(var = "...1")
colnames(protein_expr) <- as.character(colnames(protein_expr))

DX_LEVELS <- c("AD", "AsymAD")   # AD = reference, always

# Minimum fraction of non-missing values required in AsymAD or AD. It is 0, i.e. no protein is dropped
# here (proteins were filtered upstream; see README, "Presence filters").
MIN_PCT_PRESENT <- 0.00

# ── 1. one function: build pheno for a stratum, presence filter, fit limma ──────
run_deg_stratum <- function(stratum_name, stage_col, stage_values, min_pct = MIN_PCT_PRESENT) {
  pheno <- meta_cc_8reg_n323_2026_09_02 %>%
    filter(EmoryStrictDx.2019 %in% c("AD", "AsymAD"),
           .data[[stage_col]] %in% stage_values) %>%
    transmute(projid    = as.character(...1),
              diagnosis = factor(EmoryStrictDx.2019, levels = DX_LEVELS),
              age_death, msex, apoe_risk, pmi) %>%
    filter(projid %in% colnames(protein_expr)) %>%
    mutate(age_death_z = as.numeric(scale(age_death)),
           pmi_z        = as.numeric(scale(pmi)),
           msex         = factor(msex),
           apoe_risk    = factor(apoe_risk))
  
  cat("\n==", stratum_name, "==\n")
  cat("n =", nrow(pheno), "\n")
  print(table(pheno$diagnosis))
  
  if (nrow(pheno) < 15 || n_distinct(pheno$diagnosis) < 2) {
    cat("SKIPPED:", stratum_name, "— too few samples or only one diagnosis present.\n")
    return(NULL)
  }
  
  expr <- protein_expr[, pheno$projid, drop = FALSE]
  stopifnot(identical(colnames(expr), pheno$projid))
  
  # presence filter (min_pct), per group, union
  asymad_ids <- pheno$projid[pheno$diagnosis == "AsymAD"]
  ad_ids     <- pheno$projid[pheno$diagnosis == "AD"]
  pct_present <- function(x) mean(!is.na(x))
  pass_asymad <- apply(expr[, asymad_ids, drop = FALSE], 1, pct_present) >= min_pct
  pass_ad     <- apply(expr[, ad_ids,     drop = FALSE], 1, pct_present) >= min_pct
  keep_proteins <- rownames(expr)[pass_asymad | pass_ad]
  expr_filt <- expr[keep_proteins, , drop = FALSE]
  cat("Proteins kept (presence >=", min_pct, ", union):", length(keep_proteins), "/", nrow(expr), "\n")
  
  design <- model.matrix(~ diagnosis + age_death_z + msex + apoe_risk + pmi_z, data = pheno)
  
  fit <- lmFit(expr_filt, design)
  fit <- eBayes(fit)
  
  res <- topTable(fit, coef = "diagnosisAsymAD", number = Inf, adjust.method = "BH") %>%
    tibble::rownames_to_column("protein") %>%
    arrange(adj.P.Val) %>%
    mutate(stratum = stratum_name)
  
  fname <- file.path(RES_DIR, paste0("limma_AsymAD_vs_AD_", stratum_name, ".csv"))
  write.csv(res %>% select(-stratum), fname, row.names = FALSE)
  cat("Saved:", fname, "\n")
  cat("Significant (raw p<0.05):", sum(res$P.Value < 0.05, na.rm = TRUE),
      " | (BH adj.P<0.05):", sum(res$adj.P.Val < 0.05, na.rm = TRUE), "\n")
  
  res
}

# ── 2. the six strata ────────────────────────────────────────────────────────
strata_specs <- list(
  list(name = "CERAD_stage2",         stage_col = "CERAD_std", stage_values = 2),
  list(name = "CERAD_stage3",         stage_col = "CERAD_std", stage_values = 3),
  list(name = "Braak_Limbic_III_IV",  stage_col = "braaksc",   stage_values = c(3, 4)),
  list(name = "Braak_Neocortical_V",  stage_col = "braaksc",   stage_values = 5),
  list(name = "TDP43_Early_1_2",      stage_col = "tdp_st4",   stage_values = c(1, 2)),
  list(name = "TDP43_Late_3",         stage_col = "tdp_st4",   stage_values = 3)
)

all_res <- map_dfr(strata_specs, function(s)
  run_deg_stratum(s$name, s$stage_col, s$stage_values))

write.csv(all_res, file.path(RES_DIR, "limma_AsymAD_vs_AD_all_strata_combined.csv"), row.names = FALSE)

# ── 3. stacked barplot: Up/Down counts per stratum, at BOTH thresholds ──────
# raw-p and BH-adjusted counts are both shown
summary_counts <- all_res %>%
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
  labs(x = NULL, y = "Number of DE proteins",
       title = "DE protein counts across pathology strata (AsymAD vs AD, AD = reference)") +
  theme_bw(base_size = 13) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y     = element_text(face = "bold", color = "black"),
        axis.title      = element_text(face = "bold", color = "black"),
        strip.text      = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "grey92", color = NA),
        plot.title      = element_text(face = "bold", hjust = 0.5),
        legend.title    = element_text(face = "bold"),
        legend.text     = element_text(face = "bold"))

p
dir.create("figures", showWarnings = FALSE)
ggsave("figures/deg_counts_stacked_barplot.png", p, width = 10, height = 9, dpi = 600, bg = "white")


# ══════════════════════════════════════════════════════════════════════════
# 4 Venn diagrams: Amyloid vs Tau vs TDP-43 overlap of DE proteins
# (AsymAD vs AD, AD = reference, raw p<0.05), split by:
#   - severity: Low (CERAD stage2 / Braak Limbic III-IV / TDP-43 Early 1+2)
#               vs High (CERAD stage3 / Braak Neocortical V / TDP-43 Late 3)
#   - direction: Up (higher in AsymAD) vs Down (higher in AD)
#
# Requires the six result CSVs written above (results/limma_AsymAD_vs_AD_<stratum>.csv)
# ══════════════════════════════════════════════════════════════════════════

library(dplyr)
library(ggvenn)
library(ggplot2)
library(patchwork)

# ── 0. load the 6 stratified DEG results ────────────────────────────────────
files <- list(
  CERAD_stage2        = file.path(RES_DIR, "limma_AsymAD_vs_AD_CERAD_stage2.csv"),
  CERAD_stage3        = file.path(RES_DIR, "limma_AsymAD_vs_AD_CERAD_stage3.csv"),
  Braak_Limbic_III_IV = file.path(RES_DIR, "limma_AsymAD_vs_AD_Braak_Limbic_III_IV.csv"),
  Braak_Neocortical_V = file.path(RES_DIR, "limma_AsymAD_vs_AD_Braak_Neocortical_V.csv"),
  TDP43_Early_1_2     = file.path(RES_DIR, "limma_AsymAD_vs_AD_TDP43_Early_1_2.csv"),
  TDP43_Late_3        = file.path(RES_DIR, "limma_AsymAD_vs_AD_TDP43_Late_3.csv")
)
res_list <- lapply(files, read.csv)

# severity groupings — mapped to Amyloid / Tau / TDP-43
LOW_STRATA  <- c(Amyloid = "CERAD_stage2", Tau = "Braak_Limbic_III_IV", `TDP-43` = "TDP43_Early_1_2")
HIGH_STRATA <- c(Amyloid = "CERAD_stage3", Tau = "Braak_Neocortical_V", `TDP-43` = "TDP43_Late_3")

# ── 1. pull the significant protein set (raw p<0.05) for one stratum + direction
get_sig_set <- function(stratum_name, direction) {
  df <- res_list[[stratum_name]]
  if (direction == "Up") {
    df %>% filter(logFC > 0, P.Value < 0.05) %>% pull(protein)
  } else {
    df %>% filter(logFC < 0, P.Value < 0.05) %>% pull(protein)
  }
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

# ── 2. Okabe-Ito colorblind-safe palette ─
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

ggsave("figures/venn_low_upregulated.png",    p_low_up,    width = 7, height = 7, dpi = 600, bg = "white")
ggsave("figures/venn_low_downregulated.png",  p_low_down,  width = 7, height = 7, dpi = 600, bg = "white")
ggsave("figures/venn_high_upregulated.png",   p_high_up,   width = 7, height = 7, dpi = 600, bg = "white")
ggsave("figures/venn_high_downregulated.png", p_high_down, width = 7, height = 7, dpi = 600, bg = "white")

# ── 3. combined 2x2 for a quick overview (Low/High x Up/Down) ───────────────
combined <- (p_low_up | p_low_down) / (p_high_up | p_high_down) +
  plot_annotation(
    title = "DE protein overlap across Amyloid / Tau / TDP-43 (raw p<0.05, AsymAD vs AD)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

combined
ggsave("figures/venn_all_four_combined.png", combined, width = 14, height = 14, dpi = 600, bg = "white")

# ── 4. proteins consistently up/down across ALL THREE pathologies at a given
#      severity — the intersection cell of each Venn, i.e. the
#      actual gene lists behind the diagram
triple_hits <- list(
  Low_Up    = Reduce(intersect, low_up),
  Low_Down  = Reduce(intersect, low_down),
  High_Up   = Reduce(intersect, high_up),
  High_Down = Reduce(intersect, high_down)
)
print(triple_hits)

for (nm in names(triple_hits)) {
  writeLines(triple_hits[[nm]], file.path(RES_DIR, paste0(nm, "_triple_pathology_proteins.txt")))
}






# ══════════════════════════════════════════════════════════════════════════
# logFC heatmap for the 22 "robust" proteins (consistent direction at both
# Low and High pathology severity, across Amyloid/Tau/TDP-43)
#
#   Rows: the 19 robust-UP + 3 robust-DOWN proteins
#   Columns: all 6 strata, grouped Amyloid | Tau | TDP-43, Low then High
#   Fill: logFC (from each stratum's own limma fit)
#   Label: logFC value + significance stars (raw p < 0.05)
#
# This shows the actual claim directly — same-direction logFC across a whole
# row — rather than forcing raw expression into an artificial per-donor
# Low/High bucket (which doesn't exist cleanly, since a donor's severity
# differs by WHICH pathology you're asking about).
# ══════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)
library(ggplot2)

# ── 0. the 22 robust proteins ────────────────────────────────────────────────
# Robust proteins = same direction at BOTH severities for all three pathologies (intersection of the
# two triple-pathology sets). These objects were not created anywhere in the original script (they were
# defined elsewhere in the author's session); they are computed here in the same way as for metabolites
# in script 09.
robust_up   <- intersect(triple_hits$Low_Up,   triple_hits$High_Up)
robust_down <- intersect(triple_hits$Low_Down, triple_hits$High_Down)
cat("Robust UP:", length(robust_up), "| Robust DOWN:", length(robust_down),
    "(the heatmap below was built for 19 up + 3 down = 22 proteins)\n")

genes_up   <- robust_up
genes_down <- robust_down
gene_order <- c(genes_up, genes_down)   # Up block first, Down block second

# ── 1. load the 6 stratified DEG results, tag pathology + severity ─────────
stratum_info <- tibble::tribble(
  ~file,                                          ~stratum,              ~pathology, ~severity,
  file.path(RES_DIR, "limma_AsymAD_vs_AD_CERAD_stage2.csv"),          "CERAD_stage2",        "Amyloid",  "Low",
  file.path(RES_DIR, "limma_AsymAD_vs_AD_CERAD_stage3.csv"),          "CERAD_stage3",        "Amyloid",  "High",
  file.path(RES_DIR, "limma_AsymAD_vs_AD_Braak_Limbic_III_IV.csv"),   "Braak_Limbic_III_IV", "Tau",      "Low",
  file.path(RES_DIR, "limma_AsymAD_vs_AD_Braak_Neocortical_V.csv"),   "Braak_Neocortical_V", "Tau",      "High",
  file.path(RES_DIR, "limma_AsymAD_vs_AD_TDP43_Early_1_2.csv"),       "TDP43_Early_1_2",     "TDP-43",   "Low",
  file.path(RES_DIR, "limma_AsymAD_vs_AD_TDP43_Late_3.csv"),          "TDP43_Late_3",        "TDP-43",   "High"
)

heat_data <- stratum_info %>%
  mutate(data = lapply(file, function(f) {
    read.csv(f) %>% filter(protein %in% gene_order) %>%
      select(protein, logFC, P.Value, adj.P.Val)
  })) %>%
  select(-file) %>%
  tidyr::unnest(data)

# ── 2. tidy up for plotting ──────────────────────────────────────────────────
heat_data <- heat_data %>%
  mutate(
    protein  = factor(protein, levels = rev(gene_order)),   # Up block on top, Down on bottom
    pathology = factor(pathology, levels = c("Amyloid", "Tau", "TDP-43")),
    severity  = factor(severity, levels = c("Low", "High")),
    stars = case_when(P.Value < 0.001 ~ "***",
                      P.Value < 0.01  ~ "**",
                      P.Value < 0.05  ~ "*",
                      TRUE ~ ""),
    label = stars)

# x-axis: pathology (outer grouping) x severity (inner), in a fixed order
heat_data <- heat_data %>%
  mutate(x_label = factor(paste0(pathology, "\n", severity),
                          levels = paste0(rep(c("Amyloid", "Tau", "TDP-43"), each = 2),
                                          "\n", rep(c("Low", "High"), 3))))

# ── 3. plot ───────────────────────────────────────────────────────────────
divider_y <- length(genes_down) + 0.5   # boundary between Down block (bottom) and Up block (top)

p <- ggplot(heat_data, aes(x = x_label, y = protein, fill = logFC)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 6, fontface = "bold", lineheight = 0.85) +
  geom_hline(yintercept = divider_y, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "logFC\n(+ = higher in AsymAD)") +
  labs(x = NULL, y = NULL,
       title = "Robust protein\u2013pathology relationships across severity",
       subtitle = "Top: consistently UP in AsymAD  |  Bottom: consistently UP in AD  \u2014  both Low and High severity",
       caption = "Cell label: logFC (this stratum's own limma fit) and significance stars, raw p-value.\n* p<0.05  ** p<0.01  *** p<0.001") +
  theme_minimal(base_size = 12) +
  theme(axis.text     = element_text(face = "bold", color = "black"),
        axis.text.y   = element_text(size = 12),
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 10),
        plot.caption  = element_text(hjust = 0, color = "black", size = 8),
        legend.title  = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid    = element_blank())

p
dir.create("figures", showWarnings = FALSE)
ggsave("figures/heatmap_robust_proteins_logFC.png", p, width = 8, height = 12, dpi = 600, bg = "white")
ggsave("figures/heatmap_robust_proteins_logFC_1.png", p, dpi = 600, bg = "white")











# ══════════════════════════════════════════════════════════════════════════
# Proteins UNIQUE to each stratum — significant (raw p<0.05) in exactly ONE
# of the 6 strata and not shared with any of the other 5 (same direction).
#
#   6 strata x 2 directions (Up/Down) = 12 sets -> 12 .txt files, saved into
#   results/unique_protein/, plus a stacked bar plot of the counts.
# ══════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tibble)
library(ggplot2)

dir.create(file.path(RES_DIR, "unique_protein"), showWarnings = FALSE)

# ── 0. load the 6 stratified protein DEG results ────────────────────────────
files <- list(
  CERAD_stage2        = file.path(RES_DIR, "limma_AsymAD_vs_AD_CERAD_stage2.csv"),
  CERAD_stage3        = file.path(RES_DIR, "limma_AsymAD_vs_AD_CERAD_stage3.csv"),
  Braak_Limbic_III_IV = file.path(RES_DIR, "limma_AsymAD_vs_AD_Braak_Limbic_III_IV.csv"),
  Braak_Neocortical_V = file.path(RES_DIR, "limma_AsymAD_vs_AD_Braak_Neocortical_V.csv"),
  TDP43_Early_1_2     = file.path(RES_DIR, "limma_AsymAD_vs_AD_TDP43_Early_1_2.csv"),
  TDP43_Late_3        = file.path(RES_DIR, "limma_AsymAD_vs_AD_TDP43_Late_3.csv")
)
res_list     <- lapply(files, read.csv)
strata_names <- names(files)

# ── 1. one Up/Down set per stratum ──────────────────────────────────────────
get_sig_set <- function(stratum_name, direction) {
  df <- res_list[[stratum_name]]
  if (direction == "Up") df %>% filter(logFC > 0, P.Value < 0.05) %>% pull(protein)
  else                   df %>% filter(logFC < 0, P.Value < 0.05) %>% pull(protein)
}

up_sets   <- setNames(lapply(strata_names, get_sig_set, direction = "Up"),   strata_names)
down_sets <- setNames(lapply(strata_names, get_sig_set, direction = "Down"), strata_names)

# ── 2. unique-to-one-stratum: present here, absent from ALL 5 others ───────
get_unique <- function(set_list, stratum_name) {
  others <- Reduce(union, set_list[names(set_list) != stratum_name])
  setdiff(set_list[[stratum_name]], others)
}

unique_up   <- setNames(lapply(strata_names, get_unique, set_list = up_sets),   strata_names)
unique_down <- setNames(lapply(strata_names, get_unique, set_list = down_sets), strata_names)

cat("Unique UP counts per stratum:\n");   print(sapply(unique_up, length))
cat("\nUnique DOWN counts per stratum:\n"); print(sapply(unique_down, length))

# ── 3. save the 12 .txt files ───────────────────────────────────────────────
for (s in strata_names) {
  writeLines(unique_up[[s]],   file.path(RES_DIR, "unique_protein", paste0("unique_up_",   s, ".txt")))
  writeLines(unique_down[[s]], file.path(RES_DIR, "unique_protein", paste0("unique_down_", s, ".txt")))
}
cat("\nSaved 12 .txt files to unique_protein/\n")

# ── 4. stacked bar plot of unique counts ────────────────────────────────────
summary_counts <- bind_rows(
  tibble(stratum = strata_names, Direction = "Up",   n = sapply(unique_up,   length)),
  tibble(stratum = strata_names, Direction = "Down", n = sapply(unique_down, length))
) %>%
  mutate(stratum   = factor(stratum, levels = strata_names),
         Direction = factor(Direction, levels = c("Up", "Down")))

p <- ggplot(summary_counts, aes(x = stratum, y = n, fill = Direction)) +
  geom_col(position = "stack", width = 0.65) +
  geom_text(aes(label = ifelse(n > 0, n, "")),
            position = position_stack(vjust = 0.5),
            fontface = "bold", size = 6, color = "black") +
  scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue"), name = "Direction") +
  labs(x = NULL, y = "Number of UNIQUE DE proteins (not shared with any other stratum)",
       title = "Proteins uniquely DE in each pathology stratum (raw p<0.05)") +
  theme_bw(base_size = 20) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y      = element_text(face = "bold", color = "black"),
        axis.title       = element_text(face = "bold", color = "black"),
        plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.title     = element_text(face = "bold"),
        legend.text      = element_text(face = "bold"))

p
ggsave(file.path(RES_DIR, "unique_protein", "unique_protein_counts_stacked_barplot.png"), p,
       width = 10, height = 7, dpi = 600, bg = "transparent")




# pathways

# ══════════════════════════════════════════════════════════════════════════
# Pathway dot plots for the 12 "unique per stratum" protein lists —
# one Reactome enrichment result per stratum x direction, mapped by upload
# order. VERIFY this mapping against the original enrichment submissions before
# reusing it (see data/reactome/README.md).
#
#   x     = Gene ratio
#   size  = gene count
#   color = Adjusted P-value (raw, not -log10) — red gradient for Up lists,
#           blue gradient for Down lists
#   Filtered to Adjusted P-value < 0.1 first, then top N shown
# ══════════════════════════════════════════════════════════════════════════

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

TOP_N <- 10   # top pathways per file, by Adjusted P-value, after the 0.1 filter

# ── 0. the file <-> stratum/direction mapping (VERIFY before trusting output)
file_map <- tibble::tribble(
  ~file,                              ~stratum,              ~direction,
  "data/reactome/Reactome_Pathways_2024_table-7.txt",  "TDP43_Late_3",         "Up",
  "data/reactome/Reactome_Pathways_2024_table-8.txt",  "TDP43_Early_1_2",      "Up",
  "data/reactome/Reactome_Pathways_2024_table-9.txt",  "CERAD_stage3",         "Up",
  "data/reactome/Reactome_Pathways_2024_table-10.txt", "CERAD_stage2",         "Up",
  "data/reactome/Reactome_Pathways_2024_table-11.txt", "Braak_Neocortical_V",  "Up",
  "data/reactome/Reactome_Pathways_2024_table-12.txt", "Braak_Limbic_III_IV",  "Up",
  "data/reactome/Reactome_Pathways_2024_table-13.txt", "TDP43_Late_3",         "Down",
  "data/reactome/Reactome_Pathways_2024_table-14.txt", "TDP43_Early_1_2",      "Down",
  "data/reactome/Reactome_Pathways_2024_table-15.txt", "CERAD_stage3",         "Down",
  "data/reactome/Reactome_Pathways_2024_table-16.txt", "CERAD_stage2",         "Down",
  "data/reactome/Reactome_Pathways_2024_table-17.txt", "Braak_Neocortical_V",  "Down",
  "data/reactome/Reactome_Pathways_2024_table-18.txt", "Braak_Limbic_III_IV",  "Down"
)

print(file_map)

# ── 1. parse "Overlap" into gene count, pathway size, ratio ─────────────────
parse_overlap <- function(df) {
  df %>%
    mutate(n_genes        = as.numeric(str_extract(Overlap, "^[0-9]+")),
           n_pathway_total = as.numeric(str_extract(Overlap, "(?<=/)[0-9]+")),
           GeneRatio       = n_genes / n_pathway_total)
}

# ── 2. one plot builder, reused for all 12 ──────────────────────────────────
build_dotplot <- function(df, title_text, grad_low, grad_high) {
  ggplot(df, aes(x = GeneRatio, y = Term)) +
    geom_point(aes(size = n_genes, color = `Adjusted P-value`)) +
    scale_color_gradient(low = grad_low, high = grad_high, name = "Adj. P-value") +
    scale_size_continuous(name = "Gene count", range = c(2.5, 10)) +
    labs(x = "Gene ratio (overlap / pathway size)", y = NULL, title = title_text) +
    theme_bw(base_size = 13) +
    theme(axis.text.y       = element_text(face = "bold", color = "black", size = 9),
          axis.text.x       = element_text(face = "bold", color = "black"),
          axis.title        = element_text(face = "bold", color = "black"),
          plot.title        = element_text(face = "bold", hjust = 0.5, color = "black", size = 14),
          panel.grid.minor  = element_blank(),
          legend.title      = element_text(face = "bold", color = "black"),
          legend.text       = element_text(face = "bold", color = "black"))
}

dir.create("figures/unique_pathway_dotplots", recursive = TRUE, showWarnings = FALSE)

# ── 3. loop over all 12, filter to adj.p<0.1, plot + save ───────────────────
for (i in seq_len(nrow(file_map))) {
  fname     <- file_map$file[i]
  stratum   <- file_map$stratum[i]
  direction <- file_map$direction[i]
  
  df <- read_tsv(fname, show_col_types = FALSE) %>%
    parse_overlap() %>%
    filter(`Adjusted P-value` < 0.1)
  
  cat(stratum, "-", direction, ": ", nrow(df), "pathways pass adj.p<0.1\n")
  
  if (nrow(df) == 0) {
    cat("  SKIPPED (nothing passes the 0.1 filter)\n")
    next
  }
  
  df <- df %>%
    arrange(`Adjusted P-value`) %>%
    slice_head(n = TOP_N) %>%
    arrange(n_genes) %>%
    mutate(Term = factor(Term, levels = Term))   # largest dot at top
  
  is_up <- direction == "Up"
  grad_low  <- if (is_up) "#67000D" else "#08306B"
  grad_high <- if (is_up) "#FEE5D9" else "#DEEBF7"
  
  title_text <- paste0(stratum, " \u2014 ", direction, " (unique to this stratum)")
  
  p <- build_dotplot(df, title_text, grad_low, grad_high)
  print(p)
  
  out_name <- sprintf("figures/unique_pathway_dotplots/unique_%s_%s.png",
                      tolower(direction), stratum)
  ggsave(out_name, p, width = 10, height = max(4, 0.3 * nrow(df) + 2), dpi = 600, bg = "white")
  cat("  Saved:", out_name, "\n")
}
