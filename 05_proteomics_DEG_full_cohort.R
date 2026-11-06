# ══════════════════════════════════════════════════════════════════════════
# Proteomics differential abundance: AsymAD vs AD, full cohort (limma)
# (originally: proteomics_full.R)
#
#   limma on all quantified proteins, AD = reference, adjusted for age at death, sex, APOE risk
#   and postmortem interval; Emory diagnosis labels (EmoryStrictDx.2019); volcano plot.
#
# INPUTS
#   R objects: meta_cc_8reg_n323_2026_09_02, expr_filtered (proteins x donors; first column = protein id)
# OUTPUTS
#   results/limma_AsymAD_vs_AD_full_cohort.csv, results/proteins_tested_full_cohort.txt
#   results/upregulated_/downregulated_proteins_AsymAD_vs_AD_full_cohort.txt
#   figures/volcano_AsymAD_vs_AD_full_cohort.png
# RUN AFTER
#   none (script 10 reads the objects `res` and `protein_expr` created here)
#
# Run from the project root (see README.md). The analysis scripts keep their
# intermediate objects in the R session; run_all.R sources them in order.
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")
require_packages(c("dplyr", "tidyverse", "limma", "ggplot2", "ggrepel"))
load_inputs(c("meta_cc_8reg_n323_2026_09_02", "expr_filtered"))

library(dplyr)
library(tidyverse)
library(limma)
library(ggplot2)
library(ggrepel)

# =============================================================================
# 0. INPUTS — adjust names here if yours differ
# =============================================================================
# protein_expr: proteins (rows) x samples (columns), colnames = projid, values
#               may contain NA for missing/undetected proteins
# meta_cc_8reg_n323_2026_09_02: full metadata, using EmoryStrictDx.2019 (the
#               "complete", non-revised diagnosis)

pheno <- meta_cc_8reg_n323_2026_09_02 %>%
  filter(EmoryStrictDx.2019 %in% c("AD", "AsymAD")) %>%
  transmute(
    projid    = as.character(...1),   # matches expr_filtered's column IDs
    diagnosis = factor(EmoryStrictDx.2019, levels = c("AD", "AsymAD")),  # AD = reference
    age_death, msex, apoe_risk, pmi
  )

stopifnot(!anyNA(pheno$diagnosis))   # confirms every donor kept was AsymAD or AD
cat("n AsymAD:", sum(pheno$diagnosis == "AsymAD"),
    " | n AD:", sum(pheno$diagnosis == "AD"), "\n")

# Convert a column named "id" to row names
protein_expr <- expr_filtered %>%
  column_to_rownames(var = "...1")

colnames(protein_expr) <- as.character(colnames(protein_expr))
pheno <- pheno %>% filter(projid %in% colnames(protein_expr))
expr  <- protein_expr[, pheno$projid, drop = FALSE]   # enforce identical sample order

stopifnot(identical(colnames(expr), pheno$projid))
cat("Donors with matched expression data:", ncol(expr), "\n")

# =============================================================================
# 1. Presence filter (MIN_PCT_PRESENT, per group, then union) on the WHOLE protein list,
#    no subsetting to a pathology-defined donor set beforehand
# =============================================================================
asymad_ids <- pheno$projid[pheno$diagnosis == "AsymAD"]
ad_ids     <- pheno$projid[pheno$diagnosis == "AD"]

pct_present <- function(x) mean(!is.na(x))

# Minimum fraction of non-missing values required in AsymAD or AD. It is 0 here, i.e. the filter keeps
# every protein: proteins were already filtered upstream (see README, 'Presence filters').
MIN_PCT_PRESENT <- 0.00

pass_asymad <- apply(expr[, asymad_ids, drop = FALSE], 1, pct_present) >= MIN_PCT_PRESENT
pass_ad     <- apply(expr[, ad_ids,     drop = FALSE], 1, pct_present) >= MIN_PCT_PRESENT

keep_proteins <- rownames(expr)[pass_asymad | pass_ad]   # union

cat("Proteins passing the presence filter in AsymAD:", sum(pass_asymad), "\n")
cat("Proteins passing the presence filter in AD:    ", sum(pass_ad), "\n")
cat("Union (final):                 ", length(keep_proteins), "out of", nrow(expr), "\n")
writeLines(keep_proteins, file.path(RES_DIR, "proteins_tested_full_cohort.txt"))
expr_filt <- expr[keep_proteins, , drop = FALSE]

# =============================================================================
# 2. Covariates — z-score continuous ones, same convention as the DLPFC models
# =============================================================================
pheno <- pheno %>%
  mutate(age_death_z = as.numeric(scale(age_death)),
         pmi_z        = as.numeric(scale(pmi)),
         msex         = factor(msex),
         apoe_risk    = factor(apoe_risk))

# =============================================================================
# 3. limma — AsymAD vs AD, AD as reference, FULL cohort
# =============================================================================
design <- model.matrix(~ diagnosis + age_death_z + msex + apoe_risk + pmi_z, data = pheno)

fit <- lmFit(expr_filt, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "diagnosisAsymAD", number = Inf, adjust.method = "BH") %>%
  tibble::rownames_to_column("protein") %>%
  arrange(adj.P.Val)

cat("Significant proteins (BH < 0.05):", sum(res$adj.P.Val < 0.05), "\n")

write.csv(res, file.path(RES_DIR, "limma_AsymAD_vs_AD_full_cohort.csv"), row.names = FALSE)

# =============================================================================
# 4. Volcano plot
# =============================================================================
sig_full_vol <- res %>% filter(adj.P.Val < 0.05)

label_full <- unique(c(
  sig_full_vol %>% filter(logFC > 0) %>% slice_min(adj.P.Val, n = 5) %>% pull(protein),
  sig_full_vol %>% filter(logFC < 0) %>% slice_min(adj.P.Val, n = 5) %>% pull(protein),
  sig_full_vol %>% filter(logFC > 0) %>% slice_max(logFC,   n = 5) %>% pull(protein),
  sig_full_vol %>% filter(logFC < 0) %>% slice_min(logFC,   n = 5) %>% pull(protein)
))

volcano_full <- res %>%
  mutate(
    neg_log_p = -log10(adj.P.Val),
    direction = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Higher in AsymAD",
      adj.P.Val < 0.05 & logFC < 0 ~ "Higher in AD",
      TRUE ~ "ns"
    ),
    label = ifelse(protein %in% label_full, protein, NA)
  )

p_full <- ggplot(volcano_full, aes(x = logFC, y = neg_log_p, color = direction)) +
  geom_point(aes(size = ifelse(direction == "ns", 0.8, 1.2)), alpha = 0.6, show.legend = FALSE) +
  geom_point(data = filter(volcano_full, direction != "ns"), aes(size = 1.5), alpha = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey60", linewidth = 0.4) +
  geom_label_repel(aes(label = label), size = 5, fontface = "bold",
                   max.overlaps = 20, segment.color = NA,
                   fill = alpha("white", 0.8), show.legend = FALSE) +
  scale_color_manual(
    values = c("Higher in AsymAD" = "red4", "Higher in AD" = "blue4", "ns" = "grey75"),
    name = "Direction") +
  scale_size_identity() +
  labs(
    title = "AsymAD vs AD (full cohort)",
    x = "logFC (positive = higher in AsymAD)",
    y = expression(-log[10](adj.P.Val))
  ) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(face = "bold"))

p_full
ggsave(file.path(FIG_DIR, "volcano_AsymAD_vs_AD_full_cohort.png"), plot = p_full,
       dpi = 600, bg = "white", width = 12, height = 9)

up_proteins   <- res %>% filter(adj.P.Val < 0.05, logFC > 0) %>% pull(protein)
down_proteins <- res %>% filter(adj.P.Val < 0.05, logFC < 0) %>% pull(protein)

writeLines(up_proteins,   file.path(RES_DIR, "upregulated_proteins_AsymAD_vs_AD_full_cohort.txt"))
writeLines(down_proteins, file.path(RES_DIR, "downregulated_proteins_AsymAD_vs_AD_full_cohort.txt"))

cat("Upregulated (higher in AsymAD):", length(up_proteins), "\n")
cat("Downregulated (higher in AD):  ", length(down_proteins), "\n")
