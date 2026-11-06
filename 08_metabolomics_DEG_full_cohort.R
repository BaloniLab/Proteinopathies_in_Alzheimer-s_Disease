# ══════════════════════════════════════════════════════════════════════════
# Metabolomics differential abundance: AsymAD vs AD, full cohort (limma, Johnson rubric)

#
#   Cohort = Johnson-rubric AsymAD and AD donors (object `clin`, script 07) with covariates and
#   amyloid / tangle burden from dataset_1744; 80% presence filter per group (union); limma, AD =
#   reference, adjusted for age, sex, APOE genotype and PMI. Pathology burden is attached to the
#   metadata but NOT used as a covariate (it partly defines the groups). Venn diagram and volcano.
#
# INPUTS
#   R objects: clin (created by script 07), Expr_metabolomics (script 07 removes the X features),
#   dataset_1744_cross_sectional_08_30_2026
# OUTPUTS
#   results/metadata_metabolomics_johnson_with_pathology.csv
#   results/limma_metabolomics_AsymAD_vs_AD_johnson_pathology.csv
#   results/AsymAD_projid_/AD_projid_metabolomics_johnson.txt
#   figures/venn_asymad_ad_metabolomics_johnson_pathology.png, figures/volcano_AsymAD_vs_AD_metabolomics_johnson_pathology.png
# RUN AFTER
#   07_metabolomics_classification_and_JT.R (same R session)
#
# Run from the project root (see README.md). The analysis scripts keep their
# intermediate objects in the R session; run_all.R sources them in order.
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")
require_packages(c("dplyr", "ggvenn", "ggplot2", "ggrepel", "limma", "tibble"))
load_inputs(c("dataset_1744_cross_sectional_08_30_2026", "Expr_metabolomics"))

library(dplyr)
library(ggvenn)
library(ggplot2)
library(ggrepel)
library(limma)
library(tibble)

if ("...1" %in% names(Expr_metabolomics)) {
  Expr_metabolomics <- Expr_metabolomics %>% column_to_rownames(var = "...1")
}

# `clin` (Johnson-rubric classification) is created by script 07 and is required here.
stopifnot("object `clin` not found: run scripts/07_metabolomics_classification_and_JT.R first (same session)" = exists("clin"))

# Unnamed "X - #####" features are removed at the source (as in script 07). This is a no-op when script 07
# has just run in the same session; it makes the script give the same result when run on its own.
x_cols <- grep("^X\\s*-\\s*\\d+", colnames(Expr_metabolomics), value = TRUE)
Expr_metabolomics <- Expr_metabolomics[, !colnames(Expr_metabolomics) %in% x_cols, drop = FALSE]

# ── 0. cohort (Johnson rubric, clin$group) + covariates/pathology joined in
#      from dataset_1744_cross_sectional_08_30_2026 ─────────────────────────
path_source <- dataset_1744_cross_sectional_08_30_2026 %>%
  mutate(projid = as.character(projid)) %>%
  select(projid, age_death, apoe_genotype, msex, pmi,
         amylsqrt_est_8reg, tangsqrt_est_8reg)

full_meta <- clin %>%
  filter(group %in% c("AD", "AsymAD")) %>%
  mutate(projid = as.character(projid)) %>%
  transmute(projid, diagnosis = factor(group, levels = c("AD", "AsymAD"))) %>%   # AD = reference
  left_join(path_source, by = "projid")

cat("Matched to dataset_1744_cross_sectional_08_30_2026:",
    sum(!is.na(full_meta$age_death)), "/", nrow(full_meta), "\n")
cat("n AsymAD:", sum(full_meta$diagnosis == "AsymAD"),
    " | n AD:", sum(full_meta$diagnosis == "AD"), "\n")

# ── 0b. restrict to donors that actually matched dataset_1744 (326 of 344) —
#      everything downstream (Venn, DEG, saved metadata) now runs on just
#      this complete-case set, not the full 344 ──────────────────────────
n_before <- nrow(full_meta)
full_meta <- full_meta %>% filter(!is.na(age_death))
cat("\nRestricted to donors matched in dataset_1744:", nrow(full_meta), "/", n_before, "\n")
cat("n AsymAD:", sum(full_meta$diagnosis == "AsymAD"),
    " | n AD:", sum(full_meta$diagnosis == "AD"), "(after restriction)\n")

# ── 1. collect the three projid sets ────────────────────────────────────────
asymad_ids <- full_meta$projid[full_meta$diagnosis == "AsymAD"]
ad_ids     <- full_meta$projid[full_meta$diagnosis == "AD"]
metab_ids  <- rownames(Expr_metabolomics)

cat("n AsymAD:", length(asymad_ids), " | n AD:", length(ad_ids), "\n")

venn_list <- list(
  "AsymAD"       = asymad_ids,
  "AD"           = ad_ids,
  "Metabolomics" = metab_ids
)

# ── 2. Venn plot ─────────────────────────────────────────────────────────────
p_venn <- ggvenn(
  venn_list,
  fill_color      = c("#6A1B9A", "red4", "#0072B2"),
  fill_alpha      = 0.40,
  stroke_color    = "white",
  stroke_size     = 0.8,
  set_name_size   = 6,
  text_size       = 6,
  show_percentage = FALSE
) +
  theme(plot.background = element_rect(fill = "white", color = NA))

p_venn
ggsave(p_venn, filename = file.path(FIG_DIR, "venn_asymad_ad_metabolomics_johnson_pathology.png"),
       width = 8, height = 8, dpi = 600, bg = "white")

# ── 3. subset metabolomics expression to AsymAD + AD donors present in it ───
asymad_matched <- intersect(asymad_ids, metab_ids)
ad_matched     <- intersect(ad_ids, metab_ids)

cat("AsymAD in metabolomics:", length(asymad_matched), "out of", length(asymad_ids), "\n")
cat("AD in metabolomics:    ", length(ad_matched),     "out of", length(ad_ids), "\n")

keep_ids <- union(asymad_matched, ad_matched)
Expr_metabolomics_subset <- Expr_metabolomics[keep_ids, , drop = FALSE]

# ── 4. covariates + pathology, matched by projid ────────────────────────────
pheno <- full_meta %>%
  filter(projid %in% keep_ids) %>%
  mutate(age_death_z   = as.numeric(scale(age_death)),
         pmi_z          = as.numeric(scale(pmi)),
         msex           = factor(msex),
         apoe_genotype  = factor(apoe_genotype))

stopifnot(!anyNA(pheno$diagnosis))
print(table(pheno$diagnosis))

cat("\nAmyloid/tangle burden missingness among the matched cohort:\n")
cat("  missing amylsqrt_est_8reg:", sum(is.na(pheno$amylsqrt_est_8reg)), "/", nrow(pheno), "\n")
cat("  missing tangsqrt_est_8reg:", sum(is.na(pheno$tangsqrt_est_8reg)), "/", nrow(pheno), "\n")

write.csv(pheno, file.path(RES_DIR, "metadata_metabolomics_johnson_with_pathology.csv"), row.names = FALSE)

expr_metab <- t(Expr_metabolomics_subset)              # metabolites x samples
expr_metab <- expr_metab[, pheno$projid, drop = FALSE]
stopifnot(identical(colnames(expr_metab), pheno$projid))

# ── 5. 80% presence filter per group, union ─────────────────────────────────
asymad_ids_m <- pheno$projid[pheno$diagnosis == "AsymAD"]
ad_ids_m     <- pheno$projid[pheno$diagnosis == "AD"]

pct_present <- function(x) mean(!is.na(x))

pass_asymad <- apply(expr_metab[, asymad_ids_m, drop = FALSE], 1, pct_present) >= 0.80
pass_ad     <- apply(expr_metab[, ad_ids_m,     drop = FALSE], 1, pct_present) >= 0.80

keep_metabolites <- rownames(expr_metab)[pass_asymad | pass_ad]

cat("Metabolites passing 80% in AsymAD:", sum(pass_asymad), "\n")
cat("Metabolites passing 80% in AD:    ", sum(pass_ad), "\n")
cat("Union (final):                    ", length(keep_metabolites), "out of", nrow(expr_metab), "\n")

expr_metab_filt <- expr_metab[keep_metabolites, , drop = FALSE]

# ── 6. limma — diagnosis DEG (pathology NOT included as a covariate; see
#      note at the top of this script) ─────────────────────────────────────
design <- model.matrix(~ diagnosis + age_death_z + msex + apoe_genotype + pmi_z, data = pheno)
cat("Design matrix condition number (kappa):", round(kappa(design), 1), "\n")

fit <- lmFit(expr_metab_filt, design)
fit <- eBayes(fit)

res_metab <- topTable(fit, coef = "diagnosisAsymAD", number = Inf, adjust.method = "BH") %>%
  tibble::rownames_to_column("metabolite") %>%
  arrange(adj.P.Val)

cat("Significant metabolites (BH < 0.05):", sum(res_metab$adj.P.Val < 0.05), "\n")
cat("Nominal p<0.05 (uncorrected):", sum(res_metab$P.Value < 0.05), "of", nrow(res_metab), "\n")

write.csv(res_metab, file.path(RES_DIR, "limma_metabolomics_AsymAD_vs_AD_johnson_pathology.csv"), row.names = FALSE)

writeLines(asymad_ids_m, file.path(RES_DIR, "AsymAD_projid_metabolomics_johnson.txt"))
writeLines(ad_ids_m,     file.path(RES_DIR, "AD_projid_metabolomics_johnson.txt"))

# ── 7. volcano plot ──────────────────────────────────────────────────────
# drop unnamed "X - #####" features before anything else touches this —
# (convention used throughout this repository)
res_metab_named <- res_metab %>% filter(!grepl("^X\\s*-\\s*\\d+", metabolite))
res_metab_named$protein <- res_metab_named$metabolite
sig_high_vol <- res_metab_named %>% filter(adj.P.Val < 0.05)

label_high <- unique(c(
  sig_high_vol %>% filter(logFC > 0) %>% slice_min(adj.P.Val, n = 10) %>% pull(protein),
  sig_high_vol %>% filter(logFC < 0) %>% slice_min(adj.P.Val, n = 10) %>% pull(protein),
  sig_high_vol %>% filter(logFC > 0) %>% slice_max(logFC,   n = 10) %>% pull(protein),
  sig_high_vol %>% filter(logFC < 0) %>% slice_min(logFC,   n = 10) %>% pull(protein)
))

volcano_high <- res_metab_named %>%
  mutate(
    neg_log_p = -log10(adj.P.Val),
    direction = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Higher in AsymAD",
      adj.P.Val < 0.05 & logFC < 0 ~ "Higher in AD",
      TRUE ~ "ns"
    ),
    label = ifelse(protein %in% label_high, protein, NA)
  )

p_high <- ggplot(volcano_high, aes(x = logFC, y = neg_log_p, color = direction)) +
  geom_point(aes(size = ifelse(direction == "ns", 0.8, 1.2)), alpha = 0.6, show.legend = FALSE) +
  geom_point(data = filter(volcano_high, direction != "ns"), aes(size = 1.5), alpha = 0.9) +
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
    title = "AsymAD vs AD \u2014 Metabolomics (Johnson rubric + pathology-attached metadata)",
    x = "logFC (positive = higher in AsymAD)",
    y = expression(-log[10](adj.P.Val))
  ) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(face = "bold"))

p_high
ggsave(file.path(FIG_DIR, "volcano_AsymAD_vs_AD_metabolomics_johnson_pathology.png"), plot = p_high,
       dpi = 600, bg = "white", width = 16, height = 9)

