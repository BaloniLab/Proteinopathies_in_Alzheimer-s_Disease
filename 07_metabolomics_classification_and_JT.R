# ══════════════════════════════════════════════════════════════════════════
# Metabolomics donors: Johnson-rubric classification and Jonckheere-Terpstra trend test

#
#   (1) Classify the metabolomics-matched ROSMAP donors as Control / AsymAD / AD with the Johnson
#   rubric (CERAD, Braak and dementia = MMSE < 24) and split AsymAD into Resilient (cogdx 1) and
#   Vulnerable (cogdx 2-6). (2) Remove unnamed 'X - #####' features. (3) JT test, per metabolite, of a
#   monotonic order Control < Resilient < Vulnerable < AD on covariate-adjusted residuals
#   (age, sex, APOE genotype, PMI); BH across metabolites. (4) Plot of the top trending metabolites.
#   Creates the session objects `clin`, `pheno`, `expr_sub` and removes the X features from
#   `Expr_metabolomics`, which scripts 08 and 09 rely on.
#
# INPUTS
#   R objects: Expr_metabolomics (donors x metabolites; first column = projid), ROSMAP_clinical,
#   dataset_1744_cross_sectional_08_30_2026 (see data/README.md)
# OUTPUTS
#   results/metabolomics_donors_johnson_classification_cross.csv
#   results/jonckheere_trend_control_res_vol_ad_metab_johnson.csv
#   results/jt_trend_increasing_/decreasing_metab_johnson_BH05.txt
#   figures/jt_trend_top_metabolites_control_res_vol_ad_johnson.png
# RUN AFTER
#   none
#
# Run from the project root (see README.md). The analysis scripts keep their
# intermediate objects in the R session; run_all.R sources them in order.
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")
require_packages(c("clinfun", "dplyr", "tibble", "tidyr", "ggplot2", "ggpubr"))
load_inputs(c("Expr_metabolomics", "ROSMAP_clinical", "dataset_1744_cross_sectional_08_30_2026"))

library(clinfun)
library(dplyr)
library(tibble)
library(tidyr)       # pivot_longer() is used below; the script never loaded tidyr
library(ggplot2)
library(ggpubr)

if ("...1" %in% names(Expr_metabolomics)) {
  Expr_metabolomics <- Expr_metabolomics %>% column_to_rownames(var = "...1")
}


# ══════════════════════════════════════════════════════════════════════════
# Re-classify the 500 metabolomics-matched donors into Control / AsymAD / AD
# using dataset_1744_cross_sectional_08_30_2026 for CERAD/Braak/covariates,
# with MMSE (cts_mmse30_lv) JOINED IN from ROSMAP_clinical — dataset_1744 has
# no MMSE/CDR column at all, so this keeps the rubric's actual MMSE<24
# dementia criterion rather than substituting cogdx.
#
#   Rubric:
#     Control: CERAD 0-1, Braak 0-3, WITHOUT dementia (if Braak==3, CERAD
#              must be exactly 0)
#     AsymAD:  CERAD 1-3, Braak 3-6, WITHOUT dementia
#     AD:      CERAD 2-3, Braak 3-6, WITH dementia
#     Dementia = MMSE < 24 (CDR still unavailable in either source)
#
#   CERAD recoding: raw `ceradsc` is reverse-coded (1=Definite AD,
#   4=No AD) -> CERAD_std = 4 - ceradsc (ascending severity, 0-3).
# ══════════════════════════════════════════════════════════════════════════

library(dplyr)

# ── 0. the 500 metabolomics donors ──────────────────────────────────────────
metab_projids <- rownames(Expr_metabolomics)
cat("Metabolomics donors (projid):", length(metab_projids), "\n")

# ── 1. match against dataset_1744 ───────────────────────────────────────────
clin <- ROSMAP_clinical %>%
  mutate(projid = as.character(projid)) %>%
  filter(projid %in% as.character(metab_projids))

cat("Matched to ROSMAP_clinical:", nrow(clin),
    "/", length(metab_projids), "\n")
if (nrow(clin) < length(metab_projids)) {
  unmatched <- setdiff(as.character(metab_projids), clin$projid)
  cat("NOTE:", length(unmatched), "metabolomics projid(s) not found in dataset_1744.\n")
}


clin <- clin %>%
  mutate(CERAD_std = 4 - ceradsc,
         dementia  = cts_mmse30_lv < 24)


# ── 4. classify per the rubric ──────────────────────────────────────────────
clin <- clin %>%
  mutate(
    group = case_when(
      CERAD_std %in% c(0, 1) & braaksc %in% 0:3 & !dementia &
        (braaksc != 3 | CERAD_std == 0)                    ~ "Control",
      CERAD_std %in% c(1, 2, 3) & braaksc %in% 3:6 & !dementia ~ "AsymAD",
      CERAD_std %in% c(2, 3)    & braaksc %in% 3:6 &  dementia ~ "AD",
      TRUE ~ NA_character_
    )
  )
# ══════════════════════════════════════════════════════════════════════════
# Merge age_death from dataset_1744_cross_sectional_08_30_2026 into clin,
# by projid
# ══════════════════════════════════════════════════════════════════════════

library(dplyr)

clin <- clin %>%
  select(-any_of("age_death")) %>%   # drop if clin already has one, to avoid .x/.y suffixes
  left_join(
    dataset_1744_cross_sectional_08_30_2026 %>%
      mutate(projid = as.character(projid)) %>%
      select(projid, age_death),
    by = "projid"
  )

n_before <- nrow(clin)

clin <- clin %>%
  filter(!is.na(age_death), !is.na(group))

cat("clin:", nrow(clin), "/", n_before, "rows kept (dropped NA age_death or NA group)\n")
print(table(clin$group, useNA = "ifany"))
cat("age_death matched:", sum(!is.na(clin$age_death)), "/", nrow(clin), "\n")
cat("\n=== Group counts among the", nrow(clin), "metabolomics-matched donors ===\n")
print(table(clin$group, useNA = "ifany"))


write.csv(clin %>% select(projid, CERAD_std, braaksc, cts_mmse30_lv, dementia, group),
          file.path(RES_DIR, "metabolomics_donors_johnson_classification_cross.csv"), row.names = FALSE)
# (Removed: `clin <- metabolomics_donors_johnson_classification_cross`. That line replaced `clin` by an object
#  read back from the CSV just written, which only existed in the author's session and has no cogdx
#  column. `clin` already holds the same donors, with all columns, so it is used directly.)
# ── 0. build the 4-group, ordinal-scored pheno, Johnson-rubric cohort ──────
covar_source <- dataset_1744_cross_sectional_08_30_2026 %>%
  mutate(projid = as.character(projid)) %>%
  select(projid, cogdx, age_death, msex, apoe_genotype, pmi, amylsqrt_est_8reg, tangsqrt_est_8reg, tdp_st4)

control_group <- clin %>% filter(group == "Control") %>%
  mutate(projid = as.character(projid)) %>%
  transmute(projid, group = "Control")

ad_group <- clin %>% filter(group == "AD") %>%
  mutate(projid = as.character(projid)) %>%
  transmute(projid, group = "AD")

# AsymAD split into Resilient (cogdx==1, NCI) / Vulnerable (cogdx 2-6), same
# rubric as the protein version, cogdx sourced from dataset_1744 instead of
# `meta` since that's the metadata this cohort/classification is keyed to
# ══════════════════════════════════════════════════════════════════════════
# clin already has its own cogdx column (from dataset_1744) — no need to
# join it in again from covar_source, which was causing a cogdx.x/cogdx.y
# collision and leaving no plain "cogdx" column behind
# ══════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════
# Fill the one missing cogdx value in clin with 1, then rebuild
# asym_groups_metab (Resilient / Vulnerable split)
# ══════════════════════════════════════════════════════════════════════════

library(dplyr)

stopifnot("cogdx must be a column of ROSMAP_clinical" = "cogdx" %in% names(clin))
cat("NA cogdx before fill:", sum(is.na(clin$cogdx)), "\n")

# One donor has no cogdx; it was set to 1 (NCI). The guard below stops the script if the data ever
# contain more than that single donor, because the fill was only validated for one.
stopifnot("more than one donor has a missing cogdx" = sum(is.na(clin$cogdx)) <= 1)
if (any(is.na(clin$cogdx))) cat("Filling missing cogdx with 1 (NCI) for projid:", clin$projid[is.na(clin$cogdx)], "\n")

clin <- clin %>%
  mutate(cogdx = ifelse(is.na(cogdx), 1, cogdx))

cat("NA cogdx after fill: ", sum(is.na(clin$cogdx)), "\n")

asym_groups_metab <- clin %>%
  filter(group == "AsymAD", !is.na(cogdx)) %>%
  mutate(asymad_group = case_when(
    cogdx == 1     ~ "Resilient (NCI)",
    cogdx %in% 2:6 ~ "Vulnerable (MCI/Dementia)",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(asymad_group)) %>%
  transmute(projid, group = asymad_group)

cat("n Resilient:", sum(asym_groups_metab$group == "Resilient (NCI)"),
    " | n Vulnerable:", sum(asym_groups_metab$group == "Vulnerable (MCI/Dementia)"), "\n")



GROUP_LEVELS <- c("Control", "Resilient (NCI)", "Vulnerable (MCI/Dementia)", "AD")

pheno <- bind_rows(control_group, asym_groups_metab, ad_group) %>%
  mutate(group = factor(group, levels = GROUP_LEVELS),
         severity_score = as.numeric(group)) %>%
  left_join(covar_source, by = "projid") %>%
  transmute(projid, severity_score,
            age_death_z  = as.numeric(scale(age_death)),
            pmi_z        = as.numeric(scale(pmi)),
            msex         = factor(msex),
            apoe_genotype = factor(apoe_genotype)) %>%
  filter(!is.na(severity_score), projid %in% rownames(Expr_metabolomics))

cat("n per severity score (1=Control, 2=Resilient, 3=Vulnerable, 4=AD):\n")
print(table(pheno$severity_score))
# ══════════════════════════════════════════════════════════════════════════
# Remove unnamed "X - #####" metabolite features from Expr_metabolomics —
# doing it at the source matrix (not just a results table) so every
# downstream analysis (JT test, interaction test, DEG) is clean from here on
# ══════════════════════════════════════════════════════════════════════════

x_cols <- grep("^X\\s*-\\s*\\d+", colnames(Expr_metabolomics), value = TRUE)
cat("Removing", length(x_cols), "unnamed 'X - #####' metabolite column(s)\n")

Expr_metabolomics <- Expr_metabolomics[, !colnames(Expr_metabolomics) %in% x_cols, drop = FALSE]

cat("Metabolites remaining:", ncol(Expr_metabolomics), "\n")
expr_sub <- t(Expr_metabolomics[pheno$projid, , drop = FALSE])   # metabolites x samples
stopifnot(identical(colnames(expr_sub), pheno$projid))
cat("Total metabolites tested:", nrow(expr_sub), "\n")

# ── 1. one JT test per metabolite, on covariate-adjusted residuals ─────────
run_jt_one_metab <- function(metab_id) {
  d <- pheno
  d$expr <- as.numeric(expr_sub[metab_id, ])
  d <- d %>% filter(!is.na(expr))
  if (nrow(d) < 20 || n_distinct(d$severity_score) < 4) return(NULL)
  
  fit_cov <- try(lm(expr ~ age_death_z + msex + apoe_genotype + pmi_z, data = d), silent = TRUE)
  if (inherits(fit_cov, "try-error")) return(NULL)
  d$resid <- residuals(fit_cov)
  
  jt <- try(jonckheere.test(d$resid, d$severity_score, alternative = "two.sided"),
            silent = TRUE)
  if (inherits(jt, "try-error")) return(NULL)
  
  rho <- suppressWarnings(cor(d$resid, d$severity_score, method = "spearman"))
  
  data.frame(metabolite = metab_id, JT_stat = jt$statistic, p.value = jt$p.value,
             direction = ifelse(rho > 0, "Increasing", "Decreasing"),
             n = nrow(d))
}

cat("\nRunning JT test per metabolite — this loops one test at a time and may take a bit...\n")
res_jt_metab <- bind_rows(lapply(rownames(expr_sub), run_jt_one_metab)) %>%
  mutate(adj.P.Val = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value)

cat("\nSignificant trend (BH adj.P.Val < 0.05):", sum(res_jt_metab$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("Nominal p<0.05 (uncorrected):", sum(res_jt_metab$p.value < 0.05, na.rm = TRUE), "of", nrow(res_jt_metab), "\n")

write.csv(res_jt_metab, file.path(RES_DIR, "jonckheere_trend_control_res_vol_ad_metab_johnson.csv"), row.names = FALSE)

increasing_metabs_jt <- res_jt_metab %>% filter(adj.P.Val < 0.05, direction == "Increasing") %>% pull(metabolite)
decreasing_metabs_jt <- res_jt_metab %>% filter(adj.P.Val < 0.05, direction == "Decreasing") %>% pull(metabolite)

cat("Continuously INCREASING (JT, BH<0.05):", length(increasing_metabs_jt), "\n")
cat("Continuously DECREASING (JT, BH<0.05):", length(decreasing_metabs_jt), "\n")

writeLines(increasing_metabs_jt, file.path(RES_DIR, "jt_trend_increasing_metab_johnson_BH05.txt"))
writeLines(decreasing_metabs_jt, file.path(RES_DIR, "jt_trend_decreasing_metab_johnson_BH05.txt"))


# ── 3. visualize top hits, with pairwise brackets ────────────────────────
TOP_N <- 8
top_trend_jt <- res_jt_metab %>% slice_min(adj.P.Val, n = TOP_N) %>% pull(metabolite)

plot_data <- as.data.frame(t(expr_sub[top_trend_jt, , drop = FALSE])) %>%
  rownames_to_column("projid") %>%
  left_join(pheno %>% mutate(group = factor(severity_score, labels = GROUP_LEVELS)) %>%
              select(projid, group),
            by = "projid") %>%
  pivot_longer(-c(projid, group), names_to = "metabolite", values_to = "expression") %>%
  filter(!is.na(expression)) %>%
  mutate(metabolite = factor(metabolite, levels = top_trend_jt))

pal4 <- c("Control" = "skyblue3", "Resilient (NCI)" = "#0072B2",
          "Vulnerable (MCI/Dementia)" = "#CC79A7", "AD" = "red4")

comparison_list <- list(
  c("Control", "Resilient (NCI)"),
  c("Resilient (NCI)", "Vulnerable (MCI/Dementia)"),
  c("Vulnerable (MCI/Dementia)", "AD")
)

ann_jt <- res_jt_metab %>%
  filter(metabolite %in% top_trend_jt) %>%
  mutate(metabolite = factor(metabolite, levels = top_trend_jt),
         label      = sprintf("adj.p = %.3g", adj.P.Val))

# Brackets: ggpubr::stat_compare_means default (two-sided Wilcoxon rank-sum, UNADJUSTED), descriptive only;
# the JT adj.p printed in each facet is the test of record for the trend.
p_jt <- ggplot(plot_data, aes(group, expression, color = group, fill = group)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, alpha = 0.4, size = 1.3) +
  stat_compare_means(comparisons = comparison_list, label = "p.signif",
                     tip.length = 0.01, size = 6, step.increase = 0.10,
                     bracket.size = 0.4) +
  geom_text(data = ann_jt, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.4,
            size = 3, fontface = "bold", color = "black") +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  scale_color_manual(values = pal4, name = NULL) +
  scale_fill_manual(values  = pal4, name = NULL) +
  labs(x = NULL, y = NULL,
       title = paste0("Top ", TOP_N, " metabolites by Jonckheere-Terpstra trend, Control \u2192 AD (Johnson rubric)")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        strip.text  = element_text(face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold", color = "black"),
        axis.title  = element_text(face = "bold", color = "black"),
        plot.title  = element_text(face = "bold", hjust = 0.5),
        legend.text = element_text(face = "bold"))

p_jt
dir.create("figures", showWarnings = FALSE)
ggsave("figures/jt_trend_top_metabolites_control_res_vol_ad_johnson.png", p_jt,
       width = 18, height = 12, dpi = 600, bg = "white")

