# ══════════════════════════════════════════════════════════════════════════
# run_all.R -- runs the analysis scripts in dependency order, in ONE R session
# (scripts 08-10 use objects created by scripts 05 and 07, so they must share a session).
# Usage, from the project root:   source("run_all.R")
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")

scripts <- c(
  "scripts/01_fig1_cognitive_trajectories.R",
  "scripts/02_fig1_pathology_cognition_heatmap.R",
  "scripts/03_fig1_trajectory_alluvial.R",
  "scripts/04_fig2_pathology_burden_boxplots.R",
  "scripts/05_proteomics_DEG_full_cohort.R",
  "scripts/06_proteomics_DEG_stratified_by_pathology.R",
  "scripts/07_metabolomics_classification_and_JT.R",
  "scripts/08_metabolomics_DEG_full_cohort.R",
  "scripts/09_metabolomics_DEG_stratified_by_pathology.R",
  "scripts/10_protein_metabolite_correlation_network.R"
)
for (s in scripts) {
  message("\n==================== ", s, " ====================")
  source(s, echo = FALSE)
}

# Script 11 needs expr_sub / pheno / res_jt from the PROTEOMICS Jonckheere-Terpstra script
# (jonckheere_test_control_res_vol_ad.R), which is not part of this release. Source that script first
# (it must leave `pheno` with columns sample_id and group_label, and `expr_sub`, `res_jt` in the session).
if (exists("res_jt") && exists("expr_sub") && "sample_id" %in% names(pheno)) {
  source("scripts/11_proteomics_JT_pathway_heatmaps.R", echo = FALSE)
} else {
  message("\nSkipping scripts/11_proteomics_JT_pathway_heatmaps.R: the proteomics JT objects ",
          "(res_jt, expr_sub, pheno with sample_id) are not in the session.")
}

writeLines(capture.output(sessionInfo()), file.path(SRC_DIR, "sessionInfo_run_all.txt"))
