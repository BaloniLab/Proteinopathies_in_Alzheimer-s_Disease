# «Manuscript title» -- analysis code

> **Before archiving, fill in every «placeholder»** (title, authors, ORCIDs, affiliations, manuscript DOI, license
> holder) in this file, `CITATION.cff`, `.zenodo.json` and `LICENSE`, and read **"Items to verify"** below.

R code for the analyses of ROSMAP and SEA-AD data on cognitive resilience in asymptomatic Alzheimer's disease
(AsymAD): cognitive trajectories, neuropathology, proteomic and metabolomic differential abundance,
Jonckheere-Terpstra trend tests and protein-metabolite correlation networks.

**No participant-level data are included.** ROSMAP data are controlled-access; see `data/README.md`.

## Layout

```
R/00_setup.R        shared setup (folders, seed, package check, input loader)
run_all.R           runs scripts 01-10 in dependency order, in one R session
scripts/            the analysis scripts (below)
data/               inputs go here (not distributed)   data/reactome/  Enrichr exports used by scripts 06 and 11
figures/            PNG figures written by the scripts
results/            result tables, gene lists, networks written by the scripts
source_data/        tables behind the figures
CHANGES_FROM_ORIGINAL.md   every change made while preparing this release
```

## Scripts

| Script | What it does | Main outputs |
|---|---|---|
| `01_fig1_cognitive_trajectories.R` | Control / AsymAD / AD cognitive trajectories (years before death, follow-up year), mixed models (diagnosis x time + age at death, random intercept and slope), per-group slopes (emmeans) | `figures/fig_global_*`, `fig_domains_*`; `source_data/Table1_*`, `Fig1-4_*`, `FigS_*` |
| `02_fig1_pathology_cognition_heatmap.R` | Mean pathology and cognition per diagnosis; Dunn's test (BH) vs Control | `figures/pathology_heatmap_8reg*.png`; `source_data/Fig1_heatmap_*` |
| `03_fig1_trajectory_alluvial.R` | Clinical state (NCI / MCI / dementia) at -15 ... 0 years before death flowing into the final classification | `figures/trajectory_to_endpoint.png` |
| `04_fig2_pathology_burden_boxplots.R` | Amyloid / tangle / TDP-43 burden by diagnosis and by stage; Dunn's test (BH) | `figures/fig_*_by_diagnosis.png`, `fig_*_by_stage.png`; `source_data/pathology_*` |
| `05_proteomics_DEG_full_cohort.R` | limma, AsymAD vs AD, full cohort; volcano | `results/limma_AsymAD_vs_AD_full_cohort.csv` |
| `06_proteomics_DEG_stratified_by_pathology.R` | limma within CERAD / Braak / TDP-43 strata; Venn diagrams; robust and unique proteins; Reactome dot plots | `results/limma_AsymAD_vs_AD_<stratum>.csv`, `results/unique_protein/*`, `figures/*` |
| `07_metabolomics_classification_and_JT.R` | Johnson-rubric classification of the metabolomics donors; Resilient / Vulnerable split; Jonckheere-Terpstra trend test per metabolite | `results/jonckheere_trend_..._metab_johnson.csv` |
| `08_metabolomics_DEG_full_cohort.R` | limma, AsymAD vs AD (Johnson rubric); volcano | `results/limma_metabolomics_AsymAD_vs_AD_johnson_pathology.csv` |
| `09_metabolomics_DEG_stratified_by_pathology.R` | limma within pathology strata; Venn / UpSet; multi-stratum heatmaps | `results/limma_metab_AsymAD_vs_AD_<stratum>.csv`, `figures/*` |
| `10_protein_metabolite_correlation_network.R` | Spearman protein x metabolite networks (main and relaxed) on dual-omics donors; Cytoscape / GraphML files | `results/cytoscape_*`, `results/network_*.graphml` |
| `11_proteomics_JT_pathway_heatmaps.R` | One heatmap per Reactome pathway of the JT-increasing / decreasing proteins; all six pairwise group comparisons | `figures/pathway_heatmaps/`, `results/pathway_heatmap_source_data.csv` |

Manuscript figure / table numbers for each script: «to be completed by the authors».

## Requirements and installation

Developed with **R 4.4** and **ggplot2 3.5.1**; other versions have not been tested. The scripts do not install
packages. Install once:

```r
install.packages(c("tidyverse", "readxl", "conflicted", "ggalluvial", "rstatix", "ggpubr", "lmerTest", "emmeans",
                   "clinfun", "ggvenn", "ggrepel", "patchwork", "circlize", "igraph", "reshape2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "ComplexHeatmap"))
```

For a reproducible archive, run `renv::init()` / `renv::snapshot()` after a successful run and add `renv.lock` to
the release. Each script also records `sessionInfo()` where it matters (`source_data/sessionInfo*.txt`).

## Running

1. Obtain the input data and put it in `data/` (see `data/README.md`); edit `INPUT_FILES` in `R/00_setup.R` where marked.
2. From the project root: `source("run_all.R")` (or run `scripts/01_...` to `scripts/10_...` in order in one session).
   Scripts 08-10 use objects created by scripts 05 and 07, so they must share an R session; each stops with a clear
   message if its prerequisite is missing. Script 11 needs the proteomics JT objects (see "Not included").

## Analysis choices that differ between the proteomic and metabolomic analyses

* **Cohort definition.** Proteomics uses the Emory diagnosis labels (`EmoryStrictDx.2019`). Metabolomics uses the
  Johnson rubric (Control: CERAD 0-1, Braak 0-3, no dementia, and CERAD 0 if Braak = 3; AsymAD: CERAD 1-3, Braak 3-6,
  no dementia; AD: CERAD 2-3, Braak 3-6, dementia = MMSE < 24; `CERAD_std = 4 - ceradsc`). AsymAD is split into
  Resilient (`cogdx` = 1) and Vulnerable (`cogdx` 2-6). One donor without `cogdx` is set to 1 (script 07 reports and guards this).
* **Covariates.** Proteomics: age at death (z), sex, APOE risk category, PMI (z). Metabolomics: age at death (z), sex,
  APOE genotype, PMI (z). Amyloid / tangle burden is deliberately not a covariate (it partly defines the groups).
* **Presence filters.** Proteomics scripts 05 and 06 use `MIN_PCT_PRESENT <- 0`, i.e. they keep every protein of
  `expr_filtered` (filtering happened upstream). Metabolomics keeps features present in >= 80 % of AsymAD or of AD donors
  (computed in the full cohort and again within each stratum). Unnamed `X - #####` metabolites are removed.
* **Multiple testing.** Benjamini-Hochberg within each limma model across all features. Venn / UpSet / overlap sets
  and the multi-stratum heatmaps use raw p < 0.05 (exploratory, and labelled so in the figures).
* **Jonckheere-Terpstra.** `clinfun::jonckheere.test` (two-sided) on residuals of
  `expr ~ age + sex + APOE genotype + PMI`; ordering Control < Resilient < Vulnerable < AD; direction from the sign of
  the Spearman correlation with the ordering; BH across features.
* **Stage panels (script 04).** Dunn's test within each stage uses every diagnosis group with >= 3 donors in that
  stage (Control included when present) and is BH-adjusted across the pairs tested in that stage; a stage is skipped
  if any group present has fewer than 3 donors. Diagnosis-panel brackets use the BH-adjusted Dunn table.
* **Randomness.** `set.seed(2026)` affects only the jitter of points in plots.

## Items to verify before you archive

1. **`TDP-43 (8reg)` label (script 02):** the column is `tdp_cs_6reg` (6-region score) but is labelled 8-region.
2. **Presence filter in the proteomics scripts:** the original comments and messages said "80 %" while the code used 0.
   The code is unchanged (0); the wording is corrected. State in the Methods how `expr_filtered` was filtered.
3. **Unnamed metabolites in scripts 08 / 09:** these scripts now also remove the `X - #####` features (a no-op after
   script 07). If your published numbers for these models were produced with those features still present, delete the
   added lines (marked in `CHANGES_FROM_ORIGINAL.md`).
4. **Robust proteins (script 06):** `robust_up` / `robust_down` were not defined in the original script; they are now
   computed from the triple-pathology sets. Confirm that the printed counts are 19 and 3.
5. **Reactome table numbering (`data/reactome/README.md`):** the table-7 ... table-18 mapping was reconstructed from upload order.
6. **Statistics were not re-run for this release.** Re-run the pipeline and compare against the published numbers.
7. **Data sharing.** Do **not** deposit `data/` or the donor-level files listed in `.gitignore`
   (`results/metadata_metabolomics_johnson_with_pathology.csv`, `results/metabolomics_donors_johnson_classification_cross.csv`,
   `results/*projid*.txt`); check the ROSMAP data use agreement before sharing any `source_data/` or `results/` table.

## Not included in this release

* The proteomics Jonckheere-Terpstra script (`jonckheere_test_control_res_vol_ad.R`), which creates `expr_sub`, `pheno`
  and `res_jt` for script 11 and for the JT results used elsewhere.
* Scripts for other figures/tables that were not part of the files provided (e.g. the one-row ROSMAP pathology
  figure, the SEA-AD regional plots, the module heatmaps and pairwise-comparison summaries).

## License and citation

Code: «MIT License (see `LICENSE`) -- confirm with all authors and your institution». Please cite the manuscript
and this archive (`CITATION.cff`; Zenodo DOI «10.5281/zenodo.XXXXXXX»).
