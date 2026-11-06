# Input data

**No participant-level data are included in this repository.** ROSMAP data are controlled-access
(Rush Alzheimer's Disease Center Resource Sharing Hub / AD Knowledge Portal; a data use agreement is required).
Place the files you obtain in this folder (file names are registered in `R/00_setup.R`, object `INPUT_FILES`;
edit the entries marked `EDIT`).

The scripts use the following objects. Columns are those referenced in the code (other columns may be present).

| Object | File (default name) | Columns used |
|---|---|---|
| `meta_cc_8reg_n323_2026_09_02` | `meta_cc_8reg_n323_2026-09-02.csv` | `...1` (first, unnamed column = proteomics sample id), `projid`, `EmoryStrictDx.2019` (Control / AsymAD / AD), `age_death`, `msex`, `apoe_risk`, `pmi`, `CERAD_std`, `braaksc`, `tdp_st4`, `amylsqrt_est_8reg`, `tangsqrt_est_8reg`, `tdp_cs_6reg` |
| `meta_filtered` | EDIT | `projid`, `EmoryStrictDx.2019`, `age_death` |
| `final_metadata_complete` | EDIT | `EmoryStrictDx.2019`, `CERAD_std`, `amylsqrt_est_8reg`, `braaksc`, `tangsqrt_est_8reg`, `tdp_st4`, `tdp_cs_6reg`, `gpath`, `cogn_global_lv` |
| `dataset_1744_long_08_30_2026` | `dataset_1744_long_08-30-2026.xlsx` | `projid`, `age_at_visit_int`, `fu_year`, `dcfdx`, `cogn_global`, `cogn_ep`, `cogn_po`, `cogn_ps`, `cogn_se`, `cogn_wo` |
| `dataset_1744_cross_sectional_08_30_2026` | EDIT (inferred `dataset_1744_cross_sectional_08-30-2026.xlsx`) | `projid`, `age_death`, `msex`, `educ`, `apoe_genotype`, `pmi`, `ceradsc`, `braaksc`, `tdp_st4`, `cogdx`, `amylsqrt_est_8reg`, `tangsqrt_est_8reg` |
| `ROSMAP_clinical` | EDIT | `projid`, `ceradsc`, `braaksc`, `cts_mmse30_lv`, `cogdx` |
| `expr_filtered` | EDIT | proteins x donors; first column (`...1`) = protein id; remaining columns = donor ids matching `meta_cc_8reg_n323_2026_09_02$...1`; abundance values (NA allowed) |
| `Expr_metabolomics` | EDIT | donors x metabolites; first column (`...1`) = `projid`; remaining columns = metabolite names (including unnamed `X - #####` features, removed by the scripts); values z-scored per metabolite |

Coding notes used by the scripts: `CERAD_std = 4 - ceradsc` (the raw `ceradsc` is reverse-coded: 1 = definite AD,
4 = no AD); `msex` 0 = female, 1 = male.

Reactome / Enrichr enrichment exports used by scripts 06 and 11 go in `data/reactome/` (see the README there).
