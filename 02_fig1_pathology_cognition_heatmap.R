# ══════════════════════════════════════════════════════════════════════════
# Pathology burden and cognition heatmap by diagnosis (Figure 1)
#
#   Mean pathology and cognition per diagnosis (complete cases), z-scored per row for colour,
#   Dunn's test (BH-adjusted) of AsymAD and AD against Control for the stars; two layouts.
#
# INPUTS
#   R object: final_metadata_complete (see data/README.md)
# OUTPUTS
#   figures/pathology_heatmap_8reg.png, figures/pathology_heatmap_8reg_flipped.png
#   source_data/Fig1_heatmap_8reg_values.csv, source_data/Fig1_heatmap_8reg_dunn.csv
# RUN AFTER
#   none
#
# Run from the project root (see README.md). The analysis scripts keep their
# intermediate objects in the R session; run_all.R sources them in order.
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")
require_packages(c("ggplot2", "dplyr", "tidyr", "rstatix", "conflicted", "scales"))
load_inputs(c("final_metadata_complete"))

library(ggplot2); library(dplyr); library(tidyr); library(rstatix)
library(conflicted)   # needed for conflicts_prefer()
conflicts_prefer(base::unname, .quiet = TRUE)



meta_h <- final_metadata_complete

# join sanity — must be non-zero, and rows must not have grown
cat("Non-missing amyloid / tangle estimates:", sum(!is.na(meta_h$amylsqrt_est_8reg)), "/", sum(!is.na(meta_h$tangsqrt_est_8reg)), "of", nrow(meta_h), "donors\n")

# ── Column mapping (single source of truth) ──────────────────────────────────
COL_MAP <- c(
  "CERAD_std"         = "Amyloid β (CERAD)",
  "amylsqrt_est_8reg" = "Amyloid (8reg)",
  "braaksc"           = "Tau (Braak Stage)",
  "tangsqrt_est_8reg" = "Tangles (8reg)",
  "tdp_st4"           = "TDP-43 Stage",
  "tdp_cs_6reg"       = "TDP-43 (8reg)",   # NOTE: the column is tdp_cs_6reg (6-region score); the label is kept as published -- confirm before release
  "gpath"             = "Global Path",
  "cogn_global_lv"    = "Global Cognition"
)

RAW_COLS   <- names(COL_MAP)
CLEAN_LABS <- as.character(COL_MAP)
DX_LEVELS  <- c("Control", "AsymAD", "AD")

# ── how many donors survive the complete-case filter? ────────────────────────
print(meta_h %>%
  dplyr::select(EmoryStrictDx.2019, all_of(RAW_COLS)) %>%
  summarise(across(all_of(RAW_COLS), ~ sum(!is.na(.x)))) %>%
  t())

complete_n <- meta_h %>%
  filter(if_all(all_of(RAW_COLS), ~ !is.na(.))) %>%
  count(EmoryStrictDx.2019)
print(complete_n)

# ── 1. Long data ─────────────────────────────────────────────────────────────
sig_data <- meta_h %>%
  dplyr::select(EmoryStrictDx.2019, all_of(RAW_COLS)) %>%
  filter(if_all(all_of(RAW_COLS), ~ !is.na(.))) %>%
  mutate(EmoryStrictDx.2019 = factor(EmoryStrictDx.2019, levels = DX_LEVELS)) %>%
  pivot_longer(cols = all_of(RAW_COLS),
               names_to = "Pathology", values_to = "Score") %>%
  mutate(Pathology = recode(Pathology, !!!COL_MAP),
         Pathology = factor(Pathology, levels = CLEAN_LABS))

# ── 2. Dunn's test (vs Control) ──────────────────────────────────────────────
p_to_stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

dunn_res <- sig_data %>%
  group_by(Pathology) %>%
  dunn_test(Score ~ EmoryStrictDx.2019, p.adjust.method = "BH") %>%
  ungroup()

heatmap_stars <- dunn_res %>%
  filter(group1 == "Control") %>%
  mutate(stars = p_to_stars(p.adj)) %>%
  dplyr::select(Pathology, group2, stars) %>%
  rename(EmoryStrictDx.2019 = group2) %>%
  mutate(EmoryStrictDx.2019 = factor(EmoryStrictDx.2019, levels = DX_LEVELS))

# ── 3. Summary + Z-score + join stars ────────────────────────────────────────
heatmap_data <- sig_data %>%
  group_by(EmoryStrictDx.2019, Pathology) %>%
  summarise(Mean = mean(Score, na.rm = TRUE),
            SE   = sd(Score, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  group_by(Pathology) %>%
  mutate(Z_score = scale(Mean)[, 1]) %>%
  ungroup() %>%
  left_join(heatmap_stars, by = c("EmoryStrictDx.2019", "Pathology")) %>%
  mutate(stars = ifelse(is.na(stars), "", stars),
         label = sprintf("%.2f\n(\u00b1%.2f)\n%s", Mean, SE, stars))

# ── 4. Plot ──────────────────────────────────────────────────────────────────
p <- ggplot(heatmap_data,
            aes(x = EmoryStrictDx.2019, y = Pathology, fill = Z_score)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = label), size = 6, color = "black", fontface = "bold") +
  geom_hline(yintercept = 1.5, color = "grey40", linewidth = 0.8, linetype = "dashed") +
  scale_fill_gradientn(colors = c("#2166AC", "#F7F7F7", "#B2182B"),
                       values = scales::rescale(c(-2, 0, 2)),
                       name = "Z-score\n(relative burden)") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  guides(fill = guide_colorbar(barwidth = unit(8, "cm"), barheight = unit(0.8, "cm"),
                               title.position = "left", title.vjust = 0.8)) +
  labs(subtitle = "", x = NULL, y = NULL) +
  theme_minimal(base_size = 22) +
  theme(plot.subtitle   = element_text(hjust = 0.5, color = "grey40", size = 14),
        axis.text       = element_text(face = "bold", color = "black"),
        axis.title      = element_text(face = "bold", color = "black"),
        legend.position = "bottom",
        legend.text     = element_text(face = "bold", color = "black", size = 18),
        legend.title    = element_text(face = "bold", color = "black", size = 18),
        panel.grid      = element_blank())

p
ggsave(file.path(FIG_DIR, "pathology_heatmap_8reg.png"), p, dpi = 600, bg = "white", width = 10, height = 11)

# ── 5. source data ───────────────────────────────────────────────────────────
write.csv(heatmap_data, file.path(SRC_DIR, "Fig1_heatmap_8reg_values.csv"), row.names = FALSE)
write.csv(dunn_res,     file.path(SRC_DIR, "Fig1_heatmap_8reg_dunn.csv"),   row.names = FALSE)




p <- ggplot(heatmap_data,
            aes(x = Pathology, y = EmoryStrictDx.2019, fill = Z_score)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = label), size = 5.5, color = "black", fontface = "bold") +
  geom_vline(xintercept = length(CLEAN_LABS) - 0.5,
             color = "grey40", linewidth = 0.8, linetype = "dashed") +
  scale_fill_gradientn(colors = c("#2166AC", "#F7F7F7", "#B2182B"),
                       values = scales::rescale(c(-2, 0, 2)),
                       name = "Z-score\n(relative burden)") +
  scale_x_discrete(limits = CLEAN_LABS, position = "bottom") +
  scale_y_discrete(limits = rev(DX_LEVELS)) +
  guides(fill = guide_colorbar(barwidth = unit(8, "cm"), barheight = unit(0.8, "cm"),
                               title.position = "left", title.vjust = 0.8)) +
  labs(subtitle = "", x = NULL, y = NULL) +
  theme_minimal(base_size = 22) +
  theme(plot.subtitle   = element_text(hjust = 0.5, color = "grey40", size = 14),
        axis.text.x.bottom = element_text(face = "bold", color = "black",
                                       angle = 45, hjust = 1, vjust = 1),
        axis.text.y     = element_text(face = "bold", color = "black"),
        legend.position = "bottom",
        legend.text     = element_text(face = "bold", color = "black", size = 18),
        legend.title    = element_text(face = "bold", color = "black", size = 18),
        panel.grid      = element_blank())

p
ggsave(file.path(FIG_DIR, "pathology_heatmap_8reg_flipped.png"), p, width = 16, height = 7, dpi = 600, bg = "white")
