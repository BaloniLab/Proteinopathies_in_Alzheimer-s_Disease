# ══════════════════════════════════════════════════════════════════════════
# Pathology burden by diagnosis and by stage (ROSMAP): boxplots and statistics

#
#   Amyloid (8-region), tangle (8-region) and TDP-43 (6-region) burden by diagnosis (all pairs) and
#   by CERAD / Braak / TDP-43 stage (AsymAD vs AD within stage); one PNG per panel.
#   Brackets in BOTH panel types come from Dunn's test with Benjamini-Hochberg adjustment (the
#   tables written to source_data/)
#
# INPUTS
#   R object: meta_cc_8reg_n323_2026_09_02 (see data/README.md)
# OUTPUTS
#   figures/fig_<amyloid|tangles|tdp43>_by_diagnosis.png, fig_<...>_by_stage.png
#   source_data/pathology_*_dunn.csv, pathology_*_summary.csv, pathology_stage_counts.csv
# RUN AFTER
#   none
#
# Run from the project root (see README.md). The analysis scripts keep their
# intermediate objects in the R session; run_all.R sources them in order.
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")
require_packages(c("dplyr", "tidyr", "ggplot2", "rstatix"))
load_inputs(c("meta_cc_8reg_n323_2026_09_02"))

library(dplyr); library(tidyr); library(ggplot2); library(rstatix)

dx_cols   <- c(Control = "#0072B2", AsymAD = "#6A1B9A", AD = "red4")
DX_LEVELS <- c("Control", "AsymAD", "AD")

dat_p <- meta_cc_8reg_n323_2026_09_02 %>%
  filter(!is.na(EmoryStrictDx.2019)) %>%
  mutate(dx = factor(EmoryStrictDx.2019, levels = DX_LEVELS))

specs <- list(
  list(b = "amylsqrt_est_8reg", s = "CERAD_std", y = "Amyloid burden",
       x = "CERAD stage",  tag = "amyloid"),
  list(b = "tangsqrt_est_8reg", s = "braaksc",   y = "Tangle burden",
       x = "Braak stage",  tag = "tangles"),
  list(b = "tdp_cs_6reg",       s = "tdp_st4",   y = "TDP-43 burden",
       x = "TDP-43 stage", tag = "tdp43")
)

base_theme <- theme_classic(base_size = 20) +
  theme(axis.text    = element_text(face = "bold", color = "black"),
        axis.title   = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

# ── 1. stats behind panel a (computed first, because plot_a reads dx_tests) ──
dx_tests <- lapply(specs, function(s) {
  dat_p %>%
    filter(!is.na(.data[[s$b]])) %>%
    dunn_test(as.formula(paste(s$b, "~ dx")), p.adjust.method = "BH") %>%
    mutate(measure = s$y, .before = 1)
}) %>% bind_rows()

dx_summary <- lapply(specs, function(s) {
  dat_p %>%
    filter(!is.na(.data[[s$b]])) %>%
    group_by(dx) %>%
    summarise(n = n(),
              median = median(.data[[s$b]]),
              q25    = quantile(.data[[s$b]], 0.25),
              q75    = quantile(.data[[s$b]], 0.75),
              .groups = "drop") %>%
    mutate(measure = s$y, .before = 1)
}) %>% bind_rows()

stage_counts <- lapply(specs, function(s) {
  dat_p %>%
    filter(!is.na(.data[[s$b]]), !is.na(.data[[s$s]])) %>%
    count(stage = as.character(.data[[s$s]]), dx) %>%
    mutate(measure = s$y, .before = 1)
}) %>% bind_rows()

print(as.data.frame(dx_tests %>% select(measure, group1, group2, n1, n2, p, p.adj, p.adj.signif)))

# ── 2. stats behind panel b: within-stage diagnosis comparisons ──────────────
safe_dunn <- function(d, formula_str, group_col) {
  d <- d %>% filter(!is.na(.data[[group_col]])) %>% droplevels()
  if (n_distinct(d[[group_col]]) < 2) return(NULL)
  if (any(table(d[[group_col]]) < 3))  return(NULL)
  out <- try(dunn_test(d, as.formula(formula_str), p.adjust.method = "BH"),
             silent = TRUE)
  if (inherits(out, "try-error") || nrow(out) == 0) return(NULL)
  out
}

within_stage <- lapply(specs, function(s) {
  d0 <- dat_p %>%
    filter(!is.na(.data[[s$b]]), !is.na(.data[[s$s]])) %>%
    mutate(stage_f = factor(.data[[s$s]]))

  split(d0, d0$stage_f) %>%
    lapply(function(dd) {
      res <- safe_dunn(dd, paste(s$b, "~ dx"), "dx")
      if (is.null(res)) return(NULL)
      res %>% mutate(measure = s$y,
                     stage   = as.character(dd$stage_f[1]), .before = 1)
    }) %>%
    bind_rows()
}) %>% bind_rows()

# stage-to-stage comparisons within each diagnosis
across_stage <- lapply(specs, function(s) {
  d0 <- dat_p %>%
    filter(!is.na(.data[[s$b]]), !is.na(.data[[s$s]])) %>%
    mutate(stage_f = factor(.data[[s$s]]))

  split(d0, d0$dx) %>%
    lapply(function(dd) {
      res <- safe_dunn(dd, paste(s$b, "~ stage_f"), "stage_f")
      if (is.null(res)) return(NULL)
      res %>% mutate(measure = s$y,
                     dx      = as.character(dd$dx[1]), .before = 1)
    }) %>%
    bind_rows()
}) %>% bind_rows()

# medians per cell, so the tables are readable on their own
stage_summary <- lapply(specs, function(s) {
  dat_p %>%
    filter(!is.na(.data[[s$b]]), !is.na(.data[[s$s]])) %>%
    group_by(stage = as.character(.data[[s$s]]), dx) %>%
    summarise(n = n(),
              median = round(median(.data[[s$b]]), 3),
              q25    = round(quantile(.data[[s$b]], 0.25), 3),
              q75    = round(quantile(.data[[s$b]], 0.75), 3),
              .groups = "drop") %>%
    mutate(measure = s$y, .before = 1)
}) %>% bind_rows()

write.csv(dx_tests,      file.path(SRC_DIR, "pathology_by_diagnosis_dunn.csv"),    row.names = FALSE)
write.csv(dx_summary,    file.path(SRC_DIR, "pathology_by_diagnosis_summary.csv"), row.names = FALSE)
write.csv(stage_counts,  file.path(SRC_DIR, "pathology_stage_counts.csv"),         row.names = FALSE)
write.csv(within_stage,  file.path(SRC_DIR, "pathology_within_stage_dunn.csv"),    row.names = FALSE)
write.csv(across_stage,  file.path(SRC_DIR, "pathology_across_stage_dunn.csv"),    row.names = FALSE)
write.csv(stage_summary, file.path(SRC_DIR, "pathology_by_stage_summary.csv"),     row.names = FALSE)

print(as.data.frame(within_stage %>%
        select(measure, stage, group1, group2, n1, n2, p, p.adj, p.adj.signif)))

# ── 3. panel a: burden ~ diagnosis, brackets from the BH-adjusted Dunn table ──
plot_a <- function(burden, ylab, stats_tbl = dx_tests, measure_lab = ylab) {
  d    <- dat_p %>% filter(!is.na(.data[[burden]]))
  ymax <- max(d[[burden]], na.rm = TRUE)

  # bracket order bottom -> top: Control-AsymAD, Control-AD, AsymAD-AD ("ns" brackets are drawn)
  ann <- stats_tbl %>%
    filter(measure == measure_lab) %>%
    mutate(xmin = pmin(match(group1, DX_LEVELS), match(group2, DX_LEVELS)),
           xmax = pmax(match(group1, DX_LEVELS), match(group2, DX_LEVELS)),
           ord  = match(paste(xmin, xmax), c("1 2", "1 3", "2 3")),
           ybar = ymax * (1.05 + 0.11 * (ord - 1))) %>%
    filter(!is.na(ord))

  tick <- ymax * 0.02
  top  <- max(ann$ybar) + ymax * 0.12

  ggplot(d, aes(dx, .data[[burden]], fill = dx, color = dx)) +
    geom_boxplot(width = 0.62, alpha = 0.35, outlier.shape = NA, linewidth = 0.8) +
    geom_jitter(width = 0.15, size = 2.2, alpha = 0.6) +
    geom_segment(data = ann, aes(x = xmin, xend = xmax, y = ybar, yend = ybar),
                 inherit.aes = FALSE, linewidth = 0.6, color = "black") +
    geom_segment(data = ann, aes(x = xmin, xend = xmin, y = ybar, yend = ybar - tick),
                 inherit.aes = FALSE, linewidth = 0.6, color = "black") +
    geom_segment(data = ann, aes(x = xmax, xend = xmax, y = ybar, yend = ybar - tick),
                 inherit.aes = FALSE, linewidth = 0.6, color = "black") +
    geom_text(data = ann, aes(x = (xmin + xmax) / 2, y = ybar, label = p.adj.signif),
              inherit.aes = FALSE, size = 9, fontface = "bold", color = "black",
              vjust = -0.1) +
    scale_fill_manual(values = dx_cols) +
    scale_color_manual(values = dx_cols) +
    coord_cartesian(ylim = c(NA, top)) +
    labs(x = "Diagnosis", y = ylab) +
    base_theme + theme(legend.position = "none")
}

# ── 4. panel b: burden ~ stage, AsymAD vs AD within stage ────────────────────
plot_b <- function(burden, stage, ylab, xlab, min_n = 3, offset = 0.26,
                   stats_tbl = within_stage, measure_lab = NULL) {
  d <- dat_p %>%
    filter(!is.na(.data[[burden]]), !is.na(.data[[stage]])) %>%
    mutate(stage_f = factor(.data[[stage]])) %>%
    group_by(stage_f, dx) %>% filter(n() >= min_n) %>% ungroup() %>%
    mutate(xpos = as.integer(stage_f) + (as.integer(dx) - 2) * offset)

  ymax  <- max(d[[burden]], na.rm = TRUE)
  dxpos <- function(g) (match(g, DX_LEVELS) - 2) * offset

  # top of the data within each stage x pair, so brackets sit just above the boxes
  ann <- stats_tbl %>%
    filter(measure == measure_lab, p.adj.signif != "ns") %>%
    mutate(stage_f = factor(stage, levels = levels(d$stage_f))) %>%
    filter(!is.na(stage_f)) %>%
    rowwise() %>%
    mutate(xmin = as.integer(stage_f) + dxpos(group1),
           xmax = as.integer(stage_f) + dxpos(group2),
           ytop = max(d[[burden]][d$stage_f == stage_f &
                                    d$dx %in% c(group1, group2)], na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(stage_f) %>%
    arrange(p.adj, .by_group = TRUE) %>%
    mutate(ybar = ytop + ymax * (0.05 + 0.09 * (row_number() - 1))) %>%
    ungroup()

  tick <- ymax * 0.018

  ggplot(d, aes(xpos, .data[[burden]], fill = dx, color = dx)) +
    geom_boxplot(aes(group = interaction(stage_f, dx)),
                 width = 0.22, alpha = 0.35, outlier.shape = NA, linewidth = 0.8) +
    geom_jitter(width = 0.055, size = 2, alpha = 0.6) +
    geom_segment(data = ann, aes(x = xmin, xend = xmax, y = ybar, yend = ybar),
                 inherit.aes = FALSE, linewidth = 0.6) +
    geom_segment(data = ann, aes(x = xmin, xend = xmin,
                                 y = ybar, yend = ybar - tick),
                 inherit.aes = FALSE, linewidth = 0.6) +
    geom_segment(data = ann, aes(x = xmax, xend = xmax,
                                 y = ybar, yend = ybar - tick),
                 inherit.aes = FALSE, linewidth = 0.6) +
    geom_text(data = ann, aes(x = (xmin + xmax) / 2, y = ybar,
                              label = p.adj.signif),
              inherit.aes = FALSE, size = 7, vjust = -0.15) +
    scale_x_continuous(breaks = seq_along(levels(d$stage_f)),
                       labels = levels(d$stage_f)) +
    scale_fill_manual(values = dx_cols, name = "Diagnosis", drop = FALSE) +
    scale_color_manual(values = dx_cols, name = "Diagnosis", drop = FALSE) +
    coord_cartesian(ylim = c(NA, ymax * 1.28)) +
    labs(x = xlab, y = ylab) +
    base_theme + theme(legend.position = "bottom")
}

# ── 5. build and save all six ────────────────────────────────────────────────
for (s in specs) {
  ggsave(file.path(FIG_DIR, sprintf("fig_%s_by_diagnosis.png", s$tag)), plot_a(s$b, s$y),
         width = 8, height = 8, dpi = 600, bg = "white")
  ggsave(file.path(FIG_DIR, sprintf("fig_%s_by_stage.png", s$tag)),
         plot_b(s$b, s$s, s$y, s$x, measure_lab = s$y),
         width = 12, height = 8, dpi = 600, bg = "white")
}
