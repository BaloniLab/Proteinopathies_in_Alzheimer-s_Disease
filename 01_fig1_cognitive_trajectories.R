# ══════════════════════════════════════════════════════════════════════════
# Cognitive trajectories by diagnosis (Figure 1: global cognition and cognitive domains)

#
#   Control / AsymAD / AD cognitive trajectories (years before death and follow-up year),
#   linear mixed models (diagnosis x time + age at death, random intercept and slope per donor),
#   per-group slopes (emmeans), and the source-data tables for Figures 1-4.
#
# INPUTS
#   R objects: dataset_1744_long_08_30_2026, meta_filtered, dataset_1744_cross_sectional_08_30_2026 (see data/README.md)
# OUTPUTS
#   figures/fig_global_*.png, figures/fig_domains_*.png
#   source_data/Table1_*.csv, Fig1-4_*.csv, FigS_*.csv, sessionInfo.txt
# RUN AFTER
#   none
#
# Run from the project root (see README.md). The analysis scripts keep their
# intermediate objects in the R session; run_all.R sources them in order.
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")
require_packages(c("dplyr", "tidyr", "ggplot2", "conflicted", "lmerTest", "emmeans"))
load_inputs(c("dataset_1744_long_08_30_2026", "meta_filtered", "dataset_1744_cross_sectional_08_30_2026"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(conflicted)   # needed for conflicts_prefer() below


okabe <- c(Control = "#0072B2", AsymAD = "#6A1B9A", AD = "red4")


# ---- 1. data ----------------------------------------------------------------
lng  <- dataset_1744_long_08_30_2026 %>% mutate(projid = as.character(projid))
meta <- meta_filtered                %>% mutate(projid = as.character(projid))

dx_map <- meta %>%
  select(projid, diagnosis = EmoryStrictDx.2019, age_death) %>%
  filter(!is.na(diagnosis)) %>%
  distinct()

dat <- lng %>%
  inner_join(dx_map, by = "projid") %>%
  mutate(diagnosis    = factor(diagnosis, levels = c("Control", "AsymAD", "AD")),
         yrs_to_death = age_at_visit_int - age_death)

# ---- 2. global cognition: years before death --------------------------------
d1 <- dat %>% filter(!is.na(cogn_global), yrs_to_death >= -20, yrs_to_death <= 0)

p1 <- ggplot(d1, aes(yrs_to_death, cogn_global, color = diagnosis, fill = diagnosis)) +
  geom_line(aes(group = projid), alpha = 0.20, linewidth = 0.3) +
  geom_smooth(method = "loess", span = 0.9, linewidth = 1.6, alpha = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = okabe) + scale_fill_manual(values = okabe) +
  coord_cartesian(ylim = c(-2.5, 1.2)) +
  labs(x = "Years before death", y = "Global cognition (z)", color = NULL, fill = NULL) +
  theme_classic(base_size = 20) + theme(legend.position = "top") +
  theme(axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.line = element_line(color = "black", linewidth = 1),
      legend.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"))

p1
# ---- n per diagnosis, folded into the legend --------------------------------
n_dx <- d1 %>% distinct(projid, diagnosis) %>% count(diagnosis)
n_dx

dx_labs <- setNames(paste0(n_dx$diagnosis, " (n = ", n_dx$n, ")"), n_dx$diagnosis)

# ---- slopes from the mixed model, for the annotation -----------------------
library(lmerTest)
conflicts_prefer(lmerTest::lmer)

d1 <- d1 %>% mutate(age_death_c = age_death - mean(age_death, na.rm = TRUE))

m1 <- lmer(cogn_global ~ diagnosis * yrs_to_death + age_death_c +
             (yrs_to_death | projid), data = d1,
           control = lmerControl(optimizer = "bobyqa"))

fe <- summary(m1)$coefficients
b_ctrl <- fe["yrs_to_death", "Estimate"]

star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                                    ifelse(p < 0.05, "*", "ns")))

ann <- data.frame(
  diagnosis = factor(c("Control", "AsymAD", "AD"),
                     levels = c("Control", "AsymAD", "AD")),
  label = c(
    sprintf("beta == %.3f", b_ctrl),
    sprintf("beta == %.3f ~ '%s'",
            b_ctrl + fe["diagnosisAsymAD:yrs_to_death", "Estimate"],
            star(fe["diagnosisAsymAD:yrs_to_death", "Pr(>|t|)"])),
    sprintf("beta == %.3f ~ '%s'",
            b_ctrl + fe["diagnosisAD:yrs_to_death", "Estimate"],
            star(fe["diagnosisAD:yrs_to_death", "Pr(>|t|)"]))),
  x = -19, y = c(-1.3, -1.65, -2.0))

p1_ann <- p1 +
  geom_text(data = ann, aes(x, y, label = label, color = diagnosis),
            parse = TRUE, hjust = 0, size = 8, fontface = "bold",
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_color_manual(values = okabe, labels = dx_labs) +
  scale_fill_manual(values  = okabe, labels = dx_labs)
p1_ann
ggsave(file.path(FIG_DIR, "fig_global_yrs_to_death_stats.png"), p1_ann,
       width = 10, height = 9, dpi = 600, bg = "white")
ggsave(file.path(FIG_DIR, "fig_global_yrs_to_death.png"), p1, width = 10, height = 9, dpi = 600, bg = "white")

# ---- 3. global cognition: follow-up year ------------------------------------
d2 <- dat %>%
  filter(!is.na(cogn_global)) %>%
  group_by(diagnosis, fu_year) %>% filter(n() >= 10) %>% ungroup()

p2 <- ggplot(d2, aes(fu_year, cogn_global, color = diagnosis, fill = diagnosis)) +
  geom_line(aes(group = projid), alpha = 0.10, linewidth = 0.3) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.25, color = NA) +
  stat_summary(fun = mean, geom = "line",  linewidth = 1.6) +
  stat_summary(fun = mean, geom = "point", size = 3.5) +
  scale_color_manual(values = okabe) + scale_fill_manual(values = okabe) +
  coord_cartesian(ylim = c(-2.5, 1.2)) +
  labs(x = "Follow-up year", y = "Global cognition (z)", color = NULL, fill = NULL) +
  theme_classic(base_size = 20) + theme(legend.position = "top")+
  theme(axis.text = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold", color = "black"),
        axis.line = element_line(color = "black", linewidth = 1),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
p2
ggsave(file.path(FIG_DIR, "fig_global_fu_year.png"), p2, width = 10, height = 9, dpi = 600, bg = "white")
# ---- n at baseline and at the last retained visit --------------------------
n_tab <- d2 %>%
  group_by(diagnosis) %>%
  summarise(n_bl   = n_distinct(projid[fu_year == 0]),
            n_last = n_distinct(projid[fu_year == max(fu_year)]),
            yr_max = max(fu_year), .groups = "drop")
n_tab

dx_labs2 <- setNames(paste0(n_tab$diagnosis, " (n = ", n_tab$n_bl, ")"), n_tab$diagnosis)

# donors still contributing at each year, as a strip under the curves
n_strip <- d2 %>% count(diagnosis, fu_year)

p2_ann <- p2 +
  geom_text(data = n_strip,
            aes(fu_year, y = -2.15 - 0.22 * (as.integer(diagnosis) - 1),
                label = n, color = diagnosis, fontface = "bold"),
            size = 4.5, inherit.aes = FALSE, show.legend = FALSE) +
  annotate("text", x = 0, y = -1.85, label = "Donors at risk",
           hjust = 0, size = 4.5, fontface = "bold", color = "grey30") +
  scale_color_manual(values = okabe, labels = dx_labs2) +
  scale_fill_manual(values  = okabe, labels = dx_labs2) +
  coord_cartesian(ylim = c(-2.9, 1.2))
p2_ann
ggsave(file.path(FIG_DIR, "fig_global_fu_year_n.png"), p2_ann, width = 10, height = 9, dpi = 600, bg = "white")
# ---- 4. five domains, long format -------------------------------------------
dom_lab <- c(cogn_global ="Global cognition",  cogn_ep = "Episodic memory",  cogn_po = "Perceptual orientation",
             cogn_ps = "Perceptual speed", cogn_se = "Semantic memory",
             cogn_wo = "Working memory")
conflicts_prefer(base::unname)
dom <- dat %>%
  pivot_longer(all_of(names(dom_lab)), names_to = "domain", values_to = "score") %>%
  filter(!is.na(score)) %>%
  mutate(domain = factor(dom_lab[domain], levels = as.character(dom_lab)))
# ---- 5. domains: years before death -----------------------------------------
d3 <- dom %>% filter(yrs_to_death >= -20, yrs_to_death <= 0)

p3 <- ggplot(d3, aes(yrs_to_death, score, color = diagnosis, fill = diagnosis)) +
  geom_line(aes(group = projid), alpha = 0.08, linewidth = 0.25) +
  geom_smooth(method = "loess", span = 0.9, linewidth = 1.4, alpha = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ domain, nrow = 2) +
  scale_color_manual(values = okabe) + scale_fill_manual(values = okabe) +
  coord_cartesian(ylim = c(-2.5, 1.2)) +
  labs(x = "Years before death", y = "Cognitive score (z)", color = NULL, fill = NULL) +
  theme_classic(base_size = 20) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))+
  theme(axis.text = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold", color = "black"),
        axis.line = element_line(color = "black", linewidth = 1),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
p3
ggsave(file.path(FIG_DIR, "fig_domains_yrs_to_death.png"), p3, width = 16, height = 11, dpi = 600, bg = "white")
library(lmerTest)
conflicts_prefer(lmerTest::lmer)

d3 <- dom %>%
  filter(yrs_to_death >= -20, yrs_to_death <= 0) %>%
  mutate(age_death_c = age_death - mean(age_death, na.rm = TRUE))

star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                                    ifelse(p < 0.05, "*", "ns")))

# ---- one mixed model per domain --------------------------------------------
fit_domain <- function(dd) {
  m  <- lmer(score ~ diagnosis * yrs_to_death + age_death_c +
               (yrs_to_death | projid), data = dd,
             control = lmerControl(optimizer = "bobyqa"))
  fe <- summary(m)$coefficients
  b0 <- fe["yrs_to_death", "Estimate"]
  data.frame(
    diagnosis = factor(c("Control", "AsymAD", "AD"),
                       levels = c("Control", "AsymAD", "AD")),
    slope = c(b0,
              b0 + fe["diagnosisAsymAD:yrs_to_death", "Estimate"],
              b0 + fe["diagnosisAD:yrs_to_death",     "Estimate"]),
    p     = c(NA,
              fe["diagnosisAsymAD:yrs_to_death", "Pr(>|t|)"],
              fe["diagnosisAD:yrs_to_death",     "Pr(>|t|)"]))
}

stats_dom <- d3 %>%
  group_by(domain) %>%
  group_modify(~ fit_domain(.x)) %>%
  ungroup() %>%
  mutate(label = ifelse(is.na(p),
                        sprintf("beta == %.3f", slope),
                        sprintf("beta == %.3f ~ '%s'", slope, star(p))),
         x = -19.5,
         y = -1.35 - 0.32 * (as.integer(diagnosis) - 1))

stats_dom %>% select(domain, diagnosis, slope, p) %>%
  mutate(across(c(slope, p), ~ signif(.x, 3))) %>% print(n = 15)

# ---- n per group ------------------------------------------------------------
n_dom  <- d3 %>% distinct(projid, diagnosis) %>% count(diagnosis)
dx_lab3 <- setNames(paste0(n_dom$diagnosis, " (n = ", n_dom$n, ")"), n_dom$diagnosis)

# ---- annotated figure -------------------------------------------------------
p3_ann <- p3 +
  geom_text(data = stats_dom,
            aes(x, y, label = label, color = diagnosis),
            parse = TRUE, hjust = 0, size = 5, fontface = "bold",
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_color_manual(values = okabe, labels = dx_lab3) +
  scale_fill_manual(values  = okabe, labels = dx_lab3) +
  coord_cartesian(ylim = c(-2.6, 1.2))
p3_ann
ggsave(file.path(FIG_DIR, "fig_domains_yrs_to_death_stats.png"), p3_ann,
       width = 16, height = 11, dpi = 600, bg = "white")
# ---- 6. domains: follow-up year ---------------------------------------------
d4 <- dom %>% group_by(diagnosis, domain, fu_year) %>% filter(n() >= 10) %>% ungroup()

p4 <- ggplot(d4, aes(fu_year, score, color = diagnosis, fill = diagnosis)) +
  geom_line(aes(group = projid), alpha = 0.08, linewidth = 0.25) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.25, color = NA) +
  stat_summary(fun = mean, geom = "line",  linewidth = 1.4) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  facet_wrap(~ domain, nrow = 2) +
  scale_color_manual(values = okabe) + scale_fill_manual(values = okabe) +
  coord_cartesian(ylim = c(-2.5, 1.2)) +
  labs(x = "Follow-up year", y = "Cognitive score (z)", color = NULL, fill = NULL) +
  theme_classic(base_size = 20) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))+
  theme(axis.text = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold", color = "black"),
        axis.line = element_line(color = "black", linewidth = 1),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
p4
ggsave(file.path(FIG_DIR, "fig_domains_fu_year.png"), p4, width = 16, height = 11, dpi = 600, bg = "white")
d4 <- dom %>%
  group_by(diagnosis, domain, fu_year) %>% filter(n() >= 10) %>% ungroup()

# ---- baseline n per group (for the legend) ---------------------------------
n_tab4 <- d4 %>%
  group_by(diagnosis) %>%
  summarise(n_bl = n_distinct(projid[fu_year == 0]), .groups = "drop")
n_tab4

dx_lab4 <- setNames(paste0(n_tab4$diagnosis, " (n = ", n_tab4$n_bl, ")"),
                    n_tab4$diagnosis)

# ---- donors at risk per facet, per year ------------------------------------
n_strip4 <- d4 %>% count(diagnosis, domain, fu_year)

# thin the labels — every year is unreadable in a small facet
n_strip4 <- n_strip4 %>% filter(fu_year %% 3 == 0)

p4_ann <- p4 +
  geom_text(data = n_strip4,
            aes(fu_year, y = -2.15 - 0.24 * (as.integer(diagnosis) - 1),
                label = n, color = diagnosis, fontface = "bold"),
            size = 5, inherit.aes = FALSE, show.legend = FALSE) +
  scale_color_manual(values = okabe, labels = dx_lab4) +
  scale_fill_manual(values  = okabe, labels = dx_lab4) +
  coord_cartesian(ylim = c(-2.9, 1.2))
p4_ann
ggsave(file.path(FIG_DIR, "fig_domains_fu_year_n_with_total.png"), p4_ann,
       width = 16, height = 11, dpi = 600, bg = "white")



# Save data


# =============================================================================
# Source data + statistics tables for Figures 1-4
# =============================================================================
library(emmeans)
outdir <- SRC_DIR   # "source_data" (created in R/00_setup.R)
cs <- dataset_1744_cross_sectional_08_30_2026 %>%
  mutate(projid = as.character(projid))
# ---- Table 1. Cohort description -------------------------------------------
cohort <- dat %>%
  distinct(projid, diagnosis, age_death) %>%
  left_join(cs %>% select(projid, msex, educ, apoe_genotype), by = "projid") %>%
  group_by(diagnosis) %>%
  summarise(n_donors      = n(),
            age_death_mean = round(mean(age_death, na.rm = TRUE), 1),
            age_death_sd   = round(sd(age_death,   na.rm = TRUE), 1),
            female_n       = sum(msex == 0, na.rm = TRUE),
            female_pct     = round(100 * mean(msex == 0, na.rm = TRUE), 1),
            educ_mean      = round(mean(educ, na.rm = TRUE), 1),
            educ_sd        = round(sd(educ,   na.rm = TRUE), 1),
            .groups = "drop") %>%
  left_join(dat %>% filter(!is.na(cogn_global)) %>%
              group_by(diagnosis) %>%
              summarise(n_visits = n(),
                        visits_per_donor_median = median(table(projid)),
                        followup_yrs_median = round(median(tapply(fu_year, projid, max)), 1),
                        .groups = "drop"),
            by = "diagnosis")

write.csv(cohort, file.path(outdir, "Table1_cohort_characteristics.csv"), row.names = FALSE)

# ---- Fig 1. Global cognition ~ years before death ---------------------------
fe1 <- as.data.frame(summary(m1)$coefficients)
fe1 <- data.frame(term = rownames(fe1), fe1, row.names = NULL)
names(fe1) <- c("term","estimate","std_error","df","t_value","p_value")
fe1[-1] <- lapply(fe1[-1], function(x) signif(x, 4))
write.csv(fe1, file.path(outdir, "Fig1_global_LMM_fixed_effects.csv"), row.names = FALSE)

# per-group slopes and all pairwise contrasts
tr1 <- emtrends(m1, ~ diagnosis, var = "yrs_to_death")
write.csv(as.data.frame(summary(tr1, infer = TRUE)),
          file.path(outdir, "Fig1_global_slopes_by_group.csv"), row.names = FALSE)
write.csv(as.data.frame(pairs(tr1, adjust = "tukey")),
          file.path(outdir, "Fig1_global_slope_contrasts.csv"), row.names = FALSE)

# random effects + fit
vc1 <- as.data.frame(VarCorr(m1))
write.csv(vc1, file.path(outdir, "Fig1_global_LMM_random_effects.csv"), row.names = FALSE)

# ---- Fig 2. Global cognition ~ follow-up year (descriptive) -----------------
fig2_src <- d2 %>%
  group_by(diagnosis, fu_year) %>%
  summarise(n_donors = n_distinct(projid),
            mean = round(mean(cogn_global), 4),
            sd   = round(sd(cogn_global), 4),
            se   = round(sd(cogn_global) / sqrt(n()), 4),
            .groups = "drop")
write.csv(fig2_src, file.path(outdir, "Fig2_global_fu_year_summary.csv"), row.names = FALSE)

# ---- Fig 3. Domains ~ years before death ------------------------------------
dom_stats <- stats_dom %>%
  select(domain, diagnosis, slope, p) %>%
  mutate(slope = signif(slope, 4),
         p     = signif(p, 4),
         p_BH  = signif(p.adjust(p, method = "BH"), 4),
         comparison = ifelse(is.na(p), "reference", "vs Control (interaction)"))
write.csv(dom_stats, file.path(outdir, "Fig3_domain_slopes_and_tests.csv"), row.names = FALSE)

# ---- Fig 4. Domains ~ follow-up year (descriptive) --------------------------
fig4_src <- d4 %>%
  group_by(domain, diagnosis, fu_year) %>%
  summarise(n_donors = n_distinct(projid),
            mean = round(mean(score), 4),
            sd   = round(sd(score), 4),
            se   = round(sd(score) / sqrt(n()), 4),
            .groups = "drop")
write.csv(fig4_src, file.path(outdir, "Fig4_domain_fu_year_summary.csv"), row.names = FALSE)

# ---- Donors at risk (both follow-up panels) ---------------------------------
write.csv(d2 %>% count(diagnosis, fu_year, name = "n_visits"),
          file.path(outdir, "FigS_donors_at_risk_global.csv"), row.names = FALSE)
write.csv(d4 %>% count(domain, diagnosis, fu_year, name = "n_visits"),
          file.path(outdir, "FigS_donors_at_risk_domains.csv"), row.names = FALSE)

# ---- AsymAD decline heterogeneity (optional, see note below) ----------------
# NOTE: asym_slopes (per-donor AsymAD slopes) is computed in a separate analysis that is not part of
# this release; the block only runs if that object exists in the session.
if (exists("asym_slopes")) {
  write.csv(asym_slopes %>% mutate(slope = signif(slope, 4)),
            file.path(outdir, "FigS_AsymAD_per_donor_slopes.csv"), row.names = FALSE)
}

# ---- Reproducibility --------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
list.files(outdir)
