# ══════════════════════════════════════════════════════════════════════════
# Clinical-diagnosis trajectories to the neuropathological endpoint (alluvial plot)

#
#   Clinical state (NCI / MCI / dementia, from the clinical diagnosis dcfdx) at 15, 10, 5 and 0
#   years before death, flowing into the final Control / AsymAD / AD classification.
#
# INPUTS
#   data/dataset_1744_long_08-30-2026.xlsx
#   data/meta_cc_8reg_n323_2026-09-02.csv
# OUTPUTS
#   figures/trajectory_to_endpoint.png
# RUN AFTER
#   none
#
# Run from the project root (see README.md). The analysis scripts keep their
# intermediate objects in the R session; run_all.R sources them in order.
# ══════════════════════════════════════════════════════════════════════════
source("R/00_setup.R")
require_packages(c("tidyverse", "readxl", "ggalluvial"))

library(tidyverse)
library(readxl)
library(ggalluvial)

# ---------------------------------------------------------------- load ----
L <- read_excel(file.path(DATA_DIR, "dataset_1744_long_08-30-2026.xlsx")) %>%
  mutate(projid = sprintf("%08d", as.numeric(projid)))

meta <- read_csv(file.path(DATA_DIR, "meta_cc_8reg_n323_2026-09-02.csv"), show_col_types = FALSE) %>%
  mutate(projid = sprintf("%08d", as.numeric(projid))) %>%
  select(projid, age_death, dx = EmoryStrictDx.2019)

traj <- L %>%
  inner_join(meta, by = "projid") %>%
  filter(!is.na(dcfdx)) %>%
  mutate(ybd   = age_at_visit_int - age_death,
         state = case_when(dcfdx == 1          ~ "NCI",
                           dcfdx %in% c(2, 3)   ~ "MCI",
                           dcfdx %in% c(4, 5, 6) ~ "Dementia"))

# ------------------------------------------------------- build the lodes ----
axes <- c("-15", "-10", "-5", "0")
ids  <- sort(unique(traj$projid))

obs <- traj %>%
  mutate(bin = as.character(cut(ybd, c(-Inf, -12.5, -7.5, -2.5, 0), labels = axes))) %>%
  filter(!is.na(bin)) %>%
  group_by(projid, bin) %>% slice_max(ybd, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(projid, bin, state)

# bin kept as CHARACTER until the very end — avoids factor-level coercion on bind
time_lodes <- expand_grid(projid = ids, bin = axes) %>%
  mutate(bin = factor(bin, axes)) %>% arrange(projid, bin) %>%
  mutate(bin = as.character(bin)) %>%
  left_join(obs, by = c("projid", "bin")) %>%
  group_by(projid) %>% fill(state, .direction = "down") %>% ungroup() %>%
  mutate(state = replace_na(state, "Not yet enrolled"))

end_lodes <- meta %>% filter(projid %in% ids) %>%
  transmute(projid, bin = "Endpoint", state = dx)

lodes <- bind_rows(time_lodes, end_lodes) %>%
  left_join(meta %>% select(projid, dx), by = "projid") %>%
  mutate(bin = factor(bin, levels = c(axes, "Endpoint"))) %>%
  arrange(projid, bin)

# --------------------------------- flow colour = the stratum it lands in ----
lodes <- lodes %>%
  group_by(projid) %>%
  mutate(flow_fill = lead(state),                      # colour = destination
         flow_fill = coalesce(flow_fill, state)) %>%   # Endpoint has no successor
  ungroup()

lv  <- c("NCI", "MCI", "Dementia", "Control", "AsymAD", "AD", "Not yet enrolled")
pal <- c(NCI = "#0072B2", MCI = "#E69F00", Dementia = "#CC79A7",
         Control = "skyblue3", AsymAD = "#cc99ff", AD = "#ff4d4d",
         "Not yet enrolled" = "grey88")

lodes <- lodes %>%
  mutate(state     = factor(state,     levels = lv),
         flow_fill = factor(flow_fill, levels = lv))

# fail loudly rather than silently drawing grey
stopifnot(!any(is.na(lodes$state)), !any(is.na(lodes$flow_fill)),
          all(levels(lodes$state) %in% names(pal)),
          nrow(lodes) == length(ids) * (length(axes) + 1))

# ---------------------------------------------------------------- plot ----
p <- ggplot(lodes, aes(x = bin, stratum = state, alluvium = projid)) +
  geom_flow(aes(fill = flow_fill, order = flow_fill), aes.flow = "forward",
            alpha = 0.55, width = 0.30, colour = "white", linewidth = 0.12) +
  geom_stratum(aes(fill = state), width = 0.30,
               colour = "grey25", linewidth = 0.35) +
  geom_text(stat = "stratum",
            aes(label = paste0(after_stat(stratum), "\nn = ", after_stat(count))),
            size = 4, fontface = "bold", lineheight = 0.95) +
  scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
  labs(x = "Years before death", y = paste0("Donors (n = ", length(ids), ")")) +
  theme_classic(base_size = 20) +
  theme(legend.position = "top",
        panel.grid.major.y = element_line(colour = "grey93"))
p
ggsave(p, filename = file.path(FIG_DIR, "trajectory_to_endpoint.png"), width = 14, height = 9, dpi = 600)

