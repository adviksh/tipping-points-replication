# Timing ------------------------------------------------------------------
tictoc::tic()


# Libraries ---------------------------------------------------------------
library(here)
library(haven)
library(tidyverse)
library(ggdist)
library(ggthemes)
library(ggtext)
library(glue)


# Helpers -----------------------------------------------------------------
tidy_tip_method <- function(tip_method) {
  map_chr(tip_method, switch,
          struct = "structural break method",
          bigt   = "structural break method",
          fp     = "fixed point method",
          fpneg  = "fixed point method")
}


# Load Data ---------------------------------------------------------------
message("Loading data...")
obs_orig_raw    <- read_dta(here("data", "tippingsamp_withtps.dta"))
delta_orig_raw  <- read_rds(here("temp", "estimates_delta-cmr.rds"))
delta_bayes_raw <- read_rds(here("out", "bayes_deltas.rds"))


# Process Data ------------------------------------------------------------
message("Processing data...")
msa_bayes <- unique(delta_bayes_raw$msa)

msa_xwalk <- obs_orig_raw %>% 
  transmute(msa,
            msa_name = str_remove(place, " [P]{0,1}MSA")) %>% 
  count(msa, msa_name, sort = TRUE) %>% 
  mutate(msa_name = factor(msa_name, levels = msa_name))


delta_orig <- delta_orig_raw %>% 
  filter(msa %in% msa_bayes,
         y_nm == "norm_chg_white_7080") %>% 
  transmute(msa,
            tip_method = tidy_tip_method(tip_method),
            mean       = estimate)

delta_bayes <- delta_bayes_raw %>% 
  mutate(tip_method = "bayes")

delta <- bind_rows(delta_orig, delta_bayes) %>% 
  inner_join(msa_xwalk, by = "msa")

# Plot --------------------------------------------------------------------
message("Plotting...")
colors <- colorblind_pal()(3)

subtitle_text <- glue("Comparing ",
                      "**<span style='color:{colors[1]};'>Bayesian structural break model</span>**, ",
                      "**<span style='color:{colors[2]};'>CMR structural break</span>**, and ",
                      "**<span style='color:{colors[3]};'>CMR fixed point</span>**")
dodge_width <- 0.25

delta_plot <- ggplot(delta,
                     aes(x = msa_name,
                         color = tip_method)) +
  coord_flip() +
  labs(subtitle = subtitle_text,
       y = "Change in white population share (pp) at tipping point") +
  theme_bw(base_family = "Lato",
           base_size = 12) +
  theme(plot.subtitle = element_markdown(hjust = 0.5),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_manual(values = colors, guide = FALSE) +
  scale_x_discrete(limits = rev) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "darkgrey") +
  geom_point(aes(y = mean),
             position = position_dodge(width = dodge_width)) +
  geom_linerange(aes(ymin = q_10, ymax = q_90),
                 position = position_dodge(width = dodge_width))


# Save --------------------------------------------------------------------
message("Saving...")
ggsave(here("out", "extension_deltas.png"),
       delta_plot,
       height = 12, width = 10, units = "in")
