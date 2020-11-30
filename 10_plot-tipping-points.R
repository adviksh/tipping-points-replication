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
tp_orig_raw  <- read_dta(here("data", "tippingsamp_withtps.dta"))
tp_bayes_raw <- read_rds(here("out", "bayes_breakpoints.rds"))


# Process Data ------------------------------------------------------------
message("Processing data...")
msa_bayes <- unique(tp_bayes_raw$msa)

msa_xwalk <- tp_orig_raw %>% 
  transmute(msa,
            msa_name = str_remove(place, " [P]{0,1}MSA")) %>% 
  count(msa, msa_name, sort = TRUE) %>% 
  mutate(msa_name = factor(msa_name, levels = msa_name))


tp_orig <- tp_orig_raw %>% 
  filter(msa %in% msa_bayes) %>% 
  transmute(msa,
            tp_fp    = fpneg_search70_mstar,
            tp_struct = bigt_search70_mstar) %>% 
  unique() %>% 
  pivot_longer(starts_with("tp_"), 
               names_to  = "tip_method",
               values_to = "tipping_point") %>% 
  mutate(tip_method = str_remove(tip_method, "tp_"),
         tip_method = tidy_tip_method(tip_method),
         post_prob = 1)

tp_bayes <- tp_bayes_raw %>% 
  rename(tipping_point = bp_sim,
         post_prob     = prob) %>% 
  mutate(tip_method = "bayes")

tp <- bind_rows(tp_orig, tp_bayes) %>% 
  inner_join(msa_xwalk, by = "msa")
  
# Plot --------------------------------------------------------------------
message("Plotting...")
colors <- colorblind_pal()(3)

subtitle_text <- glue("Comparing ",
                      "**<span style='color:{colors[1]};'>Bayesian structural break model</span>**, ",
                      "**<span style='color:{colors[2]};'>CMR structural break</span>**, and ",
                      "**<span style='color:{colors[3]};'>CMR fixed point</span>**")

tp_plot <- ggplot(tp, aes(x = tipping_point,
               color = tip_method)) +
  facet_wrap(~ msa_name) +
  labs(subtitle = subtitle_text,
       x = "Tipping Point",
       y = "Posterior Probability") +
  theme_bw(base_family = "Lato",
           base_size = 12) +
  theme(plot.subtitle = element_markdown(hjust = 0.5),
        strip.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom") +
  scale_color_manual(values = colors, guide = FALSE) +
  geom_segment(aes(xend = tipping_point,
                   y = 0,
                   yend = post_prob)) +
  geom_point(aes(y = post_prob),
             size = 0.4)


# Save --------------------------------------------------------------------
message("Saving...")
ggsave(here("out", "extension_tipping-points.png"),
       tp_plot,
       width = 12, height = 10, units = "in")
