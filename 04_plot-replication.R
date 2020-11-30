# Timing ------------------------------------------------------------------
tictoc::tic()


# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(ggtext)
library(ggthemes)
library(glue)

# Load Data ---------------------------------------------------------------
est_orig     <- read_csv(here("data", "estimates_orig.csv"))
est_repl_cmr <- read_rds(here("temp", "estimates_repl-cmr.rds"))
est_repl_my  <- read_rds(here("temp", "estimates_repl-my-data.rds"))


# Process Data ------------------------------------------------------------
tidy_tip_method <- function(tip_method) {
  map_chr(tip_method, switch,
          struct = "structural break method",
          bigt   = "structural break method",
          fp     = "fixed point method",
          fpneg  = "fixed point method")
}

tidy_base_year <- function(base_year) {
  map_chr(base_year, switch,
          `70` = "1970-1980",
          `80` = "1980-1990",
          `90` = "1990-2000")
}

tidy_y_nm <- function(y_nm) {
  case_when(str_detect(y_nm, "white") ~ "White pop.",
            str_detect(y_nm, "min")   ~ "Minority pop.",
            str_detect(y_nm, "pop")   ~ "Total pop.")
}

tidy_regression <- function(regression) {
  map_chr(regression, switch,
          `Fully interacted` = "Fully interacted\nregression",
          `Pooled` = "Pooled\nregression")
}

estimates <- est_orig %>% 
  mutate(status = "Original",
         base_year = as.character(base_year)) %>% 
  bind_rows(est_repl_cmr, est_repl_my) %>% 
  mutate(tip_method = tidy_tip_method(tip_method),
         base_year  = tidy_base_year(base_year),
         y_nm       = tidy_y_nm(y_nm),
         regression = tidy_regression(regression)) %>% 
  mutate(y_nm = factor(y_nm,
                       levels = c("White pop.",
                                  "Minority pop.",
                                  "Total pop.")),
         status = factor(status, levels = c("Original",
                                            "Replication (CMR Data)",
                                            "Replication (My Data)")))

# Plot --------------------------------------------------------------------
message("Plotting...")

dodge_width <- 0.4
colors <- colorblind_pal()(3)

subtitle_text <- glue("Comparing ",
                      "**<span style='color:{colors[1]};'>CMR (2008)</span>**, ",
                      "**<span style='color:{colors[2]};'>Replication (CMR Tipping Points)</span>**, and ",
                      "**<span style='color:{colors[3]};'>Replication (My Tipping Points)</span>**")

replication_plot <- ggplot(estimates,
                           aes(x = regression,
                               y = estimate,
                               ymin = estimate - 1.96 * se,
                               ymax = estimate + 1.96 * se,
                               shape = tip_method,
                               color = status)) +
  facet_grid(y_nm ~ base_year,
             scales = "free") +
  labs(subtitle = subtitle_text,
       shape    = "Tipping point estimated by...",
       y = "Population change (pp)") +
  theme_bw(base_size = 16,
           base_family = "Lato") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title    = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.position = "bottom",
        strip.background = element_blank()) +
  scale_color_colorblind(guide = FALSE) +
  guides(shape = guide_legend(direction = "horizontal", nrow = 2)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed") +
  geom_linerange(position = position_dodge(width = dodge_width)) +
  geom_point(position = position_dodge(width = dodge_width),
             size = 2)


# Save --------------------------------------------------------------------
message("Saving...")
ggsave(here("out", "replication.png"), replication_plot,
       width = 10, height = 7, units = "in")

