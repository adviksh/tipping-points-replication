# Timing ------------------------------------------------------------------
tictoc::tic()


# Libraries ---------------------------------------------------------------
library(here)
library(haven)
library(tidyverse)


# Helpers -----------------------------------------------------------------
tidy_tip_method <- function(tip_method) {
  map_chr(tip_method, switch,
          struct = "Structural break method",
          bigt   = "Structural break method",
          fp     = "Fixed point method",
          fpneg  = "Fixed point method")
}

tidy_base_year <- function(base_year) {
  map_chr(base_year, switch,
          `70` = "1970-1980",
          `80` = "1980-1990",
          `90` = "1990-2000")
}

# Load Data ---------------------------------------------------------------
message("Loading data..")
tp_orig_raw  <- read_dta(here("data", "tippingsamp_withtps.dta"))
tp_repl_raw <- read_rds(here("temp", "tipping_points.rds"))


# Process Data ------------------------------------------------------------
message("Processing data...")
tp_orig <- tp_orig_raw %>% 
  select(msa,
         matches("fpneg_search[789]0_mstar"),
         matches("bigt_search[789]0_mstar")) %>% 
  unique() %>% 
  pivot_longer(contains("search"), 
               names_to  = c("method", "base_year"),
               names_pattern = "(.+)_search(.+)_mstar",
               values_to = "tipping_point") %>% 
  mutate(base_year = as.character(base_year),
         method    = tidy_tip_method(method))

tp_repl <- tp_repl_raw %>% 
  rename(tipping_point = m_star) %>% 
  mutate(method    = tidy_tip_method(method))

tp <- inner_join(rename(tp_orig, tipping_point_orig = tipping_point),
                 rename(tp_repl, tipping_point_repl = tipping_point),
                 by = c("method", "msa", "base_year")) %>% 
  mutate(base_year = tidy_base_year(base_year)) %>% 
  filter(!is.na(tipping_point_orig),
         !is.na(tipping_point_repl))


# Plot --------------------------------------------------------------------
message("Plotting...")
tp_plot <- ggplot(tp,
       aes(x = tipping_point_orig, y = tipping_point_repl)) +
  facet_grid(base_year ~ method) +
  labs(x = "CMR Tipping Point",
       y = "Re-estimated Tipping Point") +
  scale_x_continuous(expand = c(0, .05), limits = c(0, 65)) +
  scale_y_continuous(expand = c(0, .05), limits = c(0, 65)) +
  theme_bw(base_family = "Lato",
           base_size = 16) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point() +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x)


# Save --------------------------------------------------------------------
message("Saving...")
ggsave(here("out", "replication_tipping-points.pdf"),
       device = cairo_pdf,
       tp_plot,
       height = 12, width = 10, units = "in")


# Done --------------------------------------------------------------------
message("Done.")

