# Timing ------------------------------------------------------------------
tictoc::tic()

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(tidybayes)


# Helpers -----------------------------------------------------------------
summarize_deltas <- function(fit) {
  draws <- spread_draws(fit, delta)
  summarize(draws,
            mean = mean(delta),
            q_05 = quantile(delta, 0.05),
            q_10 = quantile(delta, 0.1),
            q_50 = median(delta),
            q_90 = quantile(delta, 0.9),
            q_95 = quantile(delta, 0.95))
}

extract_msa <- function(filename) {
  msa_slug <- str_extract(filename, "msa-[0-9]+")
  msa <- str_extract(msa_slug, "[0-9]+")
  as.integer(msa)
}

# Load Data ---------------------------------------------------------------
message("Loading data...")
fit_files <- list.files(here("temp"),
                        pattern = "bp-reg_msa-[0-9]+.rds",
                        full.names = TRUE)
fit_files <- head(fit_files, 2)


# Extract Deltas ----------------------------------------------------------
message("Extracting deltas...")
fit_tb <- tibble(filename = fit_files,
                 fit      = map(filename, read_rds),
                 msa      = map_int(filename, extract_msa))

delta_tb <- fit_tb %>% 
  transmute(msa,
            delta_tb = map(fit, summarize_deltas)) %>% 
  unnest(delta_tb)


# Save --------------------------------------------------------------------
message("Saving...")
write_rds(bp_tb, here("out", "bayes_deltas.rds"))