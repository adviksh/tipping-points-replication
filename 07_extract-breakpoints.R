# Timing ------------------------------------------------------------------
tictoc::tic()

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(tidybayes)


# Helpers -----------------------------------------------------------------
summarize_bps <- function(fit) {
  draws  <- spread_draws(fit, bp_sim)
  counts <- count(draws, bp_sim)
  mutate(counts, prob = n / sum(n))
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


# Extract Breakpoints -----------------------------------------------------
message("Extracting breakpoints...")
fit_tb <- tibble(filename = fit_files,
                 fit      = map(filename, read_rds),
                 msa      = map_int(filename, extract_msa))

bp_tb <- fit_tb %>% 
  transmute(msa,
            bp_tb = map(fit, summarize_bps)) %>% 
  unnest(bp_tb)


# Save --------------------------------------------------------------------
message("Saving...")
write_rds(bp_tb, here("out", "bayes_breakpoints.rds"))
