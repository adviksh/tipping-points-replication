# Seed --------------------------------------------------------------------
set.seed(150)


# Timing ------------------------------------------------------------------
tictoc::tic()


# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(tidybayes)


# Helpers -----------------------------------------------------------------
get_y_hat <- function(fit) {
  
  draws <- spread_draws(fit, y_hat[x])
  
  # only keep 100 samples
  which_keep <- sample(unique(draws$.draw), size = 100)
    
  draws %>% 
    ungroup() %>% 
    filter(.draw %in% which_keep) %>% 
    select(.draw, x, y_hat)
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


# Extract Deltas ----------------------------------------------------------
message("Extracting deltas...")
fit_tb <- tibble(filename = fit_files,
                 fit      = map(filename, read_rds),
                 msa      = map_int(filename, extract_msa))

y_hat_tb <- fit_tb %>% 
  transmute(msa,
            y_hat_tb = map(fit, get_y_hat)) %>% 
  unnest(y_hat_tb)


# Save --------------------------------------------------------------------
message("Saving...")
write_rds(y_hat_tb, here("out", "bayes_y_hats.rds"))
