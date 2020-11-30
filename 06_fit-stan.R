# Timing ------------------------------------------------------------------
tictoc::tic()


# Parallel ----------------------------------------------------------------
options(mc.cores = min(parallel::detectCores(), 8))


# Libraries ---------------------------------------------------------------
library(here)
library(glue)
library(tidyverse)
library(rstan)


# Command Line Args -------------------------------------------------------
message("Reading command line arguments...")
which_msa <- as.integer(commandArgs(trailingOnly = TRUE))


# Load Data ---------------------------------------------------------------
message("Loading data...")
tp_raw <- read_rds(here("temp", "tipping_samp.rds"))


# Process Data ------------------------------------------------------------
message("Subsetting to MSA: ", which_msa)
city <- tp_raw %>% 
  filter(samp == TRUE,
         base_year == "70",
         msa == which_msa,
         !is.na(frac_minority),
         !is.na(norm_chg_white)) %>% 
  select(frac_minority, norm_chg_white) %>% 
  arrange(frac_minority)

K <- 4
bp <- seq(2, 50, by = 0.5)
x_poly <- poly(city$frac_minority, degree = K)
x_pred <- seq(0, 100)
x_pred_mat <- predict(x_poly, x_pred)

city_data <- list(N = nrow(city),
                  K = K,
                  B = length(bp),
                  x_mat = x_poly,
                  bp = bp,
                  x = city$frac_minority,
                  y = city$norm_chg_white,
                  N_pred = length(x_pred),
                  x_pred = x_pred,
                  x_pred_mat = x_pred_mat)


# Load Stan Model ---------------------------------------------------------
message("Loading stan model...")
model <- read_rds(here("stan", "breakpoint-regression.rds"))


# Fit ---------------------------------------------------------------------
message("Fitting stan model...")
fit <- sampling(model, data = city_data,
                chains = 8,
                warmup = 1000,
                iter = 2000,
                pars = c("xb"),
                include = FALSE,
                seed = 20310)

message("There were ",
        rstan::get_num_divergent(fit),
        " divergent iterations.")


# Save --------------------------------------------------------------------
message("Saving fit...")
write_rds(fit, here("temp", glue("bp-reg_msa-{which_msa}.rds")))


# Done --------------------------------------------------------------------
message("Done.")
tictoc::toc()
