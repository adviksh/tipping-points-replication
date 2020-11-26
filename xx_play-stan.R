# Timing ------------------------------------------------------------------
tictoc::tic()


# Parallel ----------------------------------------------------------------
options(mc.cores = parallel::detectCores() - 1)


# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(rstan)
library(tidybayes)


# Load Data ---------------------------------------------------------------
message("Loading data...")

tp_raw <- read_rds(here("temp", "tipping_samp.rds"))


# Process Data ------------------------------------------------------------
tp <- tp_raw %>% 
  filter(samp == TRUE)

chi70 <- tp %>% 
  filter(msa == 1600, base_year == "70") %>% 
  select(geo2000, frac_minority, norm_chg_white)

ggplot(chi70,
       aes(x = frac_minority, y = norm_chg_white)) +
  theme_bw() +
  geom_point() +
  geom_smooth()


# Toy Data ----------------------------------------------------------------
n_toy <- 1000
toy <- tibble(x = runif(n_toy, min = 0, max = 100),
              y = rnorm(n_toy, mean = x > 50)) %>% 
  arrange(x)


# Fit Stan on Toy ---------------------------------------------------------
toy_data <- list(N = n_toy,
                 x = toy$x,
                 y = toy$y)

stan_mod <- stan_model(file = here("stan", "toy_breakpoint.stan"))

fit_toy <- sampling(stan_mod, data = toy_data,
                    chains = 6, iter = 1000, warmup = 250)

