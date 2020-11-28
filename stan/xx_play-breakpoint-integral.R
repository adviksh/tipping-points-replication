# Timing ------------------------------------------------------------------
tictoc::tic()


# Parallel ----------------------------------------------------------------
options(mc.cores = parallel::detectCores())


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
  filter(msa == 1600, base_year == "70",
         !is.na(frac_minority),
         !is.na(norm_chg_white)) %>% 
  arrange(frac_minority) %>% 
  select(geo2000, frac_minority, norm_chg_white)


# Compile Model -----------------------------------------------------------
stan_mod <- stan_model(file = here("stan", "toy_breakpoint-integral.stan"))


# Toy Data ----------------------------------------------------------------
n_toy <- 40
toy <- tibble(x = seq(30, 70, length.out = n_toy),
              y = rnorm(n_toy, mean = x + 5 * (x > 50))) %>% 
  arrange(x)
plot(toy$x, toy$y)

# Fit Stan on Toy ---------------------------------------------------------
toy_data <- list(N = n_toy,
                 K = 1,
                 x_mat = matrix(toy$x, nrow = n_toy),
                 x = toy$x,
                 y = toy$y)

fit_toy <- sampling(stan_mod, data = toy_data,
                    chains = 1, iter = 500,
                    par = c("xb"), include = FALSE)

draws_toy <- tidy_draws(fit_toy) %>% 
  mutate(s_x = toy_data$x[s])

draws_toy %>% 
  count(s_x)


# Fit Stan on Chi ---------------------------------------------------------
chi_data <- list(N = nrow(chi70),
                 x = chi70$frac_minority,
                 y = chi70$norm_chg_white)

fit_chi <- sampling(stan_mod, data = chi_data,
                    chains = 12, iter = 2000, warmup = 250,
                    par = c("lp"), include = FALSE)
