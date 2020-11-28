# Timing ------------------------------------------------------------------
tictoc::tic()


# Parallel ----------------------------------------------------------------
options(mc.cores = parallel::detectCores())


# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(rstan)
library(tidybayes)


# Helpers -----------------------------------------------------------------
poly_matrix <- function(x, degree) {
  x_cols <- map(seq(1, degree), ~x^.x)
  Reduce(cbind, x_cols)
}

# Load Data ---------------------------------------------------------------
message("Loading data...")
tp_raw <- read_rds(here("temp", "tipping_samp.rds"))


# Process Data ------------------------------------------------------------
tp <- tp_raw %>% 
  filter(samp == TRUE)

standardize <- function(x) { as.numeric(scale(x, center = TRUE, scale = FALSE)) }

chi70 <- tp %>% 
  filter(msa == 1600, base_year == "70",
         !is.na(frac_minority),
         !is.na(norm_chg_white)) %>% 
  select(msa, base_year, geo2000, frac_minority, norm_chg_white) %>% 
  arrange(frac_minority) %>% 
  group_by(msa, base_year) %>% 
  mutate(frac_minority_std = scale(frac_minority, scale = FALSE),
         norm_chg_white_std = scale(norm_chg_white, scale = FALSE)) %>% 
  ungroup()

# Compile Model -----------------------------------------------------------
stan_mod <- stan_model(file = here("stan", "toy_breakpoint-regression.stan"))


# Toy Data ----------------------------------------------------------------
n_toy <- 1000
toy <- tibble(x = sample(chi70$frac_minority, n_toy),
              mu_y = x + (1 / 1000) * x^2 + (-60) * (x > 15),
              y =  rnorm(n_toy, mean = mu_y, sd = 5)) %>% 
  arrange(x)
plot(toy$x, toy$y)

# Fit Stan on Toy ---------------------------------------------------------
toy_degree <- 2
toy_data <- list(N = n_toy,
                 K = toy_degree,
                 x_mat = poly(toy$x, degree = toy_degree),
                 x = toy$x,
                 y = toy$y)

fit_toy <- sampling(stan_mod, data = toy_data,
                    chains = 4,
                    iter = 2000,
                    pars = c("conditional_mean"),
                    include = FALSE)
fit_toy

lm(y ~ poly(x, 2) + I(x > 15), data = toy)

# Fit Stan on Chi ---------------------------------------------------------
chi_degree <- 4
chi_data <- list(N = nrow(chi70),
                 K = chi_degree + 1,
                 x = poly_matrix(chi70$frac_minority, degree = chi_degree),
                 y = chi70$norm_chg_white)

fit_chi <- sampling(stan_mod, data = chi_data,
                    chains = 8,
                    warmup = 1000,
                    iter = 3000,
                    pars = c("conditional_mean"),
                    include = FALSE)
