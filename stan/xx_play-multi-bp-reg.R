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
  filter(samp == TRUE,
         !is.na(frac_minority),
         !is.na(norm_chg_white))

# LA, Chicago
cities <- tp %>% 
  filter(base_year == "70",
         msa %in% c(4480, 1600)) %>% 
  select(msa, base_year, geo2000, frac_minority, norm_chg_white) %>% 
  arrange(base_year, msa, frac_minority)


# Compile Model -----------------------------------------------------------
stan_mod <- stan_model(file = here("stan", "toy_multi-bp-reg.stan"))


# Data Processing Helpers -------------------------------------------------
first_below <- function(x, bs) { map_int(bs, ~min(which(x < .x))) }
last_below  <- function(x, bs) { map_int(bs, ~max(which(x < .x))) }
first_above <- function(x, bs) { map_int(bs, ~min(which(x >= .x))) }
last_above  <- function(x, bs) { map_int(bs, ~max(which(x >= .x))) }

prep_data <- function(tb, K = 2, bp = seq(1, 21, by = 5), nu = 1) {
  
  city_ids <- unique(tb$msa)
  
  x_poly <- tb %>%
    group_by(msa) %>% 
    summarize(x_poly = list(poly(frac_minority, K)), .groups = "drop") %>% 
    pull(x_poly) %>% 
    map(scale) %>% 
    reduce(rbind)
  
  tb_nested <- tb %>% 
    group_by(msa) %>% 
    summarize(fm = list(frac_minority), .groups = "drop") %>% 
    mutate(start_blw = map(fm, first_below, bs = bp),
           end_blw   = map(fm, last_below, bs = bp),
           start_abv = map(fm, first_above, bs = bp),
           end_abv   = map(fm, last_above, bs = bp))
  
  list(N_obs    = nrow(tb),
       N_cit    = length(city_ids),
       N_tract  = table(tb$msa),
       y        = tb$norm_chg_white,
       K        = K,
       x_poly   = x_poly,
       B        = length(bp),
       bp       = bp,
       ct_start_blw = tb_nested$start_blw,
       ct_end_blw   = tb_nested$end_blw,
       ct_start_abv = tb_nested$start_abv,
       ct_end_abv   = tb_nested$end_abv,
       nu       = 1)
}



# Toy Data ----------------------------------------------------------------
toy <- crossing(frac_minority = seq(0, 100, by = 1),
                msa = unique(cities$msa)) %>% 
  arrange(msa, frac_minority) %>% 
  mutate(mu = frac_minority 
         + (-5) * (frac_minority > 15),
         norm_chg_white = rnorm(n = nrow(.), mu))


# Prep Toy Data -----------------------------------------------------------
toy_data <- prep_data(toy, K = 1, bp = seq(5, 25, by = 1))
toy_fit  <- sampling(stan_mod, data = toy_data,
                     chains = 4,
                     iter = 1000) 


# Prep Cities Data --------------------------------------------------------
cities_sub <- cities %>% 
  group_by(base_year, msa) %>% 
  slice_sample(n = 1000) %>% 
  arrange(msa, frac_minority) %>% 
  ungroup()

cities_data <- prep_data(cities_sub,
                         K = 1,
                         bp = seq(3, 31, by = 1))

stan_fit <- sampling(stan_mod, data = cities_data,
                     chains = 4,
                     iter = 2000)
