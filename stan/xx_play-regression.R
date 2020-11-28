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
  select(geo2000, frac_minority, norm_chg_white) %>% 
  arrange(frac_minority)


# Compile Model -----------------------------------------------------------
stan_mod <- stan_model(file = here("stan", "toy_regression.stan"))


# Toy Data ----------------------------------------------------------------
n_toy <- 1000
toy <- tibble(x = rnorm(n_toy),
              y = rnorm(n_toy, mean = x + 0.5 * x^2 - 0.3 * x^3 + 0.1 * x^4)) %>% 
  arrange(x)


# Fit Stan on Toy ---------------------------------------------------------
toy_data <- list(N = n_toy,
                 K = 4,
                 nu = 1,
                 x = cbind(toy$x, toy$x^2, toy$x^3, toy$x^4),
                 y = toy$y)

fit_toy <- sampling(stan_mod, data = toy_data,
                    chains = 6, iter = 1000, warmup = 250)

plot(fit_toy, pars = c("alpha", "beta"))

draws_toy <- tidy_draws(fit_toy) %>% 
  mutate(s_x = toy_data$x[s])

draws_toy %>% 
  count(s_x)


# Fit Stan on Chi ---------------------------------------------------------
poly_matrix <- function(x, degree) {
  x_cols <- map(seq_len(degree), ~x^.x)
  Reduce(cbind, x_cols)
}

chi_data <- list(N = nrow(chi70),
                 K = 4,
                 nu = 1,
                 x = poly_matrix(scale(chi70$frac_minority), degree = 4),
                 y = as.numeric(scale(chi70$norm_chg_white)))

fit_chi <- sampling(stan_mod, data = chi_data,
                    chains = 12, iter = 1000, warmup = 250)
plot(fit_chi, pars = c("alpha", "beta"))

draws_chi <- gather_draws(fit_chi, "y_hat")

draws_chi <- tidy_draws(fit_chi) %>% 
  rename_with(.cols = everything(), str_replace, "\\[", "_") %>% 
  rename_with(.cols = everything(), str_remove, "\\]")


chi_rows <- floor(seq(from = 1, to = chi_data$N, length.out = 10))

sample_fits_temp <- draws_chi %>% 
  slice_sample(n = 10) %>% 
  transmute(draw = .draw,
            alpha,
            beta = pmap(list(beta_1, beta_2, beta_3, beta_4), c),
            xb   = map(beta, tcrossprod, chi_data$x[chi_rows,]),
            axb  = map2(alpha, xb, `+`))

sample_fits <- sample_fits_temp %>% 
  select(draw, axb) %>% 
  mutate(x   = map(draw, ~chi70$frac_minority[chi_rows]),
         axb = map(axb, as.numeric)) %>% 
  unnest(c(x, axb))

ggplot(sample_fits,
       aes(x = x, y = axb, group = draw)) +
  theme_bw() +
  geom_line()
