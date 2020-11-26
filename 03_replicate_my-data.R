# Timing ------------------------------------------------------------------
tictoc::tic()


# Libraries ---------------------------------------------------------------
library(here)
library(haven)
library(tidyverse)
library(tidylog)
library(glue)
library(estimatr)
library(broom)


# Load Data ---------------------------------------------------------------
message("Loading data...")
tracts    <- read_rds(here("temp", "tipping_samp.rds"))
tpoints   <- read_rds(here("temp", "tipping_points.rds")) 
transport_raw <- read_dta(here("data", "transport_2_20_07.dta"))


# Process Data ------------------------------------------------------------
transport <- transport_raw %>% 
  transmute(geo2000,
            fraction_pubtran_70 = trvlpb7n / (trvlpb7d + wkhome7),
            fraction_pubtran_80 = trvlpb8n / indemp8,
            fraction_pubtran_90 = trvlpb9n / indemp9) %>% 
  pivot_longer(matches("[789]0$")) %>% 
  mutate(base_year = word(name, 3, sep = "_"),
         name = str_remove(name, "_[789]0")) %>% 
  pivot_wider()

tb <- tracts %>% 
  filter(samp == TRUE) %>% 
  filter(estsamp == TRUE) %>% 
  filter(is.finite(lnfaminc)) %>% 
  inner_join(transport, by = c("geo2000", "base_year")) %>% 
  inner_join(tpoints, by = c("msa", "base_year")) %>% 
  transmute(tip_method = method,
            msa = factor(msa),
            base_year,
            norm_chg_white, norm_chg_min, norm_chg_pop,
            frac_minority,
            past = as.integer(frac_minority > m_star),
            unemprt, lnfaminc, fr_vac, fr_rent, fr_1unit, fraction_pubtran) %>% 
  pivot_longer(starts_with("norm_chg_"),
               names_to = "y_nm",
               values_to = "y") 


# Set Up Pooled Regressions ------------------------------------------------
message("Setting up pooled regressions...")
tb_pooled <- nest(tb,
                  data = c(msa, 
                           y, past,
                           frac_minority, unemprt, lnfaminc, fr_vac, fr_rent, fr_1unit,
                           fraction_pubtran))

# Fit Pooled Regressions ---------------------------------------------------
message("Fitting pooled regressions...")
fit_pooled_model <- function(data) {
  
  controls <- c("poly(frac_minority, 4)", "unemprt", "lnfaminc", "fr_vac",
                "fr_rent", "fr_1unit", "fraction_pubtran")
  form <- reformulate(response = "y",
                      termlabels = c("past", controls))
  lm_robust(form,
            data = data,
            clusters = msa,
            se_type = "stata")
}

pooled_models <- tb_pooled %>% 
  transmute(tip_method,
            base_year,
            y_nm,
            model = map(data, fit_pooled_model))


# Tidy Pooled Models -------------------------------------------------------
message("Post-processing pooled regressions...")
tidy_pooled_model <- function(model) {
  model %>% 
    tidy() %>% 
    filter(term == "past") %>% 
    transmute(estimate = round(estimate, 1),
              se       = round(std.error, 1)) %>% 
    as_tibble()
}

pooled_estimates <- pooled_models %>% 
  transmute(tip_method,
            base_year,
            y_nm,
            model      = map(model, tidy_pooled_model)) %>% 
  unnest(model)


# Set Up By-City Regressions ----------------------------------------------
message("Setting up by-city regressions...")
tb_city <- tb %>% 
  filter(y_nm == "norm_chg_white") %>% 
  select(tip_method, msa, base_year,
         y_nm, y, past, frac_minority) %>% 
  nest(data = c(y, past, frac_minority)) %>% 
  filter(map_int(data, nrow) >= 5) %>% 
  filter(map_dbl(data, anyNA) == FALSE)


# Fit City Regressions ---------------------------------------------------
message("Fitting city regressions...")
fit_city_model <- function(data) {
  lm_robust(y ~ past + poly(frac_minority, 4),
            data = data,
            se_type = "stata")
}

city_models <- tb_city %>% 
  transmute(tip_method,
            base_year,
            msa,
            y_nm,
            model = map(data, fit_city_model)) %>% 
  filter(map_int(model, "rank") == 6)


# Tidy City Models --------------------------------------------------------
message("Processing by-city models...")
tidy_city_model <- function(model) {
  model %>% 
    tidy() %>% 
    filter(term == "past") %>% 
    transmute(estimate = estimate,
              n        = model$N) %>% 
    as_tibble()
}

reweight_city_estimates <- function(estimate, n) {
  
  rw_tb <- data.frame(est = estimate, n = n)
  mod   <- lm_robust(est ~ 1, data = rw_tb, se_type = "stata", weights = n)
  
  tidy_mod <- tidy(mod)
  
  transmute(tidy_mod,
            estimate = estimate,
            se = std.error)
}

city_estimates <- city_models %>% 
  transmute(tip_method, base_year, y_nm,
            model = map(model, tidy_city_model)) %>% 
  unnest(model) %>% 
  group_by(tip_method, base_year, y_nm) %>% 
  summarize(estimates = list(reweight_city_estimates(estimate, n)),
            .groups = "drop") %>% 
  unnest(estimates)


repl_estimates <- bind_rows(mutate(pooled_estimates, 
                                   regression = "Pooled"),
                            mutate(city_estimates, 
                                   regression = "Fully interacted")) %>% 
  mutate(status = "Replication (My Data)")


# Save --------------------------------------------------------------------
message("Saving...")
write_rds(repl_estimates, here("temp", "estimates_repl-my-data.rds"))


# Done --------------------------------------------------------------------
message("Done.")
tictoc::toc()
