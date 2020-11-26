# Timing ------------------------------------------------------------------
tictoc::tic()


# Libraries ---------------------------------------------------------------
library(here)
library(haven)
library(tidyverse)
library(glue)
library(estimatr)
library(broom)


# Helpers -----------------------------------------------------------------
first_char <- function(x) { substr(x, 1, 1) }


# Load Data ---------------------------------------------------------------
message("Loading data...")
tp_raw    <- read_dta(here("data", "tippingsamp_withtps.dta"))
transport <- read_dta(here("data", "transport_2_20_07.dta"))

tp <- tp_raw %>% 
  inner_join(transport, by = "geo2000") %>% 
  mutate(frac_minority_70 = frac_minority_70b,
         fraction_pubtran7 = trvlpb7n / (trvlpb7d + wkhome7),
         fraction_pubtran8 = trvlpb8n / indemp8,
         fraction_pubtran9 = trvlpb9n / indemp9,
         msa = factor(msa))


# Set Up Pooled Regressions ------------------------------------------------
message("Setting up pooled regressions...")
make_pooled_x_names <- function(base_year) {
  stems <- c("unemprt", "lnfaminc", "fr_vac", "fr_rent", "fr_1unit",
             "fraction_pubtran")
  
  controls <- c(glue("d{1:4}"), "past",
                glue("{stems}{first_char(base_year)}"),
                "msa")
  controls
}

make_pooled_y_names <- function(base_year) {
  stems <- c("white", "min", "pop")  
  next_year <- switch(base_year,
                      `70` = "80",
                      `80` = "90",
                      `90` = "00")
  out <- glue("norm_chg_{stems}_{base_year}{next_year}")
  as.character(out)
}

make_pooled_data <- function(tip_method, base_year, x_nm, y_nm, data) {
  
  samp_var <- glue("samp_{base_year}")
  search_var <- glue("searchsamp{base_year}")
  data_yr  <- filter(data,
                     !!sym(samp_var) == 1,
                     !!sym(search_var) == 0)
  
  frac_minority_var <- glue("frac_minority_{base_year}")
  mstar_var <- glue("{tip_method}_search{base_year}_mstar")
  
  
  data_yr[["d1"]] <- data_yr[[frac_minority_var]] - data_yr[[mstar_var]]
  data_yr <- mutate(data_yr,
                    d2 = d1^2,
                    d3 = d1^3,
                    d4 = d1^4,
                    past = as.integer(d1 > 0))
  
  select(data_yr, all_of(y_nm), all_of(x_nm), msa)
}

pooled_config <- crossing(tip_method = c("fpneg", "bigt"),
                         base_year = c("70", "80", "90")) %>% 
  mutate(x_nm  = map(base_year, make_pooled_x_names),
         y_nm  = map(base_year, make_pooled_y_names)) %>% 
  mutate(data = pmap(., make_pooled_data, data = tp))


# Fit Pooled Regressions ---------------------------------------------------
message("Fitting pooled regressions...")
fit_pooled_model <- function(y_nm, x_nm, data) {
  models <- map(y_nm,
                ~lm_robust(reformulate(response = .x,
                                       termlabels = x_nm),
                           data = data,
                           clusters = msa,
                           se_type = "stata"))
  models
}

pooled_models <- pooled_config %>% 
  transmute(tip_method,
            base_year,
            y_nm,
            model = pmap(list(y_nm, x_nm, data), fit_pooled_model)) %>% 
  unnest(c(y_nm, model))


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
  transmute(tip_method, base_year, y_nm,
            model     = map(model, tidy_pooled_model)) %>% 
  unnest(model)


# Set Up By-City Regressions ----------------------------------------------
message("Setting up by-city regressions...")
make_city_x_names <- function(base_year) { c(glue("d{1:4}"), "past") }

make_city_y_names <- function(base_year) {
  next_year <- switch(base_year,
                      `70` = "80",
                      `80` = "90",
                      `90` = "00")
  out <- glue("norm_chg_white_{base_year}{next_year}")
  as.character(out)
}

make_city_data <- function(tip_method, base_year, which_msa, x_nm, y_nm, data) {
  
  samp_var <- glue("samp_{base_year}")
  search_var <- glue("searchsamp{base_year}")
  data_yr  <- filter(data,
                     !!sym(samp_var) == 1,
                     !!sym(search_var) == 0,
                     msa == which_msa)
  
  frac_minority_var <- glue("frac_minority_{base_year}")
  mstar_var <- glue("{tip_method}_search{base_year}_mstar")
  
  
  data_yr <- mutate(data_yr,
                    d1 = !!sym(frac_minority_var),
                    d2 = d1^2,
                    d3 = d1^3,
                    d4 = d1^4,
                    past = as.integer(d1 > !!sym(mstar_var)))
  
  select(data_yr, all_of(y_nm), all_of(x_nm))
}

city_config <- crossing(tip_method = c("fpneg", "bigt"),
                        base_year = c("70", "80", "90"),
                        which_msa = unique(tp$msa)) %>% 
  mutate(x_nm  = map(base_year, make_city_x_names),
         y_nm  = map_chr(base_year, make_city_y_names)) %>% 
  mutate(data = pmap(., make_city_data, data = tp)) %>% 
  filter(map_int(data, nrow) >= 5,
         map_dbl(data, anyNA) == FALSE)


# Fit City Regressions ---------------------------------------------------
message("Fitting city regressions...")
fit_city_model <- function(y_nm, x_nm, data) {
  lm_robust(reformulate(response = y_nm,
                        termlabels = x_nm),
            data = data,
            se_type = "stata")
}

city_models <- city_config %>% 
  transmute(tip_method,
            base_year,
            msa = which_msa,
            y_nm,
            model = pmap(list(y_nm, x_nm, data), fit_city_model)) %>% 
  filter(map_int(model, "rank") == 6)


# Tidy City Models --------------------------------------------------------
message("Processing by-city models...")
tidy_city_model <- function(model) {
  model %>% 
    tidy() %>% 
    filter(term == "past") %>% 
    transmute(estimate = estimate,
              n        = model$nobs) %>% 
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
            model     = map(model, tidy_city_model)) %>% 
  unnest(model) %>% 
  group_by(tip_method, base_year, y_nm) %>% 
  summarize(estimates = list(reweight_city_estimates(estimate, n)),
            .groups = "drop") %>% 
  unnest(estimates)


# Combine Regressions -----------------------------------------------------
message("Combining regressions...")
repl_estimates <- bind_rows(mutate(pooled_estimates, 
                                   regression = "Pooled"),
                            mutate(city_estimates, 
                                   regression = "Fully interacted")) %>% 
  mutate(status = "Replication (CMR Data)")


# Save --------------------------------------------------------------------
message("Saving...")
write_rds(repl_estimates, here("temp", "estimates_repl-cmr.rds"))


# Done --------------------------------------------------------------------
message("Done.")
tictoc::toc()
