# Seed --------------------------------------------------------------------
set.seed(799002)


# Timing ------------------------------------------------------------------
tictoc::tic()


# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(tidylog)
library(haven)
library(glue)


# Load Data ---------------------------------------------------------------
message("Loading data...")
tract_raw <- read_dta(here("data", "tract_level_covars_merge_3_25_06.dta"))
msa_names <- read_dta(here("data", "msa_def.dta"))


# Name Tracts and Restrict to MSA -----------------------------------------
message("Joining MSA names to tracts")
safe_left_join <- function(x, y, by) {
  xy <- left_join(x, y, by)
  if (nrow(xy) > nrow(x)) { 
    err_msg <- glue("Join increases the number of rows from",
                    "{nrow(x)} to {nrow(xy)}")
    stop(err_msg)
  }
  xy
}

tract_named <- tract_raw %>% 
  rename(msar = msapma99) %>% 
  safe_left_join(msa_names, by = "msar") %>% 
  filter(msar != 9999) %>% 
  filter(!is.na(msar))


# Race Composition --------------------------------------------------------
message("Calculating racial composition...")
tract_races_1 <- tract_named %>% 
  transmute(geo2000,
            msa = msar,
            msa_name = place,
            # Fraction of whites each year
            # shrnhw*n is "number of non hispanic whites in year *"
            frac_white_00 = ifelse(trctpop0 == 0, NA, shrnhw0n / trctpop0),
            frac_white_90 = ifelse(trctpop9 == 0, NA, shrnhw9n / trctpop9),
            frac_white_80 = ifelse(trctpop8 == 0, NA, shrnhw8n / trctpop8),
            num_white_00 = shrnhw0n,
            num_white_90 = shrnhw9n,
            num_white_80 = shrnhw8n,
            # Different approach for 1970
            frblk7 = ifelse(trctpop7 == 0, NA, shrblk7n / trctpop7),
            frhis7 = ifelse(trctpop7 == 0, NA, shrhsp7n / trctpop7),
            frwht7 = ifelse(trctpop7 == 0, NA, shrwht7n / trctpop7),
            frblk8 = ifelse(trctpop8 == 0, NA, shrblk8n / trctpop8),
            frhis8 = ifelse(trctpop8 == 0, NA, shrhsp8n / trctpop8),
            frwht8 = ifelse(trctpop8 == 0, NA, shrwht8n / trctpop8),
            frnhw8 = ifelse(trctpop8 == 0, NA, shrnhw8n / trctpop8),
            across(starts_with("favinc"), log),
            fr_vac7 = vachu7 / tothsun7,
            fr_vac8 = vachu8 / tothsun8,
            fr_vac9 = vachu9 / tothsun9,
            fr_rent7 = rntocc7 / occhu7,
            fr_rent8 = rntocc8 / occhu8,
            fr_rent9 = rntocc9 / occhu9,
            fr_1unit7 = (ocunit17 + ocunit27) / occhu7,
            fr_1unit8 = (ocunit18 + ocunit28) / occhu8,
            fr_1unit9 = (ocunit19 + ocunit29) / occhu9,
            across(starts_with("trctpop")),
            across(starts_with("unemprt")))

frnhw8_model <- lm(frnhw ~ frblk + frhis + frwht,
                   data = rename_with(tract_races_1,
                                      .fn = str_remove,
                                      .cols = ends_with("8"),
                                      "8"),
                   model = FALSE, y = FALSE)

cap_probability <- function(x) { pmin(pmax(x, 0), 1) }

tract_races_2 <- tract_races_1 %>% 
  mutate(frac_white_80b = predict(frnhw8_model,
                                  data.frame(frblk = frblk8,
                                             frhis = frhis8,
                                             frwht = frwht8)),
         frac_white_70b = predict(frnhw8_model,
                                  data.frame(frblk = frblk7,
                                             frhis = frhis7,
                                             frwht = frwht7)),
         frac_white_80b = cap_probability(frac_white_80b),
         frac_white_70b = cap_probability(frac_white_70b),
         num_white_80b  = frac_white_80b * trctpop8,
         num_white_70b  = frac_white_70b * trctpop7,
         frac_minority_00  = (1 - num_white_00 / trctpop0) * 100,
         frac_minority_90  = (1 - num_white_90 / trctpop9) * 100,
         frac_minority_80  = (1 - num_white_80 / trctpop8) * 100,
         frac_minority_80b = (1 - frac_white_80b) * 100,
         frac_minority_70b = (1 - frac_white_70b) * 100,
         frac_minority_70  = frac_minority_70b) %>% 
  group_by(msa) %>% 
  mutate(frac_minority_msa_00  = mean(frac_minority_00, na.rm = T),
         frac_minority_msa_90  = mean(frac_minority_90, na.rm = T),
         frac_minority_msa_80  = mean(frac_minority_80, na.rm = T),
         frac_minority_msa_80b = mean(frac_minority_80b, na.rm = T),
         frac_minority_msa_70b = mean(frac_minority_70b, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(norm_chg_white_9000 = ((num_white_00-num_white_90) / (trctpop9)) * 100,
         norm_chg_white_8090 = ((num_white_90-num_white_80) / (trctpop8)) * 100,
         norm_chg_white_7080 = ((num_white_80b-num_white_70b) / (trctpop7)) * 100,
         norm_chg_min_9000   = (((trctpop0 - num_white_00) - (trctpop9 - num_white_90)) / (trctpop9)) * 100,
         norm_chg_min_8090   = (((trctpop9 - num_white_90) - (trctpop8 - num_white_80)) / (trctpop8)) * 100,
         norm_chg_min_7080   = (((trctpop8 - num_white_80b) - (trctpop7 - num_white_70b)) / (trctpop7)) * 100,
         norm_chg_pop_9000   = ((trctpop0-trctpop9)/trctpop9) * 100,
         norm_chg_pop_8090   = ((trctpop9-trctpop8)/trctpop8) * 100,
         norm_chg_pop_7080   = ((trctpop8-trctpop7)/trctpop7) * 100,
         pct_pop90 = (trctpop0 - trctpop9) / trctpop9,
         pct_pop80 = (trctpop9 - trctpop8) / trctpop8,
         pct_pop70 = (trctpop8 - trctpop7) / trctpop7)


# Sample Restrictions -----------------------------------------------------
message("Identifying sample restrictions...")
drop_if_any <- function(...) {
  d <- reduce(list(...), `|`)
  d[is.na(d)] <- FALSE
  d
}

tract_samp <- tract_races_2 %>% 
  group_by(msa) %>% 
  mutate(drop_growth_90 = as.logical(scale(pct_pop90) > 5),
         drop_growth_80 = as.logical(scale(pct_pop80) > 5),
         drop_growth_70 = as.logical(scale(pct_pop70) > 5),
         drop_small_90  = trctpop9 < 200,
         drop_small_80  = trctpop8 < 200,
         drop_small_70  = trctpop7 < 200,
         drop_xtrm_90   = norm_chg_white_9000 > 500,
         drop_xtrm_80   = norm_chg_white_8090 > 500,
         drop_xtrm_70   = norm_chg_white_7080 > 500) %>% 
  ungroup() %>% 
  mutate(drop_70 = drop_if_any(drop_growth_70, drop_small_70, drop_xtrm_70),
         drop_80 = drop_if_any(drop_growth_80, drop_small_80, drop_xtrm_80),
         drop_90 = drop_if_any(drop_growth_90, drop_small_90, drop_xtrm_90),
         samp_70 = !drop_70,
         samp_80 = !drop_80,
         samp_90 = !drop_90) %>% 
  group_by(msa) %>% 
  mutate(drop_ntrct_70 = sum(samp_70) < 100,
         drop_ntrct_80 = sum(samp_80) < 100,
         drop_ntrct_90 = sum(samp_90) < 100,
         samp_70 = samp_70 & !drop_ntrct_70,
         samp_80 = samp_80 & !drop_ntrct_80,
         samp_90 = samp_90 & !drop_ntrct_90) %>% 
  ungroup()

tract_clean <- tract_samp %>% 
  select(-starts_with("pct_pop"))

tract_tall <- tract_clean %>% 
  rename_with(str_replace, ends_with("7"), pattern = "7", replacement = "_70") %>% 
  rename_with(str_replace, ends_with("8"), pattern = "8", replacement = "_80") %>% 
  rename_with(str_replace, ends_with("9"), pattern = "9", replacement = "_90") %>% 
  rename_with(str_replace, matches("[a-z]0$"), pattern = "0", replacement = "_00") %>% 
  pivot_longer(matches(".+_[0-9]+b{0,1}"))
  
tract_wide <- tract_tall %>% 
  mutate(base_year = str_extract(name, "_[789]0"),
         base_year = str_remove(base_year, "_"),
         name      = str_remove(name, "[789][0-9]+"),
         name      = str_remove(name, "_$")) %>% 
  filter(is.na(base_year) == FALSE) %>% 
  pivot_wider(names_from = "name",
              values_from = "value") %>% 
  mutate(across(starts_with("drop"), as.logical),
         samp = as.logical(samp)) %>% 
  rename(lnfaminc = favinc)


# Sample Split ------------------------------------------------------------
message("Flagging training and test samples...")
flag_search <- function(samp) {
  
  search  <- rep(FALSE, length(samp))
  n_samp <- sum(samp)
  
  search[sample(which(samp), ceiling(n_samp * 2 / 3))] <- TRUE
  search
}

tract_split <- tract_wide %>% 
  group_by(base_year, msa) %>% 
  mutate(searchsamp = flag_search(samp),
         estsamp    = samp == TRUE & searchsamp == FALSE) %>% 
  ungroup()


# Save --------------------------------------------------------------------
message("Saving...")
write_rds(tract_split,
          here("temp", "tipping_samp.rds"))


# Timing ------------------------------------------------------------------
message("Done.")
tictoc::toc()
