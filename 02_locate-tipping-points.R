# Seed --------------------------------------------------------------------
set.seed(799002)


# Timing ------------------------------------------------------------------
tictoc::tic()


# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(furrr)
library(tidylog)
library(glue)


# Load Data ---------------------------------------------------------------
message("Loading data...")
tp_raw <- read_rds(here("temp", "tipping_samp.rds"))


# Trim Columns ------------------------------------------------------------
message("Preprocessing...")
tp_nested <- tp_raw %>% 
  filter(samp == TRUE) %>% 
  filter(searchsamp == TRUE) %>% 
  filter(frac_minority < 60) %>% 
  filter(!is.na(norm_chg_white)) %>% 
  select(msa, base_year, frac_minority, norm_chg_white) %>% 
  nest(data = c(frac_minority, norm_chg_white)) %>% 
  mutate(data = map(data, arrange, frac_minority))


# Structural Break Method -------------------------------------------------
message("Estimating breaks via structural method...")
break_ts <- function(data) {
  
  x <- data$frac_minority
  y <- data$norm_chg_white
  
  if (is.unsorted(x)) stop("x must be sorted")
  
  candidates <- head(x, -1)
  t_stats    <- map_dbl(candidates, break_t,
                        x = x,
                        y = y)
  
  tibble(breakpoint = candidates,
         t_stat     = t_stats)
}

break_t <- function(s, x, y){
  z <- x > s
  m <- lm(y ~ z, model = FALSE, y = FALSE)
  abs(broom::tidy(m)$statistic[2])
}

plan(multisession, workers = 6)
tp_candidates_break <- tp_nested %>% 
  transmute(base_year, msa,
            break_tb = future_map(data, break_ts))

tp_break <- tp_candidates_break %>% 
  transmute(method = "struct",
            msa,
            base_year,
            opt_break = map(break_tb, dplyr::top_n, 1, t_stat)) %>% 
  unnest(opt_break) %>% 
  select(-t_stat) %>% 
  rename(m_star = breakpoint)

# Fixed Point Method ------------------------------------------------------
message("Estimating breaks via fixed point method...")
fixed_point_breaks <- function(data) {
  
  x <- data$frac_minority
  y <- data$norm_chg_white
  
  if (is.unsorted(x)) stop("x must be sorted")
  
  if (length(x) <= 4) return(NA)
  
  first_pass_roots <- roots_quartic(x, y)
  
  if (nrow(first_pass_roots) == 0) return(NA)
  
  first_neg_slope <- dplyr::slice_min(first_pass_roots, slope)
  first_pass_x <- first_neg_slope$x
  
  mask <- abs(x - first_pass_x) <= 10.000001
  
  if (sum(mask) <= 4) return(NA)
  
  sec_pass_roots <- roots_quartic(x[mask], y[mask])
  
  if (nrow(sec_pass_roots) == 0) return(NA)
  
  sec_neg_slope <- dplyr::slice_min(sec_pass_roots, slope)
  sec_neg_slope$x
}

roots_quartic <- function(x, y) {
  
  model <- lm(y ~ poly(x, 4), y = FALSE, model = FALSE)
  
  fitted_tb <- tibble(x       = x,
                      y       = y,
                      y_hat   = fitted(model),
                      slope   = estimate_slope_forward(x, y_hat),
                      y_cent  = y_hat - mean(y),
                      crosses = sign(y_cent) != lead(sign(y_cent)))  
  
  dplyr::filter(fitted_tb,
                crosses == TRUE,
                lead(x) - x < 10)
}

estimate_slope <- function(x, y) {
  
  forward  <- estimate_slope_forward(x, y)
  backward <- estimate_slope_backward(x, y)
  
  average <- (forward + backward) / 2
  
  # Return the average, except where it's undefined
  coalesce(average, forward, backward)
}

estimate_slope_backward <- function(x, y) { (y - lag(y)) / (x - lag(x)) }

estimate_slope_forward  <- function(x, y) { (lead(y) - y) / (lead(x) - x) }

tp_fixed <- tp_nested %>% 
  transmute(method = "fp",
            msa, 
            base_year,
            m_star = map_dbl(data, fixed_point_breaks))


# Combine Estimates -------------------------------------------------------
message("Combining estimates...")
tp <- bind_rows(tp_break, tp_fixed)


# Save --------------------------------------------------------------------
message("Saving...")
write_rds(tp, here("temp", "tipping_points.rds"))


# Done --------------------------------------------------------------------
message("Done.")
tictoc::toc()