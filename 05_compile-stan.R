# Timing ------------------------------------------------------------------
tictoc::tic()

# Libraries ---------------------------------------------------------------
library(here)
library(rstan)


# Compile Model -----------------------------------------------------------
message("Compiling stan model")
model <- stan_model(file = here("stan", "breakpoint-regression.stan"),
                    auto_write = TRUE,
                    save_dso = TRUE,
                    verbose = TRUE)



# Save Model --------------------------------------------------------------
saveRDS(model, here("temp", "breakpoint-regression.rds"))

# Done --------------------------------------------------------------------
message("Done.")
tictoc::toc()
