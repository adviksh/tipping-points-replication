include makefile_helpers


# Target groups ----------------------------------------------------------------
dirs = temp out log

# MSAs to model ----------------------------------------------------------------
msas = 5600 4480 1600 6160 2160 8840 3360 5120 1680 1120
msas += 6280 1920 720 520 5380 7320 6780 6200 5945 7040

msa_models = $(foreach m,$(msas),temp/bp-reg_msa-$(m).rds)

# All --------------------------------------------------------------------------
all: $(dirs)
all: temp/estimates_repl-cmr.rds
all: temp/tipping_samp.rds
all: temp/tipping_points.rds
all: temp/estimates_repl-my-data.rds
all: out/replication.png
all: stan/breakpoint-regression.rds
all: $(msa_models)
all: out/bayes_breakpoints.rds
all: out/bayes_deltas.rds
all: out/bayes_y_hats.rds

# Recipes ----------------------------------------------------------------------
$(dirs):
	mkdir -p $@


temp/estimates_repl-cmr.rds: data/tippingsamp_withtps.dta
temp/estimates_repl-cmr.rds: data/transport_2_20_07.dta
temp/estimates_repl-cmr.rds: 00_replicate_with-cmr-data.R
	$(call r, $<)


temp/tipping_samp.rds: data/msa_def.dta
temp/tipping_samp.rds: data/tract_level_covars_merge_3_25_06.dta
temp/tipping_samp.rds: 01_make-sample.R
	$(call r, $<)
	

temp/tipping_points.rds: temp/tipping_samp.rds
temp/tipping_points.rds: 02_locate-tipping-points.R
	$(call r, $<)

	
temp/estimates_repl-my-data.rds: data/transport_2_20_07.dta
temp/estimates_repl-my-data.rds: temp/tipping_points.rds
temp/estimates_repl-my-data.rds: temp/tipping_points.rds
temp/estimates_repl-my-data.rds: 03_replicate_my-data.R
	$(call r, $<)

	
out/replication.png: data/estimates_orig.csv
out/replication.png: temp/estimates_repl-cmr.rds
out/replication.png: temp/estimates_repl-my-data.rds
out/replication.png: 04_plot-replication.R
	$(call r, $<)
	
	
stan/breakpoint-regression.rds: stan/breakpoint-regression.stan
stan/breakpoint-regression.rds: 05_compile-stan.R
	$(call r, $<)
	

temp/bp-reg_msa-%.rds: 06_fit-stan.R stan/breakpoint-regression.rds
	Rscript --vanilla --verbose 06_fit-stan.R $* > log/06_fit-stan_msa-$*.log 2>&1
	
	
out/bayes_breakpoints.rds: $(msa_models)
out/bayes_breakpoints.rds: 07_extract-breakpoints.R
	$(call r, $<)


out/bayes_deltas.rds: $(msa_models)
out/bayes_deltas.rds: 08_extract-deltas.R
	$(call r, $<)
	
out/bayes_y_hats.rds: $(msa_models)
out/bayes_y_hats.rds: 09_extract-y-hats.R
	$(call r, $<)