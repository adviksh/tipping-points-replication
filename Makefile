include makefile_helpers


# Target groups ----------------------------------------------------------------
dirs = temp out log


# All --------------------------------------------------------------------------
all: $(dirs)
all: temp/estimates_repl-cmr.rds
all: temp/tipping_samp.rds
all: temp/tipping_points.rds
all: temp/estimates_repl-my-data.rds
all: out/replication.png


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