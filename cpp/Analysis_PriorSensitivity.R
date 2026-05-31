#!/usr/bin/env Rscript
## Analysis_PriorSensitivity.R — Produces Table S5: loading, psi, beta BIAS/RMS/CP
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Analysis_PriorSensitivity.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "analysis_common.R"))
analysis_prior_sensitivity()
