#!/usr/bin/env Rscript
## Analysis_ErrorSensitivity.R — Produces Table S4: loading, psi, beta BIAS/RMS/CP
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Analysis_ErrorSensitivity.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "analysis_common.R"))
analysis_error_sensitivity()
