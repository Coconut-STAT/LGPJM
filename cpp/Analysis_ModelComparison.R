#!/usr/bin/env Rscript
## Analysis_ModelComparison.R — Produces Table S3: BIAS/RMS/CP + DIC/WAIC
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Analysis_ModelComparison.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "analysis_common.R"))
analysis_model_comparison()
