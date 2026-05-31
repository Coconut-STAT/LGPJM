#!/usr/bin/env Rscript
## Analysis_VaryingCoefAlpha.R — Produces Table 1 style: RMSE Median (IQR) for 4 settings
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Analysis_VaryingCoefAlpha.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "analysis_common.R"))
analysis_varying_coef_alpha()
