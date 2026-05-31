#!/usr/bin/env Rscript
## Analysis_Imputation.R — Produces Table 1 style: RMSE Median (IQR) for 2 sample sizes
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Analysis_Imputation.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "analysis_common.R"))
analysis_imputation()
