#!/usr/bin/env Rscript
## Analysis_Main.R — Produces Table 1: RMSE Median (IQR) for 8 settings
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Analysis_Main.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "analysis_common.R"))
analysis_main()
