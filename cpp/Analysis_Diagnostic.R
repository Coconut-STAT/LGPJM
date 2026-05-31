#!/usr/bin/env Rscript
## Analysis_Diagnostic.R — Produces traceplots (PDF) and R-hat (CSV)
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Analysis_Diagnostic.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "analysis_common.R"))
analysis_diagnostic()
