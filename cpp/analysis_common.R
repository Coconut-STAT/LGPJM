## analysis_common.R — Reads <SimName>_<i>.RData files, produces paper-specific CSV/PDF outputs

resolve_repo_root <- function() {
  argv_file <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(argv_file) == 0) argv_file <- "--file=analysis_common.R"
  this_file <- sub("^--file=", "", argv_file[[1]])
  normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = FALSE)
}

load_simulation_results <- function(repo_root, simulation_name, expected_count = NULL) {
  sim_dir <- file.path(repo_root, "output", simulation_name)
  if (!dir.exists(sim_dir)) stop(sprintf("No output directory: %s", sim_dir))
  files <- sort(list.files(sim_dir, pattern = sprintf("^%s_\\d+\\.RData$", simulation_name), full.names = TRUE))
  if (length(files) == 0) stop(sprintf("No RData files found in %s", sim_dir))
  if (!is.null(expected_count) && length(files) < expected_count)
    cat(sprintf("WARNING: expected %d settings but found only %d.\n", expected_count, length(files)))
  results <- list()
  for (f in files) {
    e <- new.env(parent = emptyenv()); load(f, envir = e)
    results[[length(results) + 1L]] <- list(
      simulation_name = e$simulation_name, setting = e$setting,
      idx = e$idx, def = e$def, out = e$out
    )
  }
  list(sim_dir = sim_dir, results = results)
}

bind_rows_flexible <- function(dfs) {
  dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, dfs)
  if (length(dfs) == 0) return(data.frame())
  cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  aligned <- lapply(dfs, function(df) {
    miss <- setdiff(cols, names(df))
    if (length(miss) > 0) for (m in miss) df[[m]] <- NA
    df[cols]
  })
  do.call(rbind, aligned)
}

collect_main_df <- function(results) bind_rows_flexible(lapply(results, function(r) r$out$main_results))
collect_param_df <- function(results) bind_rows_flexible(lapply(results, function(r) r$out$parameter_results))

## Table 1 style: RMSE median (IQR)
make_table1_csv <- function(main_df, file_path) {
  fmt <- function(x) sprintf("%.3f (%.3f)", stats::median(x, na.rm = TRUE), stats::IQR(x, na.rm = TRUE))
  tbl <- stats::aggregate(
    cbind(measurement_rms, cox_rms) ~ setting + method + fit_model,
    data = main_df, FUN = fmt
  )
  names(tbl)[names(tbl) == "measurement_rms"] <- "Measurement_RMSE_Median_IQR"
  names(tbl)[names(tbl) == "cox_rms"] <- "Cox_RMSE_Median_IQR"
  write.csv(tbl, file = file_path, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", file_path))
  invisible(tbl)
}

## Table S2 style: BIAS, RMS, CP for loading, psi, beta per setting
make_tableS2_csv <- function(param_df, file_path) {
  df <- param_df[grepl("^(loading|psi|beta)", param_df$parameter), ]
  if (nrow(df) == 0) { cat("  No parameters for Table S2.\n"); return(invisible(NULL)) }
  key <- interaction(df$setting, df$method, df$fit_model, df$parameter, drop = TRUE)
  parts <- split(df, key)
  out <- lapply(parts, function(d) {
    data.frame(
      setting = d$setting[1], method = d$method[1], fit_model = d$fit_model[1],
      parameter = d$parameter[1], truth = d$truth[1],
      BIAS = mean(d$bias, na.rm = TRUE),
      RMS = sqrt(mean(d$sq_error, na.rm = TRUE)),
      CP = mean(d$covered, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  tbl <- do.call(rbind, out); rownames(tbl) <- NULL
  write.csv(tbl, file = file_path, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", file_path))
  invisible(tbl)
}

## Parameter summary (generic, with filter)
make_param_summary_csv <- function(param_df, file_path, param_filter = "^(loading|psi|beta)") {
  df <- param_df[grepl(param_filter, param_df$parameter), ]
  if (nrow(df) == 0) { cat("  No matching parameters.\n"); return(invisible(NULL)) }
  key <- interaction(df$setting, df$method, df$fit_model, df$parameter, drop = TRUE)
  parts <- split(df, key)
  out <- lapply(parts, function(d) {
    data.frame(
      setting = d$setting[1], method = d$method[1], fit_model = d$fit_model[1],
      parameter = d$parameter[1], truth = d$truth[1],
      BIAS = mean(d$bias, na.rm = TRUE),
      RMS = sqrt(mean(d$sq_error, na.rm = TRUE)),
      CP = mean(d$covered, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  tbl <- do.call(rbind, out); rownames(tbl) <- NULL
  write.csv(tbl, file = file_path, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", file_path))
  invisible(tbl)
}

## Model criteria (DIC, WAIC)
make_criteria_csv <- function(main_df, file_path) {
  out <- stats::aggregate(
    cbind(DIC, WAIC) ~ setting + method + fit_model,
    data = main_df, FUN = function(x) mean(x, na.rm = TRUE)
  )
  write.csv(out, file = file_path, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", file_path))
  invisible(out)
}

## Diagnostic: traceplot + R-hat
make_diagnostic_output <- function(results, out_dir) {
  r <- results[[1]]
  chains <- r$out$chain_objects
  if (is.null(chains) || length(chains) == 0) {
    cat("  No multi-chain objects for diagnostics.\n"); return(invisible(NULL))
  }
  chain_list <- chains[[1]]
  rhat <- compute_rhat(chain_list)
  # Output loading, psi, and beta trajectory summaries
  keep <- grepl("^(loading|psi|beta)", names(rhat))
  rhat_filtered <- rhat[keep]
  rhat_df <- data.frame(parameter = names(rhat_filtered), Rhat = as.numeric(rhat_filtered), stringsAsFactors = FALSE)
  write.csv(rhat_df, file = file.path(out_dir, "rhat.csv"), row.names = FALSE)
  cat(sprintf("  Saved: %s\n", file.path(out_dir, "rhat.csv")))
  trace_file <- file.path(out_dir, "traceplots.pdf")
  save_trace_plots(chain_list, trace_file)
  cat(sprintf("  Saved: %s\n", trace_file))
}

## ============================================================
## Analysis runners
## ============================================================

analysis_main <- function() {
  repo_root <- resolve_repo_root()
  source(file.path(repo_root, "cpp", "R", "helpers.R"))
  source(file.path(repo_root, "cpp", "R", "diagnostics_metrics.R"))
  loaded <- load_simulation_results(repo_root, "Simulation_Main", expected_count = 8)
  out_dir <- loaded$sim_dir
  main_df <- collect_main_df(loaded$results)
  param_df <- collect_param_df(loaded$results)
  cat("[Analysis_Main] Table 1: RMSE Median (IQR)\n")
  make_table1_csv(main_df, file.path(out_dir, "table1_rmse.csv"))
  cat("[Analysis_Main] Table S2: BIAS, RMS, CP for loading, psi, beta\n")
  make_tableS2_csv(param_df, file.path(out_dir, "tableS2_param_summary.csv"))
}

analysis_main_gp_kernels <- function() {
  repo_root <- resolve_repo_root()
  source(file.path(repo_root, "cpp", "R", "helpers.R"))
  source(file.path(repo_root, "cpp", "R", "diagnostics_metrics.R"))
  loaded <- load_simulation_results(repo_root, "Simulation_MainGPKernels", expected_count = 8)
  out_dir <- loaded$sim_dir
  main_df <- collect_main_df(loaded$results)
  param_df <- collect_param_df(loaded$results)
  cat("[Analysis_MainGPKernels] Table 1: RMSE Median (IQR)\n")
  make_table1_csv(main_df, file.path(out_dir, "table1_rmse.csv"))
  cat("[Analysis_MainGPKernels] Table S2: BIAS, RMS, CP for loading, psi, beta\n")
  make_tableS2_csv(param_df, file.path(out_dir, "tableS2_param_summary.csv"))
}

analysis_diagnostic <- function() {
  repo_root <- resolve_repo_root()
  source(file.path(repo_root, "cpp", "R", "helpers.R"))
  source(file.path(repo_root, "cpp", "R", "diagnostics_metrics.R"))
  loaded <- load_simulation_results(repo_root, "Simulation_Diagnostic", expected_count = 1)
  out_dir <- loaded$sim_dir
  cat("[Analysis_Diagnostic] traceplots + R-hat\n")
  make_diagnostic_output(loaded$results, out_dir)
}

analysis_error_sensitivity <- function() {
  repo_root <- resolve_repo_root()
  source(file.path(repo_root, "cpp", "R", "helpers.R"))
  source(file.path(repo_root, "cpp", "R", "diagnostics_metrics.R"))
  loaded <- load_simulation_results(repo_root, "Simulation_ErrorSensitivity", expected_count = 2)
  out_dir <- loaded$sim_dir
  param_df <- collect_param_df(loaded$results)
  cat("[Analysis_ErrorSensitivity] Table S4: BIAS, RMS, CP\n")
  make_param_summary_csv(param_df, file.path(out_dir, "tableS4_param_summary.csv"), "^(loading|psi|beta)")
}

analysis_imputation <- function() {
  repo_root <- resolve_repo_root()
  source(file.path(repo_root, "cpp", "R", "helpers.R"))
  source(file.path(repo_root, "cpp", "R", "diagnostics_metrics.R"))
  loaded <- load_simulation_results(repo_root, "Simulation_Imputation", expected_count = 4)
  out_dir <- loaded$sim_dir
  main_df <- collect_main_df(loaded$results)
  cat("[Analysis_Imputation] Table 1 style: RMSE Median (IQR)\n")
  make_table1_csv(main_df, file.path(out_dir, "table1_rmse.csv"))
}

analysis_model_comparison <- function() {
  repo_root <- resolve_repo_root()
  source(file.path(repo_root, "cpp", "R", "helpers.R"))
  source(file.path(repo_root, "cpp", "R", "diagnostics_metrics.R"))
  loaded <- load_simulation_results(repo_root, "Simulation_ModelComparison", expected_count = 4)
  out_dir <- loaded$sim_dir
  main_df <- collect_main_df(loaded$results)
  param_df <- collect_param_df(loaded$results)
  cat("[Analysis_ModelComparison] Table S3: BIAS, RMS, CP + DIC, WAIC\n")
  make_param_summary_csv(param_df, file.path(out_dir, "tableS3_param_summary.csv"), "^(loading|psi|beta)")
  make_criteria_csv(main_df, file.path(out_dir, "tableS3_criteria.csv"))
}

analysis_prior_sensitivity <- function() {
  repo_root <- resolve_repo_root()
  source(file.path(repo_root, "cpp", "R", "helpers.R"))
  source(file.path(repo_root, "cpp", "R", "diagnostics_metrics.R"))
  loaded <- load_simulation_results(repo_root, "Simulation_PriorSensitivity", expected_count = 1)
  out_dir <- loaded$sim_dir
  param_df <- collect_param_df(loaded$results)
  cat("[Analysis_PriorSensitivity] Table S5: BIAS, RMS, CP (Prior II only)\n")
  make_param_summary_csv(param_df, file.path(out_dir, "tableS5_param_summary.csv"), "^(loading|psi|beta)")
}

analysis_subinterval <- function() {
  repo_root <- resolve_repo_root()
  source(file.path(repo_root, "cpp", "R", "helpers.R"))
  source(file.path(repo_root, "cpp", "R", "diagnostics_metrics.R"))
  loaded <- load_simulation_results(repo_root, "Simulation_Subinterval", expected_count = 6)
  out_dir <- loaded$sim_dir
  main_df <- collect_main_df(loaded$results)
  cat("[Analysis_Subinterval] Table 1 style: RMSE Median (IQR)\n")
  make_table1_csv(main_df, file.path(out_dir, "table1_rmse.csv"))
}

analysis_varying_coef_alpha <- function() {
  repo_root <- resolve_repo_root()
  source(file.path(repo_root, "cpp", "R", "helpers.R"))
  source(file.path(repo_root, "cpp", "R", "diagnostics_metrics.R"))
  loaded <- load_simulation_results(repo_root, "Simulation_VaryingCoefAlpha", expected_count = 4)
  out_dir <- loaded$sim_dir
  main_df <- collect_main_df(loaded$results)
  cat("[Analysis_VaryingCoefAlpha] Table 1 style: RMSE Median (IQR)\n")
  make_table1_csv(main_df, file.path(out_dir, "table1_rmse.csv"))
}

analysis_missing <- function() {
  repo_root <- resolve_repo_root()
  source(file.path(repo_root, "cpp", "R", "helpers.R"))
  source(file.path(repo_root, "cpp", "R", "diagnostics_metrics.R"))
  loaded <- load_simulation_results(repo_root, "Simulation_Missing", expected_count = 6)
  out_dir <- loaded$sim_dir
  main_df <- collect_main_df(loaded$results)
  cat("[Analysis_Missing] Table 1 style: RMSE Median (IQR)\n")
  make_table1_csv(main_df, file.path(out_dir, "table1_rmse.csv"))
}
