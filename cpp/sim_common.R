## sim_common.R — Shared bootstrap for all Simulation_*.R scripts
## Usage: Rscript Simulation_Main.R <setting_index> [--datasets 500 --seed 31 ...]
## Output: output/<SimName>/<SimName>_<index>.RData (one file per run)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

resolve_repo_root <- function() {
  argv_file <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(argv_file) == 0) argv_file <- "--file=sim_common.R"
  this_file <- sub("^--file=", "", argv_file[[1]])
  normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = FALSE)
}

source_core_files <- function(repo_root) {
  suppressPackageStartupMessages({
    library(MASS)
  })
  # Try to load C++ acceleration (optional)
  cpp_path <- file.path(repo_root, "cpp", "cpp", "mcmc_core.cpp")
  if (file.exists(cpp_path)) {
    tryCatch({
      suppressPackageStartupMessages({ library(Rcpp); library(RcppArmadillo) })
      Rcpp::sourceCpp(cpp_path, verbose = FALSE)
      cat("[C++ acceleration enabled]\n")
    }, error = function(e) cat("[C++ not available, using pure R]\n"))
  }
  src_dir <- file.path(repo_root, "cpp", "R")
  source(file.path(src_dir, "helpers.R"))
  source(file.path(src_dir, "settings.R"))
  source(file.path(src_dir, "data_generation.R"))
  source(file.path(src_dir, "imputation.R"))
  source(file.path(src_dir, "mcmc_helpers.R"))
  source(file.path(src_dir, "model_mcmc.R"))
  source(file.path(src_dir, "diagnostics_metrics.R"))
  source(file.path(src_dir, "run_pipeline.R"))
}

parse_setting_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  pos <- NULL
  if (length(args) > 0 && !startsWith(args[[1]], "--")) {
    pos <- suppressWarnings(as.integer(args[[1]]))
    args <- args[-1]
  }
  cli <- parse_cli_args(args)
  cli$setting_index <- pos
  cli
}

build_setting <- function(overrides = list()) {
  st <- modifyList(get_common_setting(), overrides)
  apply_setting_functions(st)
}

default_cli_controls <- function(cli, def) {
  env <- function(var, default) as.character(Sys.getenv(var, unset = default))
  if (is.null(cli$datasets))   cli$datasets   <- env("DS",      def$datasets %||% "500")
  if (is.null(cli$iterations)) cli$iterations <- env("IR",      "10000")
  if (is.null(cli$burnin))     cli$burnin     <- env("BI",      "6000")
  if (is.null(cli$chains))     cli$chains     <- env("CHAINS",  def$chains %||% "1")
  if (is.null(cli$n_cores))    cli$n_cores    <- env("N_CORES", "500")
  if (is.null(cli$seed))       cli$seed       <- env("SEED",    "31")
  if (is.null(cli$parallel))   cli$parallel   <- "true"
  cli
}

print_setting_catalog <- function(defs) {
  cat("Available settings:\n")
  for (i in seq_along(defs)) {
    cat(sprintf("  %d: %s\n", i, defs[[i]]$label %||% defs[[i]]$tag %||% paste("setting", i)))
  }
}

## Main entry point: run one setting, save one RData file
run_simulation_script <- function(simulation_name, defs) {
  repo_root <- resolve_repo_root()
  source_core_files(repo_root)
  cli <- parse_setting_args()

  if (as_flag(cli$help, FALSE) || is.null(cli$setting_index) || is.na(cli$setting_index)) {
    script_name <- basename(commandArgs()[grep("^--file=", commandArgs())][1] %||% simulation_name)
    cat(sprintf("Usage: Rscript %s <setting_index> [--datasets 500 --n_cores 500 --seed 31]\n", script_name))
    print_setting_catalog(defs)
    quit(save = "no", status = 0)
  }
  idx <- cli$setting_index
  if (idx < 1 || idx > length(defs)) stop(sprintf("Invalid setting index %s. Choose 1..%s", idx, length(defs)))

  def <- defs[[idx]]
  cli <- default_cli_controls(cli, def)

  setting <- build_setting(def$setting_overrides %||% list())
  setting$setting_id <- def$tag %||% sprintf("%s_setting_%d", simulation_name, idx)

  # CLI override hooks
  if (!is.null(cli$missing_rate))      setting$missing_rate <- as_num(cli$missing_rate, setting$missing_rate)
  if (!is.null(cli$block_missing_prob)) setting$block_missing_prob <- as_num(cli$block_missing_prob, setting$block_missing_prob)
  if (!is.null(cli$missing_mechanism)) setting$missing_mechanism <- toupper(as.character(cli$missing_mechanism))
  setting <- apply_setting_functions(setting)

  methods <- def$methods %||% "internal"
  fit_model <- def$fit_model %||% "gp"
  cli$methods <- methods
  cli$fit_model <- fit_model
  if (!is.null(def$chains))   cli$chains <- as.character(def$chains)
  if (!is.null(def$datasets)) cli$datasets <- as.character(def$datasets)
  if (!is.null(def$cli_overrides)) {
    for (nm in names(def$cli_overrides)) cli[[nm]] <- as.character(def$cli_overrides[[nm]])
  }

  cat(sprintf("[%s] running setting %d: %s\n", simulation_name, idx, def$label %||% def$tag))
  cat(sprintf("  N=%d, baseline=%s, CR=%.2f, alpha=%s, omega=%s, methods=%s\n",
    setting$n_subject, setting$baseline_mode, setting$censoring_rate,
    setting$alpha_mode, setting$omega_mode, methods))

  # Run pipeline
  out <- run_setting_pipeline(setting, cli)

  # Save single RData file: output/<SimName>/<SimName>_<index>.RData
  out_dir <- file.path(repo_root, "output", simulation_name)
  ensure_dir(out_dir)
  rdata_file <- file.path(out_dir, sprintf("%s_%d.RData", simulation_name, idx))
  save(simulation_name, setting, idx, def, out, file = rdata_file)
  cat(sprintf("[%s] done. saved: %s\n", simulation_name, rdata_file))
}
