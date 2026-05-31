## Latent Gaussian Process Joint Model for Integrative Analysis of Multi-Type Biomarkers and Initiation of Medication of Parkinson’s Disease

###### Junxuan Chen #, Zijian Ye #, Xiangnan Feng, and Kai Kang*

## Overview

This repository contains all code needed to reproduce the simulation results in the paper and supplement.

## Project Structure

```
LGPJM/

├── README.md

│

├── cpp/ # All source code

│ ├── R/ # Core engine (sourced by all scripts)

│ │ ├── helpers.R # Utility functions

│ │ ├── settings.R # Simulation parameter definitions

│ │ ├── data_generation.R # Data generation (ordinal + continuous + survival)

│ │ ├── imputation.R # LOCF imputation for missing data

│ │ ├── model_mcmc.R # MCMC sampler (C++ accelerated with pure R fallback)

│ │ ├── diagnostics_metrics.R # R-hat, WAIC, DIC, traceplots

│ │ └── run_pipeline.R # Pipeline: data gen → MCMC → parameter extraction

│ │

│ ├── cpp/ # Optional C++ acceleration (RcppArmadillo)

│ │ └── mcmc_core.cpp # C++ implementations of MCMC inner loops

│ │

│ ├── sim_common.R # Shared bootstrap for Simulation_*.R scripts

│ ├── analysis_common.R # Shared functions for Analysis_*.R scripts

│ │

│ ├── Simulation_Main.R # Table 1

│ ├── Simulation_MainGPKernels.R # Table S6, Table S7

│ ├── Simulation_Diagnostic.R # Figure S2

│ ├── Simulation_ErrorSensitivity.R # Table S4

│ ├── Simulation_Imputation.R # Table S9

│ ├── Simulation_ModelComparison.R # Table S3

│ ├── Simulation_PriorSensitivity.R # Table S5

│ ├── Simulation_Subinterval.R # Table S8

│ ├── Simulation_VaryingCoefAlpha.R # Table 1

│ ├── Simulation_Missing.R # Table S10

│ │

│ ├── Analysis_Main.R 

│ ├── Analysis_MainGPKernels.R

│ ├── Analysis_Diagnostic.R

│ ├── Analysis_ErrorSensitivity.R 

│ ├── Analysis_Imputation.R 

│ ├── Analysis_ModelComparison.R 

│ ├── Analysis_PriorSensitivity.R 

│ ├── Analysis_Subinterval.R 

│ ├── Analysis_VaryingCoefAlpha.R 

│ ├── Analysis_Missing.R 

│ │

│ └── run_all.sh # Batch runner for all simulations + analyses

```

## Quick Start

### Run one simulation setting

```bash
cd code

Rscript Simulation_Main.R 1 # run setting 1

Rscript Simulation_Main.R 2 # run setting 2

# ... up to setting 8
```

### Run all settings for a simulation

```bash
cd code

for i in $(seq 1 8); do Rscript Simulation_Main.R $i; done
```

### Generate analysis results (after all settings are done)

```bash
cd code

Rscript Analysis_Main.R # reads all RData, produces table1_rmse.csv
```

### Run everything

```bash
cd code

bash run_all.sh # all simulations + analyses

bash run_all.sh Simulation_Main # just Simulation_Main + its analysis
```

### Environment variable overrides

```bash
DS=100 IR=5000 BI=3000 N_CORES=4 Rscript Simulation_Main.R 1

# DS = number of replications (default: 500)

# IR = MCMC iterations (default: 10000)

# BI = burn-in (default: 6000)

# N_CORES = parallel cores (default: 500, auto-capped to available)

# SEED = random seed (default: 31)
```

## Simulation Details

### Simulation_Main / Simulation_MainGPKernels (Table 1, Table S6, Table S7)

- 8 settings: baseline ∈ {constant, linear, nonlinear} × N ∈ {200, 500} × CR ∈ {30%, 50%}
  
- Output: `table1_rmse.csv` — RMSE Median (IQR) for measurement model and Cox model
  

### Simulation_Diagnostic (Figure S2)

- 1 setting: nonlinear baseline, N=500, 3 chains
  
- Output: `rhat.csv` (R-hat per parameter), `traceplots.pdf`
  

### Simulation_ErrorSensitivity (Table S4)

- 2 settings: error ∈ {scaled t(3), shifted gamma}
  
- Output: `tableS4_param_summary.csv` — BIAS, RMS, CP for loading, ψ, β
  

### Simulation_Imputation (Table S9)

- 2 settings: N ∈ {200, 500}, internal vs LOCF imputation
  
- Output: `table1_rmse.csv` — RMSE Median (IQR)
  

### Simulation_ModelComparison (Table S3)

- 2 settings: case (a) linear ω truth, case (b) nonlinear ω truth
  
- Output: `tableS3_param_summary.csv` — BIAS, RMS, CP for loading, ψ, β
  
- Output: `tableS3_criteria.csv` — DIC, WAIC
  

### Simulation_PriorSensitivity (Table S5)

- 2 settings: Prior(I) and Prior(II)
  
- Output: `tableS5_param_summary.csv` — BIAS, RMS, CP for loading, ψ, β
  

### Simulation_Subinterval (Table S8)

- 6 settings: G ∈ {5, 10, 15} × N ∈ {200, 500}
  
- Output: `table1_rmse.csv` — RMSE Median (IQR)
  

### Simulation_VaryingCoefAlpha (Table 1)

- 4 settings: α mode ∈ {linear_tv, nonlinear_tv} × N ∈ {200, 500}
  
- Output: `table1_rmse.csv` — RMSE Median (IQR)
  

### Simulation_Missing (Table S10)

- 6 settings: missing rate ∈ {10%, 30%, 50%} × N ∈ {200, 500}
  
- Output: `table1_rmse.csv` — RMSE Median (IQR)
  

## Model Specification

The LGPJM consists of:

1. **Measurement model**: p markers (3 ordinal + 6 continuous), 2 latent factors via factor loading matrix Λ
  
2. **Latent process**: ω_k(t) modeled by Gaussian process with SE kernel
  
3. **Survival model**: piecewise constant baseline hazard with Cox-type regression
  

- h(t|u,ω) = λ_g · exp(β'u + α'ω(t))

## Requirements

- **Operating System**: Linux (tested on Ubuntu 24.04 LTS)
  
- **R**: 4.3.3
  
- **C++ Compiler**: g++ 13.3.0
  
- R packages: `MASS`, `zoo`, `Rcpp`, `RcppArmadillo`
  
- **Rscript** must be available on `PATH` (verify with `which Rscript`)
  

> **Note on reproducibility across platforms**: Main results were generated on Linux. Running on Windows/macOS may yield small numerical differences due to floating-point, parallel scheduling, and RNG backend differences.
