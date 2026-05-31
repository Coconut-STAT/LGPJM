#!/bin/bash
## run_all.sh — Run all simulations then all analyses
## Usage:
##   bash run_all.sh                     # run everything
##   bash run_all.sh Simulation_Main     # run only Simulation_Main (all settings)
##
## Helper: run_range <script> <from> <to>
##   run_range Simulation_Main.R 1 8     # run settings 1 through 8

set -e

run_range() {
  local script=$1 from=$2 to=$3
  echo "=== Running $script settings $from..$to ==="
  for i in $(seq $from $to); do
    echo "--- $script $i ---"
    Rscript "$script" "$i"
  done
}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

TARGET="${1:-all}"

if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_Main" ]; then
  run_range Simulation_Main.R 1 8
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_MainGPKernels" ]; then
  run_range Simulation_MainGPKernels.R 1 8
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_Diagnostic" ]; then
  run_range Simulation_Diagnostic.R 1 1
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_ErrorSensitivity" ]; then
  run_range Simulation_ErrorSensitivity.R 1 2
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_Imputation" ]; then
  run_range Simulation_Imputation.R 1 4
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_ModelComparison" ]; then
  run_range Simulation_ModelComparison.R 1 4
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_PriorSensitivity" ]; then
  run_range Simulation_PriorSensitivity.R 1 1
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_Subinterval" ]; then
  run_range Simulation_Subinterval.R 1 6
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_VaryingCoefAlpha" ]; then
  run_range Simulation_VaryingCoefAlpha.R 1 4
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_Missing" ]; then
  run_range Simulation_Missing.R 1 6
fi

echo ""
echo "=== Running Analysis scripts ==="

if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_Main" ]; then
  Rscript Analysis_Main.R
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_MainGPKernels" ]; then
  Rscript Analysis_MainGPKernels.R
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_Diagnostic" ]; then
  Rscript Analysis_Diagnostic.R
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_ErrorSensitivity" ]; then
  Rscript Analysis_ErrorSensitivity.R
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_Imputation" ]; then
  Rscript Analysis_Imputation.R
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_ModelComparison" ]; then
  Rscript Analysis_ModelComparison.R
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_PriorSensitivity" ]; then
  Rscript Analysis_PriorSensitivity.R
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_Subinterval" ]; then
  Rscript Analysis_Subinterval.R
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_VaryingCoefAlpha" ]; then
  Rscript Analysis_VaryingCoefAlpha.R
fi
if [ "$TARGET" = "all" ] || [ "$TARGET" = "Simulation_Missing" ]; then
  Rscript Analysis_Missing.R
fi

echo ""
echo "=== All done ==="
