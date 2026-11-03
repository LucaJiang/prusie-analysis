#!/bin/bash
#SBATCH --job-name=prusie-paired
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
# Set partition/account as required by your installation.
set -euo pipefail
: "${ANALYSIS_CONFIG:?Set ANALYSIS_CONFIG to your populated YAML configuration}"
: "${ANALYSIS_PYTHON:?Set ANALYSIS_PYTHON to the pinned Python executable}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RAYON_NUM_THREADS=1 PRUSIE_NUM_THREADS=1
exec "$ANALYSIS_PYTHON" "$(dirname "$0")/campaign.py" --config "$ANALYSIS_CONFIG" --workers 1
