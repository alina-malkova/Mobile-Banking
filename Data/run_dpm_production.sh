#!/bin/bash
#SBATCH --job-name=dpm_production
#SBATCH --output=dpm_production_%j.out
#SBATCH --error=dpm_production_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amalkova@fit.edu

# ==============================================================================
# Bayesian DPM Production Run — SLURM Job Script
#
# Usage:
#   Option 1 (full sample):   sbatch run_dpm_production.sh
#   Option 2 (subsample):     sbatch run_dpm_production.sh 50000
#
# Requirements:
#   - Copy Data/ directory to cluster scratch space
#   - Set DATA_DIR below to the cluster path
#   - Ensure R 4.x with rstan, haven, dplyr is available (module load)
#
# Estimated runtime:
#   N=50K:  ~12-24 hours, ~32GB RAM
#   N=95K:  ~24-48 hours, ~64GB RAM
# ==============================================================================

# --- Adjust these for your cluster ---
# module load R/4.3.0        # Uncomment and adjust for your cluster
# module load gcc/12.2.0     # Stan needs a C++ compiler

# Set data directory (copy Data/ folder to cluster scratch)
export DATA_DIR="${SLURM_TMPDIR:-/scratch/${USER}/mobile_banking}/Data"

# If data not on scratch, use local path
if [ ! -d "$DATA_DIR" ]; then
  export DATA_DIR="/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/_Research/Mobile_Money_Banking/Mobile banking USA/Data"
fi

echo "============================================"
echo "DPM Production Run"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "DATA_DIR: $DATA_DIR"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-$(nproc)}"
echo "Memory: ${SLURM_MEM_PER_NODE:-unknown}MB"
echo "============================================"

# Run R script (pass subsample size as argument, if any)
Rscript "${DATA_DIR}/phase8_bayesian_dpm_production.R" $1

echo ""
echo "Job completed: $(date)"
