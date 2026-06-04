#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=probe_filtering
#SBATCH --output=logs/06_probe_filtering_%j.out
#SBATCH --error=logs/06_probe_filtering_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

# =============================================================
# Script: run_06_probe_filtering.sh
# Description: SLURM wrapper for 06_probe_filtering.R
# Usage: sbatch slurm/run_06_probe_filtering.sh
# =============================================================

module purge
source ~/.bashrc
conda activate methylation_env
export PATH="/home/gsd67/miniconda3/envs/methylation_env/bin:$PATH"
export PROJECT_DIR=/home/gsd67/Projects/Epigenetics/methylation

Rscript 06_probe_filtering.R
