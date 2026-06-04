#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=methylation_normalization
#SBATCH --output=logs/05_normalization_%j.out
#SBATCH --error=logs/05_normalization_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

# =============================================================
# Script: run_05_normalization.sh
# Description: SLURM wrapper for 05_normalization.R
# Usage: sbatch slurm/run_05_normalization.sh
# =============================================================

module purge
source ~/.bashrc
conda activate methylation_env
export PATH="/home/gsd67/miniconda3/envs/methylation_env/bin:$PATH"
export PROJECT_DIR=/home/gsd67/Projects/Epigenetics/methylation

Rscript 05_normalization.R
