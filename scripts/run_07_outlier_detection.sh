#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=outlier_detection
#SBATCH --output=logs/07_outlier_%j.out
#SBATCH --error=logs/07_outlier_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

# =============================================================
# Script: run_07_outlier_detection.sh
# Description: SLURM wrapper for 07_outlier_detection.R
# Usage: sbatch slurm/run_07_outlier_detection.sh
# =============================================================

module purge
source ~/.bashrc
conda activate methylation_env
export PATH="/home/gsd67/miniconda3/envs/methylation_env/bin:$PATH"
export PROJECT_DIR=/home/gsd67/Projects/Epigenetics/methylation

Rscript 07_outlier_detection.R
