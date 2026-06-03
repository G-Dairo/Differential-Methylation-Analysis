#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=load_methylation_data
#SBATCH --output=logs/03_load_data_%j.out
#SBATCH --error=logs/03_load_data_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# =============================================================
# Script: run_03_load_data.sh
# Description: SLURM wrapper for 03_load_data.R
# Usage: sbatch slurm/run_03_load_data.sh
# =============================================================

module purge
source ~/.bashrc
conda activate methylation_env
export PATH="/home/gsd67/miniconda3/envs/methylation_env/bin:$PATH"
export PROJECT_DIR=/home/gsd67/Projects/Epigenetics/methylation

Rscript 03_load_data.R
