#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=differential_methylation
#SBATCH --output=logs/08_dmp_%j.out
#SBATCH --error=logs/08_dmp_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

# =============================================================
# Script: run_08_differential_methylation.sh
# Description: SLURM wrapper for 08_differential_methylation.R
# Usage: sbatch slurm/run_08_differential_methylation.sh
# =============================================================

module purge
source ~/.bashrc
conda activate methylation_env
export PATH="/home/gsd67/miniconda3/envs/methylation_env/bin:$PATH"
export PROJECT_DIR=/home/gsd67/Projects/Epigenetics/methylation

Rscript 08_differential_methylation.R
