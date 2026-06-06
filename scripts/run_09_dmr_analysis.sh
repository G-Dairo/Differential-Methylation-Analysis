#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=dmr_analysis
#SBATCH --output=logs/09_dmr_%j.out
#SBATCH --error=logs/09_dmr_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

# =============================================================
# Script: run_09_dmr_analysis.sh
# Description: SLURM wrapper for 09_dmr_analysis.R
# Usage: sbatch slurm/run_09_dmr_analysis.sh
# =============================================================

module purge
source ~/.bashrc
conda activate methylation_env
export PATH="/home/gsd67/miniconda3/envs/methylation_env/bin:$PATH"
export PROJECT_DIR=/home/gsd67/Projects/Epigenetics/methylation

Rscript 09_dmr_analysis.R
