#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=methylation_qc
#SBATCH --output=logs/04_qc_%j.out
#SBATCH --error=logs/04_qc_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

# =============================================================
# Script: run_04_quality_control.sh
# Description: SLURM wrapper for 04_quality_control.R
# Usage: sbatch slurm/run_04_quality_control.sh
# =============================================================

module purge
source ~/.bashrc
conda activate methylation_env
export PATH="/home/gsd67/miniconda3/envs/methylation_env/bin:$PATH"
export PROJECT_DIR=/home/gsd67/Projects/Epigenetics/methylation

Rscript 04_quality_control.R
