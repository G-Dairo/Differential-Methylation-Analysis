#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=create_samplesheet
#SBATCH --output=logs/02_samplesheet_%j.out
#SBATCH --error=logs/02_samplesheet_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# =============================================================
# Script: run_02_create_samplesheet.sh
# Description: SLURM wrapper for 02_create_samplesheet.R
# Usage: sbatch slurm/run_02_create_samplesheet.sh
# =============================================================

module purge
source ~/.bashrc
conda activate methylation_env
export PATH="/home/gsd67/miniconda3/envs/methylation_env/bin:$PATH"
export PROJECT_DIR=/home/gsd67/Projects/Epigenetics/methylation

Rscript 02_create_samplesheet.R
