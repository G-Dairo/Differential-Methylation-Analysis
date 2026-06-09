#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=enrichment_viz
#SBATCH --output=logs/11_enrichment_viz_%j.out
#SBATCH --error=logs/11_enrichment_viz_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# =============================================================
# Script: run_11_enrichment_visualization.sh
# Description: SLURM wrapper for 11_enrichment_visualization.R
# Usage: sbatch slurm/run_11_enrichment_visualization.sh
# =============================================================

module purge
source ~/.bashrc
conda activate methylation_env
export PATH="/home/gsd67/miniconda3/envs/methylation_env/bin:$PATH"
export PROJECT_DIR=/home/gsd67/Projects/Epigenetics/methylation

Rscript 11_enrichment_visualization.R
