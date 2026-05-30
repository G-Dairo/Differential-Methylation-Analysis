#!/bin/bash
#SBATCH --partition=main 
#SBATCH --job-name=methylation_setup
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --output=logs/00_setup_%j.out
#SBATCH --error=logs/00_setup_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

module purge
# Activate conda environment with R 4.5.3
source ~/.bashrc
conda activate methylation_env

# Explicitly use conda's R
export PATH="/home/gsd67/miniconda3/envs/methylation_env/bin:$PATH"

which R
R --version
which Rscript
Rscript --version
  

# Run R script
Rscript 00_setup_environment.R
