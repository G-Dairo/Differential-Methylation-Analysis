#!/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=download_GSE240469
#SBATCH --output=logs/01_download_%j.out
#SBATCH --error=logs/01_download_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# =============================================================
# Script: 01_download_data.sh
# Description: Download raw IDAT files for GSE240469 from GEO
#              40 samples, EPIC v2.0 array, ~580MB
# Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240469
# Usage: sbatch scripts/01_download_data.sh
# =============================================================

PROJECT_DIR=/home/gsd67/Projects/Epigenetics/methylation

cd $PROJECT_DIR/data/raw

echo "Downloading GSE240469 raw IDAT files..."

# Download the raw tar file containing all IDAT files
wget -c \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE240469&format=file" \
    -O GSE240469_RAW.tar

echo "Download complete. Extracting IDAT files..."

# Extract the tar file
tar -xf GSE240469_RAW.tar

# Decompress any gzipped idat files
gunzip -f *.idat.gz 2>/dev/null || true

echo "Extraction complete. Cleaning up..."

# Remove tar file to save space (already extracted)
rm GSE240469_RAW.tar

echo "Done! IDAT files are in $PROJECT_DIR/data/raw"
ls *.idat | wc -l
echo "idat files found"
