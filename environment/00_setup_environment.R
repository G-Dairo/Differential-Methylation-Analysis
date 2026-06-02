# =============================================================
# Script: 00_setup_environment.R
# Description: Install and load required Bioconductor packages
#              for EPIC v2.0 DNA methylation analysis from Lei Guo
#	       NGS 101.com tutorial	
# Note: Run environment/conda_setup.sh FIRST before this script
#       to install core Bioconductor dependencies via conda
# Author: Gbenga Dairo
# Date: 2026-06-02
# Usage: Rscript environment/00_setup_environment.R
# =============================================================

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = "3.21", ask = FALSE)

# NOTE: The following packages are already installed via conda_setup.sh
# and do NOT need to be installed here:
#   - GenomicRanges, SummarizedExperiment, AnnotationDbi, Biostrings
#   - XVector, Rsamtools, GenomicAlignments, GenomicFeatures
#   - rtracklayer, DelayedArray, HDF5Array, DelayedMatrixStats
#   - bumphunter, genefilter, annotate, minfi, GEOquery, DMRcate

# Install remaining higher-level Bioconductor packages
remaining_bioc <- c(
    "missMethyl",                                    # Specialized methylation functions
    "limma",                                         # Linear models for differential analysis
    "sva",                                           # Surrogate variable analysis
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38",  # EPIC v2.0 annotations
    "IlluminaHumanMethylationEPICv2manifest",         # EPIC v2.0 manifest
    "TxDb.Hsapiens.UCSC.hg38.knownGene",             # Gene annotations
    "org.Hs.eg.db"                                   # Gene symbol mappings
)
BiocManager::install(remaining_bioc, ask = FALSE, update = FALSE)

# Install CRAN packages
cran_packages <- c(
    "ggplot2",      # General plotting
    "pheatmap",     # Heatmap visualization
    "RColorBrewer", # Color palettes
    "dplyr",        # Data manipulation
    "reshape2"      # Data reshaping
)
install.packages(cran_packages)

# Load and verify all key libraries
cat("Loading and verifying packages...\n")

library(minfi)
library(limma)
library(DMRcate)
library(sva)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

cat("All packages loaded successfully!\n")

# Configure R session options for optimal performance
options(
    stringsAsFactors = FALSE,
    scipen = 999,
    max.print = 100,
    width = 120
)

set.seed(12345)  # Ensure reproducible results

# Set up working directory structure
# Uses current directory if running from project root via SLURM
# Override by setting PROJECT_DIR environment variable
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
setwd(project_dir)

# Create organized directory structure for analysis
directories <- c(
    "data/raw",               # Raw .idat files
    "data/processed",         # Processed data objects
    "results/qc",             # Quality control outputs
    "results/differential",   # Differential methylation results
    "results/annotation",     # Annotation results
    "plots",                  # All visualization outputs
    "reports"                 # Analysis reports
)
sapply(directories, function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))
cat("Directory structure ready!\n")
