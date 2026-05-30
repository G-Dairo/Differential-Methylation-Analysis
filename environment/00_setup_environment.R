# =============================================================
# Script: 00_setup_environment.R
# Description: Install and load required Bioconductor packages
#              for EPIC v2.0 DNA methylation analysis
# Original tutorial can be found on ngs101.com website by Lei Guo
# Author: Gbenga Dairo
# Date: 05-29-2026
# Usage: Rscript environment/00_setup_environment.R
# =============================================================

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install essential Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 
BiocManager::install(version = "3.21", ask = FALSE)
 
# Core methylation analysis packages
essential_packages <- c(
    "minfi",              # Core package for Illumina methylation arrays
    "missMethyl",         # Specialized functions for methylation analysis
    "limma",              # Linear models for differential analysis
    "DMRcate",            # Differential methylated regions detection
    "sva",                # Surrogate variable analysis for batch effects
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38", # EPIC v2.0 annotations
    "IlluminaHumanMethylationEPICv2manifest",        # EPIC v2.0 manifest
    "TxDb.Hsapiens.UCSC.hg38.knownGene",             # Gene annotations
    "org.Hs.eg.db",                                  # Gene symbol mappings
    "ggplot2",            # General plotting
    "pheatmap",           # Heatmap visualization
    "RColorBrewer",       # Color palettes
    "dplyr",              # Data manipulation
    "reshape2"            # Data reshaping
)
 
BiocManager::install(essential_packages, update = TRUE, ask = FALSE)
 
# Load essential libraries
library(minfi)
library(limma)
library(DMRcate)
library(sva)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
 
# Set up working directory structure
project_dir <- "~/methylation"
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
setwd(project_dir)
 
# Create organized directory structure for analysis
directories <- c(
    "data/raw",           # Raw .idat files
    "data/processed",     # Processed data objects
    "results/qc",         # Quality control outputs
    "results/differential", # Differential methylation results
    "results/annotation", # Annotation results
    "plots",              # All visualization outputs
    "reports"             # Analysis reports
)
 
sapply(directories, function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))
 
# Configure R session options for optimal performance
options(
    stringsAsFactors = FALSE,
    scipen = 999,
    max.print = 100,
    width = 120
)
 
set.seed(12345)  # Ensure reproducible results

