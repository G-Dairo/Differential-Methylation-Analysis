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
 
# Stage 1: Core Bioconductor infrastructure (must come first)
stage1 <- c(
    "BiocGenerics",
    "S4Vectors",
    "IRanges",
    "GenomeInfoDb",
    "GenomicRanges",
    "SummarizedExperiment",
    "Biostrings",
    "XVector",
    "AnnotationDbi",
    "Biobase"
)
BiocManager::install(stage1, ask = FALSE, update = FALSE)

# Stage 2: Annotation and genomic feature packages
stage2 <- c(
    "GenomicFeatures",
    "rtracklayer",
    "GenomicAlignments",
    "Rsamtools",
    "BSgenome",
    "biomaRt",
    "AnnotationHub",
    "ExperimentHub"
)
BiocManager::install(stage2, ask = FALSE, update = FALSE)

# Stage 3: Illumina array infrastructure
stage3 <- c(
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)
BiocManager::install(stage3, ask = FALSE, update = FALSE)

# Stage 4: minfi and direct dependents
BiocManager::install("minfi", ask = FALSE, update = FALSE)

# Stage 5: Packages that depend on minfi
stage5 <- c(
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    "IlluminaHumanMethylationEPICv2manifest",
    "missMethyl",
    "DMRcate"
)
BiocManager::install(stage5, ask = FALSE, update = FALSE)

# Stage 6: Remaining analysis packages
stage6 <- c(
    "limma",
    "sva",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db",
    "bsseq",
    "Gviz"
)
BiocManager::install(stage6, ask = FALSE, update = FALSE)

# Stage 7: CRAN packages
cran_packages <- c(
    "ggplot2",
    "pheatmap",
    "RColorBrewer",
    "dplyr",
    "reshape2"
)
install.packages(cran_packages)

# Load and verify key libraries
library(minfi)
library(limma)
library(DMRcate)
library(sva)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

cat("All packages loaded successfully!\n")
 
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

