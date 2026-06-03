# =============================================================
# Script: 03_load_data.R
# Description: Import raw IDAT files into R as RGChannelSet
#              Reads both Red and Green channels for all 40 samples
# Input:  data/raw/*.idat + data/raw/sample_sheet.csv
# Output: data/processed/rgSet.rds
# Usage:  Rscript scripts/03_load_data.R
# =============================================================

library(minfi)

# Set working directory
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
setwd(project_dir)

# Load sample sheet created in previous step
cat("Loading sample sheet...\n")
sample_info <- read.csv("data/raw/sample_sheet.csv", stringsAsFactors = FALSE)
cat("Samples to load:", nrow(sample_info), "\n")

# Import methylation data using the sample information
# This function reads both red and green .idat files for each sample
cat("Reading IDAT files - this may take several minutes...\n")
rgSet <- read.metharray.exp(targets = sample_info, recursive = TRUE, extended = TRUE, force = TRUE)

# Print summary of imported data
cat("Data import complete!\n")
print(rgSet)

# Get probe information to understand array content
cat("\nExtracting probe information...\n")
probe_info <- getProbeInfo(rgSet)
cat("Probe info dimensions:", nrow(probe_info), "probes x", ncol(probe_info), "columns\n")

# Save rgSet object for use in subsequent steps
cat("Saving rgSet to data/processed/rgSet.rds...\n")
saveRDS(rgSet, "data/processed/rgSet.rds")
saveRDS(probe_info, "data/processed/probe_info.rds")

cat("Done!\n")
