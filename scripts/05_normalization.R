# =============================================================
# Script: 05_normalization.R
# Description: Methylation quantification and normalization
#              Applies functional normalization to correct for
#              technical variation, then extracts beta and M values
# Input:  data/processed/rgSet_qc.rds
#         data/processed/detP.rds
#         data/processed/sample_info_qc.csv
# Output: data/processed/normalized_data.rds
# Usage:  Rscript scripts/05_normalization.R
# =============================================================

library(minfi)

# Set working directory
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
setwd(project_dir)

# Load QC-cleaned data from previous step
# IMPORTANT: We load rgSet_qc.rds NOT rgSet.rds
# rgSet_qc.rds has the failed sample (LNCAP_250_1) already removed
cat("Loading QC-cleaned data...\n")
rgSet <- readRDS("data/processed/rgSet_qc.rds")
detP <- readRDS("data/processed/detP.rds")
sample_info <- read.csv("data/processed/sample_info_qc.csv", stringsAsFactors = FALSE)

cat("Samples loaded:", ncol(rgSet), "\n")
cat("Probes loaded:", nrow(rgSet), "\n")

# Recalculate failed probes from detP
# Only keep detP columns matching remaining samples after QC
detP_threshold <- 0.01
probe_failure_rate <- 0.10
detP_qc <- detP[, colnames(rgSet)]  # align detP to QC-passing samples only
failed_probes <- rowMeans(detP_qc > detP_threshold, na.rm = TRUE) > probe_failure_rate
cat("Failed probes to remove:", sum(failed_probes), "\n")

# Apply functional normalization to correct for technical variation
# nPCs = 2: uses 2 principal components from control probes
# bgCorr = TRUE: corrects for background fluorescence
# dyeCorr = TRUE: corrects for Red/Green dye bias
cat("Applying functional normalization...\n")
mSet_norm <- preprocessFunnorm(rgSet, nPCs = 2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE)
cat("Normalization complete!\n")
print(mSet_norm)

# Remove failed probes after normalization
if(any(failed_probes)) {
    common_probes <- intersect(rownames(mSet_norm), names(failed_probes))
    failed_probes_subset <- failed_probes[common_probes]
    mSet_norm <- mSet_norm[!failed_probes_subset, ]
    cat("Removed", sum(failed_probes_subset), "failed probes\n")
} else {
    cat("All probes passed quality control\n")
}

cat("Probes remaining after filtering:", nrow(mSet_norm), "\n")

# Extract beta-values and M-values
# Beta values: proportion of methylation (0-1), intuitive for interpretation
# M values: log2 ratio of methylated/unmethylated, better for statistical testing
cat("Extracting beta and M values...\n")
beta_values <- getBeta(mSet_norm)
m_values <- getM(mSet_norm)

cat("Beta values dimensions:", nrow(beta_values), "probes x", ncol(beta_values), "samples\n")
cat("M values dimensions:", nrow(m_values), "probes x", ncol(m_values), "samples\n")

# Quick sanity check - beta values should be between 0 and 1
cat("\nBeta value range:", round(min(beta_values, na.rm = TRUE), 4),
    "to", round(max(beta_values, na.rm = TRUE), 4), "\n")

# Save all processed data for downstream analysis
cat("Saving normalized data...\n")
saveRDS(list(
    mSet_norm   = mSet_norm,
    beta_values = beta_values,
    m_values    = m_values,
    sample_info = sample_info
), "data/processed/normalized_data.rds")

cat("\nNormalization complete!\n")
cat("Output saved to data/processed/normalized_data.rds\n")
