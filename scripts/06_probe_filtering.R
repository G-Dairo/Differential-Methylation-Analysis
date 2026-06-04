# =============================================================
# Script: 06_probe_filtering.R
# Description: Additional probe filtering to remove problematic
#              probes that could confound downstream analysis.
#              Removes cross-reactive, SNP-affected, and sex
#              chromosome probes, then ensures complete data.
# Input:  data/processed/normalized_data.rds
# Output: data/processed/filtered_data.rds
# Usage:  Rscript scripts/06_probe_filtering.R
# =============================================================

library(minfi)

# Set working directory
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
setwd(project_dir)

# Load normalized data from previous step
cat("Loading normalized data...\n")
normalized_data <- readRDS("data/processed/normalized_data.rds")
mSet_norm    <- normalized_data$mSet_norm
beta_values  <- normalized_data$beta_values
m_values     <- normalized_data$m_values
sample_info  <- normalized_data$sample_info

cat("Starting probes:", nrow(mSet_norm), "\n")
cat("Samples:", ncol(mSet_norm), "\n")

# Get comprehensive annotation data for all probes
cat("\nFetching probe annotations...\n")
annotation <- getAnnotation(mSet_norm)

# ── Filter 1: Cross-reactive and SNP-affected probes ──────────────────────────
# Cross-reactive probes can bind to multiple genomic locations giving
# misleading methylation readings.
# SNP-affected probes sit on or near genetic variants — a SNP under
# the probe changes the binding, making the signal reflect genotype
# rather than methylation status.
cross_reactive <- rownames(annotation)[!is.na(annotation$Probe_rs) |
                                       !is.na(annotation$CpG_rs)   |
                                       !is.na(annotation$SBE_rs)]
cat("Cross-reactive / SNP-affected probes:", length(cross_reactive), "\n")

# ── Filter 2: Sex chromosome probes ───────────────────────────────────────────
# Probes on chrX and chrY behave differently between males and females
# due to X-inactivation and Y-chromosome absence in females.
# Removing them prevents biological sex from confounding the analysis.
sex_chr_probes <- rownames(annotation)[annotation$chr %in% c("chrX", "chrY")]
cat("Sex chromosome probes:", length(sex_chr_probes), "\n")

# ── Combine and apply all filters ─────────────────────────────────────────────
problematic_probes <- unique(c(cross_reactive, sex_chr_probes))
keep_probes <- !rownames(mSet_norm) %in% problematic_probes
cat("Total problematic probes to remove:", sum(!keep_probes), "\n")

mSet_filtered <- mSet_norm[keep_probes, ]
beta_filtered <- getBeta(mSet_filtered)
m_filtered    <- getM(mSet_filtered)
cat("Probes after filtering:", nrow(mSet_filtered), "\n")

# ── Filter 3: Incomplete and non-finite probes ────────────────────────────────
# Some probes may have NA or infinite M-values after normalization.
# These arise from probes where methylated or unmethylated signal is
# zero — log2(0) = -Inf. They must be removed before statistical testing.
cat("\nRemoving probes with incomplete or non-finite values...\n")
complete_probes <- complete.cases(m_filtered) &
                   apply(m_filtered, 1, function(x) all(is.finite(x)))

m_complete    <- m_filtered[complete_probes, ]
beta_complete <- beta_filtered[complete_probes, ]
cat("Probes removed (incomplete/non-finite):", sum(!complete_probes), "\n")
cat("Final probes remaining:", nrow(m_complete), "\n")

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n── Filtering Summary ──────────────────────────────\n")
cat("Started with:               ", nrow(mSet_norm), "probes\n")
cat("Removed (cross-reactive/SNP):", length(cross_reactive), "probes\n")
cat("Removed (sex chromosomes):   ", length(sex_chr_probes), "probes\n")
cat("Removed (incomplete/nonfinite):", sum(!complete_probes), "probes\n")
cat("Final clean probes:          ", nrow(m_complete), "probes\n")
cat("Samples:                     ", ncol(m_complete), "\n")
cat("───────────────────────────────────────────────────\n")

# Save filtered data
cat("\nSaving filtered data...\n")
saveRDS(list(
    mSet_filtered = mSet_filtered,
    beta_complete = beta_complete,
    m_complete    = m_complete,
    sample_info   = sample_info,
    annotation    = annotation[rownames(m_complete), ]
), "data/processed/filtered_data.rds")

cat("Saved to data/processed/filtered_data.rds\n")
cat("\nProbe filtering complete!\n")
