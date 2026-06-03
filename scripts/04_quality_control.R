# =============================================================
# Script: 04_quality_control.R
# Description: Quality control assessment of raw methylation data
#              Detects and removes failed probes and samples
# Input:  data/processed/rgSet.rds + data/raw/sample_sheet.csv
# Output: results/qc/, plots/qc_summary.pdf, data/processed/rgSet_qc.rds
# Usage:  Rscript scripts/04_quality_control.R
# =============================================================

library(minfi)

# Set working directory
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
setwd(project_dir)

# Load data from previous step
cat("Loading rgSet...\n")
rgSet <- readRDS("data/processed/rgSet.rds")
sample_info <- read.csv("data/raw/sample_sheet.csv", stringsAsFactors = FALSE)
cat("Loaded:", ncol(rgSet), "samples,", nrow(rgSet), "probes\n")

# Calculate detection p-values for all probes and samples
# These p-values indicate whether each probe's signal is above background
cat("Calculating detection p-values...\n")
detP <- detectionP(rgSet)

# Set quality control thresholds
detP_threshold <- 0.01          # Probes with p > 0.01 are considered "failed"
sample_failure_rate <- 0.05     # Remove samples with >5% failed probes
probe_failure_rate <- 0.10      # Remove probes that fail in >10% of samples

# Identify failed probes and samples using these criteria
failed_probes <- rowMeans(detP > detP_threshold, na.rm = TRUE) > probe_failure_rate
failed_samples <- colMeans(detP > detP_threshold, na.rm = TRUE) > sample_failure_rate

cat("Failed probes:", sum(failed_probes), "out of", nrow(detP), "\n")
cat("Failed samples:", sum(failed_samples), "out of", ncol(detP), "\n")

# Create quality control plots to visualize data quality
cat("Generating QC plots...\n")
pdf("plots/qc_summary.pdf", width = 12, height = 8)
par(mfrow = c(2, 2))

# Plot sample failure rates
sample_failure_rates <- colMeans(detP > detP_threshold, na.rm = TRUE)
barplot(sample_failure_rates, main = "Sample Quality Assessment",
        ylab = "Proportion of Failed Probes", las = 2, cex.names = 0.7,
        col = ifelse(sample_failure_rates > sample_failure_rate, "red", "lightblue"))
abline(h = sample_failure_rate, col = "red", lwd = 2, lty = 2)

# Plot probe failure distribution
probe_failure_rates <- rowMeans(detP > detP_threshold, na.rm = TRUE)
hist(probe_failure_rates, main = "Probe Quality Distribution",
     xlab = "Proportion of Failed Samples", breaks = 50, col = "lightblue")
abline(v = probe_failure_rate, col = "red", lwd = 2, lty = 2)

dev.off()
cat("QC plots saved to plots/qc_summary.pdf\n")

# Generate comprehensive quality control report
cat("Generating QC report...\n")
qcReport(rgSet, sampNames = sample_info$Sample_Name,
         sampGroups = sample_info$Group,
         pdf = "results/qc/qc_report.pdf")

# Remove poor quality samples if any are detected
if(any(failed_samples)) {
    cat("Removing", sum(failed_samples), "poor quality samples...\n")
    cat("Failed samples:", colnames(rgSet)[failed_samples], "\n")
    rgSet <- rgSet[, !failed_samples]
    sample_info <- sample_info[!failed_samples, ]
} else {
    cat("All samples passed quality control\n")
}

# Save QC results for reference
cat("Saving QC results...\n")
saveRDS(rgSet, "data/processed/rgSet_qc.rds")
saveRDS(detP, "data/processed/detP.rds")
write.csv(sample_info, "data/processed/sample_info_qc.csv", row.names = FALSE)

# Save QC summary table
qc_summary <- data.frame(
    Sample = colnames(detP),
    Failure_Rate = sample_failure_rates,
    Passed = !failed_samples
)
write.csv(qc_summary, "results/qc/qc_summary.csv", row.names = FALSE)
cat("QC summary saved to results/qc/qc_summary.csv\n")

cat("\nQC complete!\n")
cat("Samples remaining:", ncol(rgSet), "\n")
cat("Probes remaining:", nrow(rgSet), "\n")
