# =============================================================
# Script: 07_outlier_detection.R
# Description: Sample outlier detection using PCA and Mahalanobis
#              distance before differential methylation analysis.
#              Identifies samples that deviate from the expected
#              methylation patterns of their group.
# Input:  data/processed/filtered_data.rds
# Output: plots/sample_qc_pca.pdf
#         results/qc/outlier_summary.csv
#         data/processed/final_data.rds
# Usage:  Rscript scripts/07_outlier_detection.R
# =============================================================

library(minfi)
library(pheatmap)
library(RColorBrewer)

# Set working directory
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
setwd(project_dir)

# Load filtered data from previous step
cat("Loading filtered data...\n")
filtered_data <- readRDS("data/processed/filtered_data.rds")
m_complete    <- filtered_data$m_complete
beta_complete <- filtered_data$beta_complete
sample_info   <- filtered_data$sample_info
annotation    <- filtered_data$annotation

cat("Samples:", ncol(m_complete), "\n")
cat("Probes:", nrow(m_complete), "\n")

# ── PCA ───────────────────────────────────────────────────────────────────────
# PCA reduces 719,922 probe dimensions down to a handful of components
# that capture the major patterns of variation across samples.
# Samples that cluster far from their group in PC space are likely outliers.
cat("\nRunning PCA...\n")
pca_result  <- prcomp(t(m_complete), center = TRUE, scale. = FALSE)
var_explained <- summary(pca_result)$importance[2, 1:10]

cat("Variance explained by top 5 PCs:\n")
for(i in 1:5) {
    cat("  PC", i, ":", round(var_explained[i] * 100, 1), "%\n")
}

# ── Plots ─────────────────────────────────────────────────────────────────────
cat("\nGenerating QC plots...\n")
pdf("plots/sample_qc_pca.pdf", width = 12, height = 8)
par(mfrow = c(2, 2))

# Color palette for sample types
sample_colors <- as.integer(as.factor(sample_info$Sample_Type))
color_palette <- brewer.pal(length(unique(sample_info$Sample_Type)), "Set1")

# Plot 1: PC1 vs PC2 colored by sample type
plot(pca_result$x[, 1], pca_result$x[, 2],
     main = "PCA: Sample Clustering by Methylation Patterns",
     xlab = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)"),
     col  = color_palette[sample_colors],
     pch  = 16, cex = 1.2)
text(pca_result$x[, 1], pca_result$x[, 2],
     labels = sample_info$Sample_Type, cex = 0.5, pos = 3)
legend("topright",
       legend = unique(sample_info$Sample_Type),
       col    = color_palette[1:length(unique(sample_info$Sample_Type))],
       pch = 16, cex = 0.8)

# Plot 2: PC1 vs PC3
plot(pca_result$x[, 1], pca_result$x[, 3],
     main = "PCA: PC1 vs PC3",
     xlab = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
     ylab = paste0("PC3 (", round(var_explained[3] * 100, 1), "%)"),
     col  = color_palette[sample_colors],
     pch  = 16, cex = 1.2)
legend("topright",
       legend = unique(sample_info$Sample_Type),
       col    = color_palette[1:length(unique(sample_info$Sample_Type))],
       pch = 16, cex = 0.8)

# Plot 3: Scree plot — how much variance each PC explains
barplot(var_explained * 100,
        main = "Scree Plot: Variance Explained by Each PC",
        xlab = "Principal Component",
        ylab = "Variance Explained (%)",
        names.arg = paste0("PC", 1:10),
        col = "steelblue")

dev.off()
cat("PCA plots saved to plots/sample_qc_pca.pdf\n")

# ── Correlation heatmap ───────────────────────────────────────────────────────
cat("Generating correlation heatmap...\n")
sample_cors <- cor(beta_complete, use = "complete.obs")

annotation_col <- data.frame(
    Type      = sample_info$Sample_Type,
    Treatment = sample_info$Treatment,
    row.names = colnames(sample_cors)
)

pdf("plots/sample_correlation_heatmap.pdf", width = 12, height = 10)
pheatmap(sample_cors,
         annotation_col  = annotation_col,
         show_rownames   = FALSE,
         show_colnames   = FALSE,
         main            = "Sample-to-Sample Correlations",
         color           = colorRampPalette(brewer.pal(9, "Blues"))(100),
         breaks          = seq(0.95, 1, length.out = 101))
dev.off()
cat("Correlation heatmap saved to plots/sample_correlation_heatmap.pdf\n")

# ── Mahalanobis distance outlier detection ────────────────────────────────────
# Mahalanobis distance measures how far each sample is from the group
# centre in PC space, accounting for correlations between PCs.
# Samples beyond the 97.5th percentile of a chi-squared distribution
# are flagged as statistical outliers.
cat("\nRunning Mahalanobis distance outlier detection...\n")
pc_scores        <- pca_result$x[, 1:5]
mahal_dist       <- mahalanobis(pc_scores, colMeans(pc_scores), cov(pc_scores))
outlier_threshold <- qchisq(0.975, df = 5)

cat("Outlier threshold (chi-sq 97.5%, df=5):", round(outlier_threshold, 2), "\n")

outliers <- names(mahal_dist)[mahal_dist > outlier_threshold]

if(length(outliers) > 0) {
    cat("Outlier samples detected:", paste(outliers, collapse = ", "), "\n")
    cat("Consider investigating these samples for technical issues\n")
} else {
    cat("No outlier samples detected - data quality looks good\n")
}

# ── Save outlier summary ──────────────────────────────────────────────────────
outlier_summary <- data.frame(
    Sample           = rownames(pc_scores),
    Mahal_Distance   = mahal_dist,
    Is_Outlier       = mahal_dist > outlier_threshold,
    Sample_Type      = sample_info$Sample_Type,
    Treatment        = sample_info$Treatment
)
outlier_summary <- outlier_summary[order(outlier_summary$Mahal_Distance,
                                         decreasing = TRUE), ]
write.csv(outlier_summary, "results/qc/outlier_summary.csv", row.names = FALSE)
cat("Outlier summary saved to results/qc/outlier_summary.csv\n")

# ── Save final clean data ─────────────────────────────────────────────────────
# Remove confirmed outliers if any, then save as final_data.rds
if(length(outliers) > 0) {
    keep_samples  <- !colnames(m_complete) %in% outliers
    m_complete    <- m_complete[, keep_samples]
    beta_complete <- beta_complete[, keep_samples]
    sample_info   <- sample_info[keep_samples, ]
    cat("Removed", length(outliers), "outlier samples\n")
}

saveRDS(list(
    m_complete    = m_complete,
    beta_complete = beta_complete,
    sample_info   = sample_info,
    annotation    = annotation
), "data/processed/final_data.rds")

cat("\nFinal dataset saved to data/processed/final_data.rds\n")
cat("Final samples:", ncol(m_complete), "\n")
cat("Final probes:", nrow(m_complete), "\n")
cat("\nOutlier detection complete!\n")
