# =============================================================
# Script: 09_dmr_analysis.R
# Description: Identify Differentially Methylated Regions (DMRs)
#              using DMRcate. Unlike DMPs which test individual
#              CpG sites, DMRs identify genomic regions where
#              multiple nearby CpGs show coordinated methylation
#              differences between groups.
# Input:  data/processed/dmp_results.rds
# Output: results/differential/dmr_results.csv
#         data/processed/dmr_results.rds
# Usage:  Rscript scripts/09_dmr_analysis.R
# =============================================================

library(minfi)
library(limma)
library(DMRcate)

# Set working directory
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
setwd(project_dir)

# Load DMP results from previous step
# We need m_comparison, design, contrast_matrix and sample_info_comparison
cat("Loading DMP results...\n")
dmp_data <- readRDS("data/processed/dmp_results.rds")
m_comparison           <- dmp_data$m_comparison
beta_comparison        <- dmp_data$beta_comparison
sample_info_comparison <- dmp_data$sample_info_comparison

cat("Samples:", ncol(m_comparison), "\n")
cat("Probes:", nrow(m_comparison), "\n")

# Recreate design and contrast matrices
# These must match exactly what was used in Step 8
cat("\nRecreating design and contrast matrices...\n")
design <- model.matrix(~ 0 + Comparison_Group, data = sample_info_comparison)
colnames(design) <- gsub("Comparison_Group", "", colnames(design))

contrast_matrix <- makeContrasts(
    LNCAP_vs_PREC = paste0(colnames(design)[1], "-", colnames(design)[2]),
    levels = design
)
cat("Contrast: LNCAP vs PREC\n")

# ── DMRcate annotation ────────────────────────────────────────────────────────
# cpg.annotate() runs the same limma model as Step 8 but formats the
# results in a way DMRcate understands. It attaches statistical information
# (t-statistics, p-values) to each CpG's genomic coordinates so that
# nearby significant CpGs can be grouped into regions.
cat("\nAnnotating CpGs for DMR analysis...\n")
cat("This fits the linear model and prepares spatial correlation...\n")

dmr_annotation <- cpg.annotate(
    object        = m_comparison,         # M-values for statistical testing
    datatype      = "array",              # Array-based methylation data
    what          = "M",                  # Use M-values
    arraytype     = "EPICv2",             # EPIC v2.0 array
    analysis.type = "differential",       # Differential methylation
    design        = design,               # Design matrix
    contrasts     = TRUE,                 # Use contrast matrix
    cont.matrix   = contrast_matrix,      # Contrast matrix
    coef          = colnames(contrast_matrix)[1],  # Which contrast to test
    fdr           = 0.05                  # FDR threshold
)

cat("Annotation complete!\n")
cat("Significant CpGs for DMR analysis:",
    sum(dmr_annotation@ranges$is.sig), "\n")

# ── DMRcate region finding ────────────────────────────────────────────────────
# dmrcate() groups nearby significant CpGs into regions using a
# Gaussian kernel smoother.
# lambda = 1000: smoothing bandwidth in base pairs — CpGs within
#                1000bp of each other are considered for grouping
# C = 2:         scaling factor for lambda
# pcutoff = 0.05: only include regions with p < 0.05
#
# Think of it like this: instead of asking "is THIS CpG different?"
# DMRcate asks "is this NEIGHBORHOOD of CpGs consistently different?"
# This is more robust because it requires coordinated changes across
# multiple nearby sites rather than a single CpG fluctuation.
cat("\nFinding differentially methylated regions...\n")
cat("Parameters: lambda=1000bp, C=2, pcutoff=0.05\n")

dmrs <- dmrcate(dmr_annotation, lambda = 1000, C = 2, pcutoff = 0.05)

cat("DMR identification complete!\n")
cat("Total DMRs found:", length(dmrs@coord), "\n")

# ── Extract genomic ranges ────────────────────────────────────────────────────
# Convert DMRcate results to a GRanges object with genomic coordinates
# then to a data frame for saving and downstream analysis
cat("\nExtracting genomic ranges...\n")
dmr_results <- extractRanges(dmrs, genome = "hg38")
dmr_df      <- as.data.frame(dmr_results)

# Print summary statistics
cat("\n── DMR Summary ─────────────────────────────────────\n")
cat("Total DMRs identified:      ", nrow(dmr_df), "\n")
cat("Median DMR width (bp):      ", median(dmr_df$width), "\n")
cat("Median CpGs per DMR:        ", median(dmr_df$no.cpgs), "\n")
cat("DMRs with >= 5 CpGs:        ", sum(dmr_df$no.cpgs >= 5), "\n")
cat("Hypermethylated DMRs (cancer):", sum(dmr_df$meandiff > 0), "\n")
cat("Hypomethylated DMRs (cancer): ", sum(dmr_df$meandiff < 0), "\n")
cat("────────────────────────────────────────────────────\n")

# Show top 10 DMRs by significance
cat("\nTop 10 DMRs by significance:\n")
top_dmrs <- dmr_df[order(dmr_df$Stouffer), ][1:10, 
    c("seqnames", "start", "end", "width", "no.cpgs", 
      "meandiff", "Stouffer", "HMFDR", "overlapping.genes")]
print(top_dmrs)

# ── Save results ──────────────────────────────────────────────────────────────
cat("\nSaving DMR results...\n")
write.csv(dmr_df, "results/differential/dmr_results.csv", row.names = FALSE)

saveRDS(list(
    dmr_df         = dmr_df,
    dmr_results    = dmr_results,
    dmrs           = dmrs,
    dmr_annotation = dmr_annotation
), "data/processed/dmr_results.rds")

cat("Results saved to results/differential/dmr_results.csv\n")
cat("\nDMR analysis complete!\n")
