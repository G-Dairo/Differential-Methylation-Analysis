# =============================================================
# Script: 08_differential_methylation.R
# Description: Identify Differentially Methylated Positions (DMPs)
#              between LNCAP (prostate cancer) and PREC (normal
#              prostate epithelial) cell lines using limma.
# Input:  data/processed/final_data.rds
# Output: results/differential/dmp_results_annotated.csv
#         data/processed/dmp_results.rds
# Usage:  Rscript scripts/08_differential_methylation.R
# =============================================================

library(minfi)
library(limma)

# Set working directory
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
setwd(project_dir)

# Load final clean data from previous step
cat("Loading final data...\n")
final_data    <- readRDS("data/processed/final_data.rds")
m_complete    <- final_data$m_complete
beta_complete <- final_data$beta_complete
sample_info   <- final_data$sample_info
annotation    <- final_data$annotation

cat("Samples:", ncol(m_complete), "\n")
cat("Probes:", nrow(m_complete), "\n")

# ── Define comparison groups ──────────────────────────────────────────────────
# We compare LNCAP (prostate cancer) vs PREC (normal prostate epithelial)
# This is a biologically meaningful comparison — same tissue type but
# cancer vs normal, so differences reflect cancer-related methylation changes
cat("\nDefining comparison groups...\n")
sample_info$Comparison_Group <- ifelse(
    sample_info$Sample_Type %in% c("LNCAP", "PREC"),
    sample_info$Sample_Type,
    "Other"
)

# Filter to only include LNCAP and PREC samples
comparison_samples      <- sample_info$Comparison_Group %in% c("LNCAP", "PREC")
m_comparison            <- m_complete[, comparison_samples]
beta_comparison         <- beta_complete[, comparison_samples]
sample_info_comparison  <- sample_info[comparison_samples, ]

cat("LNCAP samples:", sum(sample_info_comparison$Sample_Type == "LNCAP"), "\n")
cat("PREC samples:", sum(sample_info_comparison$Sample_Type == "PREC"), "\n")
cat("Total samples in comparison:", ncol(m_comparison), "\n")

# ── Design matrix ─────────────────────────────────────────────────────────────
# The design matrix tells limma which samples belong to which group.
# ~ 0 + Comparison_Group means no intercept — each group gets its own
# coefficient, making contrasts easier to define and interpret.
cat("\nCreating design matrix...\n")
design <- model.matrix(~ 0 + Comparison_Group, data = sample_info_comparison)
colnames(design) <- gsub("Comparison_Group", "", colnames(design))
cat("Design matrix groups:", paste(colnames(design), collapse = " vs "), "\n")
cat("Design matrix dimensions:", nrow(design), "samples x", ncol(design), "groups\n")

# ── Contrast matrix ───────────────────────────────────────────────────────────
# The contrast defines exactly what comparison we want to make.
# LNCAP - PREC means: positive values = higher methylation in LNCAP (cancer)
#                     negative values = higher methylation in PREC (normal)
cat("Defining contrasts...\n")
contrast_matrix <- makeContrasts(
    LNCAP_vs_PREC = paste0(colnames(design)[1], "-", colnames(design)[2]),
    levels = design
)
cat("Contrast:", rownames(contrast_matrix)[1], "-",
    rownames(contrast_matrix)[2], "\n")

# ── Linear model fitting ──────────────────────────────────────────────────────
# limma fits a linear model to each of the 719,922 CpG sites simultaneously.
# eBayes applies empirical Bayes moderation — it borrows information across
# all probes to stabilize variance estimates, giving more reliable results
# than a standard t-test especially with small sample sizes.
cat("\nFitting linear models across", nrow(m_comparison), "CpG sites...\n")
fit           <- lmFit(m_comparison, design)
fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_eb        <- eBayes(fit_contrasts)
cat("Model fitting complete!\n")

# ── Extract results ───────────────────────────────────────────────────────────
cat("Extracting results...\n")
results_all <- topTable(fit_eb,
                        coef    = colnames(contrast_matrix),
                        n       = Inf,
                        p.value = 1,
                        lfc     = 0,
                        sort.by = "p")

cat("Total CpG sites tested:", nrow(results_all), "\n")

# Count significant DMPs at different thresholds
sig_fdr05  <- sum(results_all$adj.P.Val < 0.05, na.rm = TRUE)
sig_fdr01  <- sum(results_all$adj.P.Val < 0.01, na.rm = TRUE)
sig_lfc1   <- sum(results_all$adj.P.Val < 0.05 & abs(results_all$logFC) > 1, na.rm = TRUE)

cat("\n── DMP Summary ────────────────────────────────────\n")
cat("Significant DMPs (FDR < 0.05):              ", sig_fdr05, "\n")
cat("Significant DMPs (FDR < 0.01):              ", sig_fdr01, "\n")
cat("Significant DMPs (FDR < 0.05 & |M| > 1):   ", sig_lfc1, "\n")
cat("───────────────────────────────────────────────────\n")

# Direction of methylation changes
sig_results  <- results_all[results_all$adj.P.Val < 0.05, ]
hypermethylated <- sum(sig_results$logFC > 0, na.rm = TRUE)
hypomethylated  <- sum(sig_results$logFC < 0, na.rm = TRUE)
cat("Hypermethylated in LNCAP (cancer):", hypermethylated, "\n")
cat("Hypomethylated in LNCAP (cancer): ", hypomethylated, "\n")

# ── Add genomic annotation ────────────────────────────────────────────────────
# Merge statistical results with genomic coordinates and gene information
# so we know WHERE in the genome the methylation differences occur
cat("\nAdding genomic annotation...\n")
annotation_subset <- annotation[rownames(results_all), ]
results_annotated <- cbind(results_all, annotation_subset)

# ── Save results ──────────────────────────────────────────────────────────────
cat("Saving results...\n")
write.csv(results_annotated,
          "results/differential/dmp_results_annotated.csv",
          row.names = TRUE)

saveRDS(list(
    results_all       = results_all,
    results_annotated = results_annotated,
    fit_eb            = fit_eb,
    m_comparison      = m_comparison,
    beta_comparison   = beta_comparison,
    sample_info_comparison = sample_info_comparison
), "data/processed/dmp_results.rds")

cat("Results saved to results/differential/dmp_results_annotated.csv\n")
cat("\nDifferential methylation analysis complete!\n")
