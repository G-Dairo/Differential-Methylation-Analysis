# =============================================================
# Script: 10_functional_enrichment.R
# Description: Functional enrichment analysis of differentially
#              methylated positions and regions. Tests whether
#              genes near DMPs are enriched for specific biological
#              processes (GO) and pathways (KEGG), and extracts
#              genes overlapping DMRs for region-based enrichment.
# Input:  data/processed/dmp_results.rds
#         data/processed/dmr_results.rds
# Output: results/annotation/go_enrichment_results.csv
#         results/annotation/kegg_enrichment_results.csv
#         results/annotation/dmr_genes.csv
#         results/annotation/dmr_kegg_results.csv
#         data/processed/enrichment_results.rds
# Usage:  Rscript scripts/10_functional_enrichment.R
# =============================================================

library(missMethyl)
library(limma)

# Set working directory
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
setwd(project_dir)

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading DMP and DMR results...\n")
dmp_data <- readRDS("data/processed/dmp_results.rds")
dmr_data <- readRDS("data/processed/dmr_results.rds")

results_all <- dmp_data$results_all
dmr_df      <- dmr_data$dmr_df

cat("Total CpGs tested:", nrow(results_all), "\n")
cat("Significant DMPs (FDR<0.05):", sum(results_all$adj.P.Val < 0.05), "\n")
cat("Total DMRs:", nrow(dmr_df), "\n")

# ── Part 1: CpG-based enrichment (gometh) ────────────────────────────────────
# gometh is specifically designed for methylation array data.
# Standard enrichment tools assume each gene is equally likely to be
# selected, but methylation arrays have MORE probes in gene-rich regions
# and near CpG islands — creating a BIAS toward finding those genes enriched.
# gometh corrects for this probe-count bias, making results more reliable
# than using a standard hypergeometric test.
cat("\n── Part 1: CpG-based GO and KEGG Enrichment ────────────────────────────\n")

sig_cpgs <- rownames(results_all[results_all$adj.P.Val < 0.05, ])
all_cpgs <- rownames(results_all)
cat("Significant CpGs for enrichment:", length(sig_cpgs), "\n")
cat("Background CpGs:", length(all_cpgs), "\n")

# Gene Ontology enrichment
# Tests three GO categories:
# BP = Biological Process (e.g. "cell proliferation")
# MF = Molecular Function (e.g. "transcription factor binding")
# CC = Cellular Component (e.g. "nucleus")
cat("\nRunning GO enrichment analysis...\n")
go_results <- gometh(
    sig.cpg    = sig_cpgs,
    all.cpg    = all_cpgs,
    collection = "GO",
    array.type = "EPICv2"
)
go_results <- go_results[order(go_results$P.DE), ]

cat("GO terms tested:", nrow(go_results), "\n")
cat("Significant GO terms (FDR<0.05):", sum(go_results$FDR < 0.05, na.rm = TRUE), "\n")

cat("\nTop 10 GO terms:\n")
print(head(go_results[, c("ONTOLOGY", "TERM", "N", "DE", "P.DE", "FDR")], 10))

write.csv(go_results, "results/annotation/go_enrichment_results.csv", row.names = TRUE)
cat("GO results saved to results/annotation/go_enrichment_results.csv\n")

# KEGG pathway enrichment
# Tests whether genes near DMPs are overrepresented in specific
# biological pathways like "Prostate cancer", "DNA methylation",
# "Cell cycle" etc. Useful for identifying affected signaling networks.
cat("\nRunning KEGG pathway enrichment analysis...\n")
kegg_results <- gometh(
    sig.cpg    = sig_cpgs,
    all.cpg    = all_cpgs,
    collection = "KEGG",
    array.type = "EPICv2"
)
kegg_results <- kegg_results[order(kegg_results$P.DE), ]

cat("KEGG pathways tested:", nrow(kegg_results), "\n")
cat("Significant KEGG pathways (FDR<0.05):", sum(kegg_results$FDR < 0.05, na.rm = TRUE), "\n")

cat("\nTop 10 KEGG pathways:\n")
print(head(kegg_results[, c("ONTOLOGY", "TERM", "N", "DE", "P.DE", "FDR")], 10))

write.csv(kegg_results, "results/annotation/kegg_enrichment_results.csv", row.names = TRUE)
cat("KEGG results saved to results/annotation/kegg_enrichment_results.csv\n")

# ── Part 2: DMR-based gene extraction ────────────────────────────────────────
# While gometh works at the individual CpG level, here we take a
# complementary region-based approach. We extract all unique genes
# that physically overlap with our DMRs — these are the genes most
# likely to have their expression affected by the methylation changes.
cat("\n── Part 2: DMR Gene Extraction and Enrichment ──────────────────────────\n")

# Extract unique genes overlapping with DMRs
dmr_genes <- unique(unlist(strsplit(dmr_df$overlapping.genes, ", ")))
dmr_genes  <- dmr_genes[!is.na(dmr_genes) & dmr_genes != ""]

cat("Total unique genes overlapping DMRs:", length(dmr_genes), "\n")

# Split by methylation direction
hyper_dmrs <- dmr_df[dmr_df$meandiff > 0, ]
hypo_dmrs  <- dmr_df[dmr_df$meandiff < 0, ]

hyper_genes <- unique(unlist(strsplit(hyper_dmrs$overlapping.genes, ", ")))
hypo_genes  <- unique(unlist(strsplit(hypo_dmrs$overlapping.genes, ", ")))
hyper_genes <- hyper_genes[!is.na(hyper_genes) & hyper_genes != ""]
hypo_genes  <- hypo_genes[!is.na(hypo_genes)  & hypo_genes  != ""]

cat("Genes in hypermethylated DMRs:", length(hyper_genes), "\n")
cat("Genes in hypomethylated DMRs: ", length(hypo_genes), "\n")

# Save gene lists
dmr_gene_df <- data.frame(
    Gene      = dmr_genes,
    Direction = ifelse(dmr_genes %in% hyper_genes & dmr_genes %in% hypo_genes,
                       "Both",
                       ifelse(dmr_genes %in% hyper_genes,
                              "Hypermethylated", "Hypomethylated"))
)
dmr_gene_df <- dmr_gene_df[order(dmr_gene_df$Direction), ]
write.csv(dmr_gene_df, "results/annotation/dmr_genes.csv", row.names = FALSE)
cat("DMR gene list saved to results/annotation/dmr_genes.csv\n")

# High confidence DMR genes (DMRs with >= 5 CpGs only)
hc_dmrs      <- dmr_df[dmr_df$no.cpgs >= 5, ]
hc_genes     <- unique(unlist(strsplit(hc_dmrs$overlapping.genes, ", ")))
hc_genes     <- hc_genes[!is.na(hc_genes) & hc_genes != ""]
cat("High confidence DMR genes (>=5 CpGs):", length(hc_genes), "\n")

# ── Save all results ──────────────────────────────────────────────────────────
cat("\nSaving all enrichment results...\n")
saveRDS(list(
    go_results   = go_results,
    kegg_results = kegg_results,
    dmr_genes    = dmr_genes,
    hyper_genes  = hyper_genes,
    hypo_genes   = hypo_genes,
    hc_genes     = hc_genes
), "data/processed/enrichment_results.rds")

cat("\n── Enrichment Analysis Complete ────────────────────────────────────────\n")
cat("GO results:      results/annotation/go_enrichment_results.csv\n")
cat("KEGG results:    results/annotation/kegg_enrichment_results.csv\n")
cat("DMR genes:       results/annotation/dmr_genes.csv\n")
cat("────────────────────────────────────────────────────────────────────────\n")
