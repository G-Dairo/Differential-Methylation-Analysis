# =============================================================
# Script: 11_enrichment_visualization.R
# Description: Functional enrichment visualization for DNA
#              methylation analysis. Uses clusterProfiler for
#              GSEA and ORA approaches on GO and KEGG databases.
#              Adapted from NGS101 enrichment tutorial for
#              human prostate cancer methylation data.
# Input:  data/processed/dmp_results.rds (results_annotated)
# Output: plots/go_gsea_dotplot.pdf
#         plots/kegg_gsea_dotplot.pdf
#         plots/gsea_top_kegg_pathway.pdf
#         plots/go_ora_hyper_dotplot.pdf
#         plots/go_ora_hypo_dotplot.pdf
#         results/annotation/go_gsea_results.csv
#         results/annotation/kegg_gsea_results.csv
#         results/annotation/go_ora_hyper_results.csv
#         results/annotation/go_ora_hypo_results.csv
# Usage:  Rscript scripts/11_enrichment_visualization.R
# =============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(DOSE)
library(enrichplot)

# Set working directory
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
setwd(project_dir)

# ── Load data ─────────────────────────────────────────────────────────────────
# Use results_annotated which contains both limma statistics
# AND genomic annotation including gene symbols (UCSC_RefGene_Name)
cat("Loading DMP annotated results...\n")
dmp_data    <- readRDS("data/processed/dmp_results.rds")
results_all <- dmp_data$results_annotated
cat("Total CpG sites:", nrow(results_all), "\n")
cat("Columns available:", paste(colnames(results_all), collapse = ", "), "\n")

# ── Extract gene symbols from annotation ─────────────────────────────────────
# UCSC_RefGene_Name can contain multiple genes per probe separated by
# semicolons e.g. "TP53;TP53" or "BRCA1;BRCA2"
# We take the first gene for each probe
cat("\nExtracting gene symbols from UCSC_RefGene_Name...\n")
gene_symbols <- sapply(results_all$UCSC_RefGene_Name, function(x) {
    if (is.na(x) || x == "") return(NA)
    strsplit(as.character(x), ";")[[1]][1]
    genes <- strsplit(as.character(x), ";")[[1]]
    genes <- unique(genes)                    
    genes <- genes[genes != ""]           
    if (length(genes) == 0) return(NA)
    return(genes[1]) 
})
cat("Probes with gene annotation:",
    sum(!is.na(gene_symbols) & gene_symbols != ""), "\n")
cat("Probes without gene annotation:",
    sum(is.na(gene_symbols) | gene_symbols == ""), "\n")
cat("Example gene symbols:\n")
print(head(unique(na.omit(gene_symbols)), 10))

# ── Prepare ranked gene list for GSEA ────────────────────────────────────────
# GSEA requires ALL genes ranked by logFC — not just significant ones.
# Positive logFC = hypermethylated in LNCAP (cancer)
# Negative logFC = hypomethylated in LNCAP (cancer)
cat("\nPreparing ranked gene list...\n")
results_sorted        <- results_all[order(-results_all$logFC), ]
logfc                 <- results_sorted$logFC
names(logfc)          <- gene_symbols[rownames(results_sorted)]

# Remove probes with no gene annotation and duplicates
# Keep only one probe per gene (highest absolute logFC already at top
# after sorting)
logfc <- logfc[!is.na(names(logfc)) & names(logfc) != ""]
logfc <- logfc[!duplicated(names(logfc))]
cat("Total unique genes in ranked list:", length(logfc), "\n")

# Convert gene symbols to Entrez IDs for KEGG
# KEGG requires Entrez IDs rather than gene symbols
cat("Converting to Entrez IDs for KEGG...\n")
entrez_ids <- mapIds(
    x         = org.Hs.eg.db,
    keys      = names(logfc),
    keytype   = "SYMBOL",
    column    = "ENTREZID",
    multiVals = "first"
)
logfc_entrez        <- logfc
names(logfc_entrez) <- entrez_ids
logfc_entrez        <- logfc_entrez[!is.na(names(logfc_entrez))]
logfc_entrez        <- logfc_entrez[!duplicated(names(logfc_entrez))]
cat("Genes mapped to Entrez IDs:", length(logfc_entrez), "\n")

# ── Part 1: GO GSEA ───────────────────────────────────────────────────────────
# GSEA tests whether genes associated with methylation changes are
# systematically enriched at the top or bottom of the ranked list.
# Captures coordinated biological signals rather than just counting
# significant genes — more sensitive than ORA for subtle effects.
cat("\n── Part 1: GO GSEA ─────────────────────────────────────────────────────\n")
cat("Running GO GSEA (this may take a few minutes)...\n")

enrich_go_gsea <- gseGO(
    geneList     = logfc,
    OrgDb        = org.Hs.eg.db,
    ont          = "ALL",
    pvalueCutoff = 0.05,
    keyType      = "SYMBOL",
    verbose      = FALSE
)

enrich_go_gsea_df <- enrich_go_gsea@result
cat("Significant GO terms:", nrow(enrich_go_gsea_df), "\n")

write.csv(enrich_go_gsea_df,
          "results/annotation/go_gsea_results.csv",
          row.names = FALSE)
cat("Saved: results/annotation/go_gsea_results.csv\n")

if (nrow(enrich_go_gsea_df) > 0) {
    cat("Generating GO GSEA dotplot...\n")
    dotplot_go_gsea <- dotplot(enrich_go_gsea,
                               showCategory = 15,
                               orderBy      = "GeneRatio",
                               split        = ".sign") +
        facet_grid(. ~ .sign) +
        ggtitle("GO Enrichment: LNCAP vs PREC Methylation") +
        theme(plot.title = element_text(hjust = 0.5, size = 12))

    ggsave("plots/go_gsea_dotplot.pdf",
           dotplot_go_gsea,
           device = "pdf",
           units  = "cm",
           width  = 28,
           height = 22)
    cat("Saved: plots/go_gsea_dotplot.pdf\n")
} else {
    cat("No significant GO terms found for GSEA dotplot\n")
}

# ── Part 2: KEGG GSEA ────────────────────────────────────────────────────────
# KEGG GSEA identifies dysregulated signaling and metabolic pathways.
# Uses Entrez IDs because KEGG requires them for pathway mapping.
cat("\n── Part 2: KEGG GSEA ───────────────────────────────────────────────────\n")
cat("Running KEGG GSEA...\n")

enrich_kegg_gsea <- gseKEGG(
    geneList     = logfc_entrez,
    organism     = "hsa",
    pvalueCutoff = 0.05,
    verbose      = FALSE
)

enrich_kegg_gsea_df <- enrich_kegg_gsea@result
cat("Significant KEGG pathways:", nrow(enrich_kegg_gsea_df), "\n")

write.csv(enrich_kegg_gsea_df,
          "results/annotation/kegg_gsea_results.csv",
          row.names = FALSE)
cat("Saved: results/annotation/kegg_gsea_results.csv\n")

if (nrow(enrich_kegg_gsea_df) > 0) {
    cat("Generating KEGG GSEA dotplot...\n")
    dotplot_kegg_gsea <- dotplot(enrich_kegg_gsea,
                                 showCategory = 15,
                                 orderBy      = "GeneRatio",
                                 split        = ".sign") +
        facet_grid(. ~ .sign) +
        ggtitle("KEGG Pathway Enrichment: LNCAP vs PREC Methylation") +
        theme(plot.title = element_text(hjust = 0.5, size = 12))

    ggsave("plots/kegg_gsea_dotplot.pdf",
           dotplot_kegg_gsea,
           device = "pdf",
           units  = "cm",
           width  = 28,
           height = 22)
    cat("Saved: plots/kegg_gsea_dotplot.pdf\n")

    # GSEA plot for the single most significant KEGG pathway
    top_pathway <- enrich_kegg_gsea_df$ID[1]
    top_name    <- enrich_kegg_gsea_df$Description[1]
    cat("Generating GSEA plot for top pathway:", top_name, "\n")

    gsea_plot <- gseaplot2(enrich_kegg_gsea,
                           geneSetID = top_pathway,
                           title     = top_name)

    ggsave("plots/gsea_top_kegg_pathway.pdf",
           gsea_plot,
           device = "pdf",
           units  = "cm",
           width  = 20,
           height = 16)
    cat("Saved: plots/gsea_top_kegg_pathway.pdf\n")
} else {
    cat("No significant KEGG pathways found\n")
}

# ── Part 3: GO ORA — directional analysis ────────────────────────────────────
# ORA separately tests hypermethylated and hypomethylated genes.
# This reveals whether silenced genes (hypermethylated promoters) and
# activated genes (hypomethylated promoters) affect different pathways —
# a key biological question in cancer epigenetics.
cat("\n── Part 3: GO ORA (Directional) ────────────────────────────────────────\n")

# Extract gene symbols for hypermethylated probes
hyper_genes <- unique(na.omit(gene_symbols[
    rownames(results_all[!is.na(results_all$adj.P.Val) &
                         results_all$adj.P.Val < 0.05  &
                         results_all$logFC > 1, ])
]))
hyper_genes <- hyper_genes[hyper_genes != ""]

# Extract gene symbols for hypomethylated probes
hypo_genes  <- unique(na.omit(gene_symbols[
    rownames(results_all[!is.na(results_all$adj.P.Val) &
                         results_all$adj.P.Val < 0.05  &
                         results_all$logFC < -1, ])
]))
hypo_genes <- hypo_genes[hypo_genes != ""]

cat("Unique genes hypermethylated (FDR<0.05, logFC>1):", length(hyper_genes), "\n")
cat("Unique genes hypomethylated (FDR<0.05, logFC<-1):", length(hypo_genes), "\n")

# ORA for hypermethylated genes
cat("\nRunning GO ORA for hypermethylated genes...\n")
enrich_go_hyper <- enrichGO(
    gene          = hyper_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "ALL",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
)

if (!is.null(enrich_go_hyper) && nrow(enrich_go_hyper@result) > 0) {
    cat("Significant GO terms (hypermethylated):",
        nrow(enrich_go_hyper@result), "\n")

    dotplot_hyper <- dotplot(enrich_go_hyper,
                             showCategory = 15,
                             orderBy      = "GeneRatio") +
        ggtitle("GO ORA: Hypermethylated in LNCAP (Cancer)") +
        theme(plot.title = element_text(hjust = 0.5, size = 12))

    ggsave("plots/go_ora_hyper_dotplot.pdf",
           dotplot_hyper,
           device = "pdf",
           units  = "cm",
           width  = 20,
           height = 20)
    cat("Saved: plots/go_ora_hyper_dotplot.pdf\n")

    write.csv(enrich_go_hyper@result,
              "results/annotation/go_ora_hyper_results.csv",
              row.names = FALSE)
    cat("Saved: results/annotation/go_ora_hyper_results.csv\n")
} else {
    cat("No significant GO terms for hypermethylated genes\n")
}

# ORA for hypomethylated genes
cat("\nRunning GO ORA for hypomethylated genes...\n")
enrich_go_hypo <- enrichGO(
    gene          = hypo_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "ALL",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
)

if (!is.null(enrich_go_hypo) && nrow(enrich_go_hypo@result) > 0) {
    cat("Significant GO terms (hypomethylated):",
        nrow(enrich_go_hypo@result), "\n")

    dotplot_hypo <- dotplot(enrich_go_hypo,
                            showCategory = 15,
                            orderBy      = "GeneRatio") +
        ggtitle("GO ORA: Hypomethylated in LNCAP (Cancer)") +
        theme(plot.title = element_text(hjust = 0.5, size = 12))

    ggsave("plots/go_ora_hypo_dotplot.pdf",
           dotplot_hypo,
           device = "pdf",
           units  = "cm",
           width  = 20,
           height = 20)
    cat("Saved: plots/go_ora_hypo_dotplot.pdf\n")

    write.csv(enrich_go_hypo@result,
              "results/annotation/go_ora_hypo_results.csv",
              row.names = FALSE)
    cat("Saved: results/annotation/go_ora_hypo_results.csv\n")
} else {
    cat("No significant GO terms for hypomethylated genes\n")
}

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n── Enrichment Visualization Complete ───────────────────────────────────\n")
cat("Plots saved to plots/:\n")
cat("  go_gsea_dotplot.pdf\n")
cat("  kegg_gsea_dotplot.pdf\n")
cat("  gsea_top_kegg_pathway.pdf\n")
cat("  go_ora_hyper_dotplot.pdf\n")
cat("  go_ora_hypo_dotplot.pdf\n")
cat("Results saved to results/annotation/\n")
cat("────────────────────────────────────────────────────────────────────────\n")
