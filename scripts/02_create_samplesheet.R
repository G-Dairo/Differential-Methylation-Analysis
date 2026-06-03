# =============================================================
# Script: 02_create_samplesheet.R
# Description: Parse IDAT filenames to create sample metadata
#              sheet for GSE240469 (40 samples, EPIC v2.0)
# Input:  data/raw/*.idat files
# Output: data/raw/sample_sheet.csv
# Usage:  Rscript scripts/02_create_samplesheet.R
# =============================================================

# Set working directory
project_dir <- Sys.getenv("PROJECT_DIR", unset = "/home/gsd67/Projects/Epigenetics/methylation")
setwd(project_dir)

# Extract sample information from filenames
# The file naming convention contains important sample metadata
idat_files <- list.files("data/raw", pattern = "*.idat", full.names = TRUE)
red_files <- idat_files[grepl("Red.idat", idat_files)]
basenames <- gsub("_Red.idat", "", red_files)

# Create sample information dataframe
sample_info <- data.frame(
    Basename = basenames,
    Sample_Name = basename(basenames),
    stringsAsFactors = FALSE
)

# Parse sample types from filenames (extract biological information)
sample_info$Sample_Type <- sapply(sample_info$Sample_Name, function(x) {
    if(grepl("PREC", x)) return("PREC")      # Prostate cancer cell line
    if(grepl("LNCAP", x)) return("LNCAP")    # Prostate cancer cell line
    if(grepl("SYN", x)) return("SYN")        # Synthetic control
    if(grepl("MCF7", x)) return("MCF7")      # Breast cancer cell line
    if(grepl("TAMR", x)) return("TAMR")      # Tamoxifen-resistant cell line
    if(grepl("FD", x)) return("FD")          # Control condition
    return("Unknown")
})

# Extract treatment information from filenames
sample_info$Treatment <- sapply(sample_info$Sample_Name, function(x) {
    if(grepl("500", x)) return("500")        # High concentration treatment
    if(grepl("250", x)) return("250")        # Medium concentration treatment
    if(grepl("125", x)) return("125")        # Low concentration treatment
    if(grepl("Aza", x)) return("Aza")        # 5-azacytidine treatment
    if(grepl("NoAza", x)) return("NoAza")    # No treatment control
    return("Control")
})

# Create simplified grouping for statistical analysis
sample_info$Group <- paste(sample_info$Sample_Type, sample_info$Treatment, sep = "_")
sample_info$Group[sample_info$Group == "FD_Control"] <- "FD"
sample_info$Group[sample_info$Group == "MCF7_Control"] <- "MCF7"

# Add technical variables that may affect data quality
sample_info$Slide <- "GSE240469"    # Experimental batch identifier
sample_info$Array <- "EPICv2"       # Array type

# Save sample sheet for future reference
write.csv(sample_info, "data/raw/sample_sheet.csv", row.names = FALSE)

# Print summary to log
cat("Sample sheet created successfully!\n")
cat("Total samples:", nrow(sample_info), "\n")
cat("\nSample type breakdown:\n")
print(table(sample_info$Sample_Type))
cat("\nGroup breakdown:\n")
print(table(sample_info$Group))
