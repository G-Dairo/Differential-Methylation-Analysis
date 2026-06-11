# DNA Methylation Analysis — EPIC v2.0 Array

## Overview
End-to-end analysis of DNA methylation data using Illumina EPIC v2.0 arrays,
following the NGS101 tutorial pipeline on the Amarel HPC cluster (Rutgers University).
Dataset: [GSE240469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240469) —
40 samples from prostate and breast cancer cell lines profiled on the EPIC v2.0 BeadChip.

## Background
DNA methylation is a chemical modification where a methyl group is added to cytosine
bases in the genome, typically at CpG sites. It plays a critical role in gene regulation
— heavily methylated regions are usually silenced, while unmethylated regions are
typically active. The Illumina EPIC v2.0 array measures methylation at over 900,000
CpG sites across the human genome, making it one of the most comprehensive tools
for genome-wide methylation profiling.

## Dataset
| Property | Detail |
|----------|--------|
| GEO Accession | [GSE240469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240469) |
| Array Platform | Illumina EPIC v2.0 (hg38) |
| Total Samples | 40 |
| Cell Lines | LNCaP, PrEC (prostate), MCF7, TAMR (breast) |
| Controls | Synthetic (SYN), Primary (FD) |
| Publication | [BMC Genomics 2024](https://link.springer.com/article/10.1186/s12864-024-10027-5) |

## Pipeline

### Step 0: Environment Setup
**Scripts:** `environment/conda_setup.sh` → `scripts/run_00_setup.sh`

The EPIC v2.0 analysis requires a specific set of R and Bioconductor packages.
Because the Amarel cluster only provides R 3.4.1 (too old), I use conda to
install R 4.5.3 in a personal environment. Core Bioconductor infrastructure
packages (minfi, GenomicRanges, etc.) are installed via conda to handle
system-level compilation dependencies, while higher-level packages are installed
via BiocManager on top.

**Note:** Having the correct package versions ensures reproducibility.
Bioconductor 3.21 and R 4.5.3 are pinned so that anyone running this pipeline
gets identical results regardless of when they run it.

```bash
# Run once to create the conda environment
bash environment/conda_setup.sh

# Then install R packages
sbatch slurm/run_00_setup.sh
```

---

### Step 1: Download Raw Data
**Script:** `scripts/01_download_data.sh`

Downloads the raw IDAT files for GSE240469 directly from NCBI GEO (~580MB).
IDAT is Illumina's proprietary binary format that stores the raw fluorescence
intensity signals from the array — one file for the Green channel (unmethylated)
and one for the Red channel (methylated) per sample, giving 80 files total for
40 samples.


```bash
sbatch scripts/01_download_data.sh
```

---

### Step 2: Create Sample Sheet
**Script:** `scripts/02_create_samplesheet.R` | **SLURM:** `slurm/run_02_create_samplesheet.sh`

Parses the IDAT filenames to extract sample metadata — cell line, DNA input
concentration, replicate number, and treatment condition. This produces a
structured CSV (`data/raw/sample_sheet.csv`) that links each sample to its
biological context. minfi requires this sample sheet to correctly pair Red
and Green IDAT files for each sample.

```bash
sbatch scripts/run_02_create_samplesheet.sh
```

**Output summary (40 samples):**

| Group | N | Description |
|-------|---|-------------|
| PREC_500/250/125 | 9 | Prostate epithelial, 3 DNA concentrations |
| LNCAP_500/250/125 | 9 | Prostate cancer, 3 DNA concentrations |
| TAMR_Aza / TAMR_Control | 7 | Breast cancer, decitabine treated vs untreated |
| FD | 8 | Primary tissue controls |
| SYN | 6 | Synthetic controls for technical QC |
| MCF7 | 1 | Breast cancer cell line |

---

### Step 3: Load Raw Data
**Script:** `scripts/03_load_data.R` | **SLURM:** `slurm/run_03_load_data.sh`

Reads all 80 IDAT files into R as an `RGChannelSetExtended` object using minfi.
This object stores the raw Red and Green fluorescence intensities for all
1,105,209 probes across all 40 samples, along with quality metrics including
signal standard deviations (`GreenSD`, `RedSD`) and bead counts (`NBeads`) —
extra information captured because we used `extended = TRUE`.

**Note:** The `RGChannelSetExtended` is the starting point for all
downstream analysis. Loading extended data gives us richer QC information that
helps identify unreliable probes before any biological analysis begins.

```bash
sbatch scripts/run_03_load_data.sh
```

**Result:** 1,105,209 probes × 40 samples loaded successfully.

---

### Step 4: Quality Control
**Script:** `scripts/04_quality_control.R` | **SLURM:** `slurm/run_04_quality_control.sh`

Performs quality control using detection p-values — a statistical measure of
whether each probe's signal is distinguishable from background noise. The
detection p-value is an empirical measure that compares each probe's signal
intensity against the distribution of negative control probes built into the
array. A probe fails if its p-value exceeds 0.01, meaning we cannot confidently
distinguish its signal from background.

Samples with more than 5% failed probes are removed, and probes failing in
more than 10% of samples are flagged for removal in the normalization step.

**Note:** Poor quality probes and samples introduce noise that can
lead to false discoveries in differential methylation analysis. Removing them
early ensures that downstream results reflect true biological signal rather
than technical artifacts.

```bash
sbatch scripts/run_04_quality_control.sh
```

**QC Results:**

| Metric | Result |
|--------|--------|
| Samples removed | 1 (GSM7698439_LNCAP_250_1 — poor hybridization) |
| Probes flagged | 8,331 out of 936,990 (< 1%) |
| Samples remaining | 39 |

---

### Step 5: Normalization
**Script:** `scripts/05_normalization.R` | **SLURM:** `scripts/run_05_normalization.sh`

Applies functional normalization (`preprocessFunnorm`) to correct for technical
variation introduced during sample preparation and array hybridization. This method
uses control probes built into the EPIC array to estimate and remove technical
effects including background fluorescence (`bgCorr = TRUE`) and dye bias —
the tendency for Red and Green channels to have different baseline intensities
(`dyeCorr = TRUE`). After normalization, failed probes identified in Step 4 are
removed, and methylation is extracted as two complementary value types:

- **Beta values** (0–1): the proportion of methylation at each CpG site.
  Intuitive and interpretable but statistically suboptimal due to
  heteroscedasticity at the extremes.
- **M values** (log2 ratio): better statistical properties for linear modelling
  and differential analysis, used in all downstream statistical tests.

**Note:** Without normalization, technical differences between samples
can masquerade as biological differences, leading to false discoveries. Functional
normalization is the recommended approach for datasets with mixed cell types or
multiple biological groups.

```bash
sbatch scripts/run_05_normalization.sh
```

> **HPC Note:** `preprocessCore` must be reinstalled with threading disabled on
> Amarel before running this step:
> ```bash
> R -e 'BiocManager::install("preprocessCore", configure.args = "--disable-threading", force = TRUE, ask = FALSE)'
> ```

**Result:** 921,370 probes × 39 samples. Beta values range: 0 to 0.99.

---

### Step 6: Probe Filtering
**Script:** `scripts/06_probe_filtering.R` | **SLURM:** `scripts/run_06_probe_filtering.sh`

Removes three categories of probes that could introduce noise or confounding
into downstream differential methylation analysis:

1. **Cross-reactive probes** — probes that can bind to multiple locations in
   the genome, meaning their signal is a mixture of methylation from different
   regions rather than a single CpG site.

2. **SNP-affected probes** — probes where a Single Nucleotide Polymorphism
   (SNP) sits at or near the CpG being measured. Since SNPs change the DNA
   sequence, they affect probe binding and make the signal reflect genetic
   variation rather than true methylation differences. Three SNP positions
   are checked: within the probe body (`Probe_rs`), at the CpG site itself
   (`CpG_rs`), and at the single base extension position (`SBE_rs`).

3. **Sex chromosome probes** — probes on chrX and chrY behave differently
   between males and females due to X-chromosome inactivation and Y-chromosome
   absence in females. Removing them prevents biological sex from confounding
   the analysis.

Finally, probes with missing or non-finite M-values are removed to ensure
complete data for statistical testing.

**Note:** Including these probes would introduce false positives
in differential methylation analysis — apparent methylation differences that
reflect technical artifacts or genetic variation rather than true epigenetic
changes.

```bash
sbatch slurm/run_06_probe_filtering.sh
```

**Filtering Summary:**

| Filter | Probes Removed |
|--------|---------------|
| Cross-reactive / SNP-affected | 180,909 |
| Sex chromosomes (chrX, chrY) | 23,538 |
| Incomplete / non-finite | 19 |
| **Final clean probes** | **719,922** |

---

### Step 7: Sample Outlier Detection
**Script:** `scripts/07_outlier_detection.R` | **SLURM:** `scripts/run_07_outlier_detection.sh`

Performs a final sample-level quality check before differential methylation
analysis using two complementary approaches:

**Principal Component Analysis (PCA)** reduces the 719,922-dimensional
methylation data to a small number of components that capture the major
patterns of variation. Samples are plotted in PC space — samples that cluster
far from their biological group are candidates for removal. The scree plot
shows how much biological structure is captured by the top PCs.

**Mahalanobis distance** provides a statistical measure of how far each sample
sits from the group centre in the top 5 PC space, accounting for correlations
between components. Samples beyond the 97.5th percentile of a chi-squared
distribution (threshold = 12.83) are flagged as outliers.

A **sample-to-sample correlation heatmap** is also generated using beta values,
providing a visual overview of how similar samples are to each other — samples
of the same cell line should cluster together with high correlation.

**Note:** Outlier samples that passed earlier QC filters can still
have aberrant methylation patterns that distort group comparisons. Detecting
them before statistical testing prevents a single bad sample from generating
false discoveries across thousands of CpG sites.

```bash
sbatch scripts/run_07_outlier_detection.sh
```

**Results:**

| Metric | Result |
|--------|--------|
| Outliers detected | 0 |
| Mahalanobis threshold | 12.83 |
| PC1 variance explained | 32.4% |
| PC2 variance explained | 19.7% |
| Top 5 PCs combined | 85.8% |
| Final samples | 39 |
| Final probes | 719,922 |

---

### Step 8: Differential Methylation Analysis
**Script:** `scripts/08_differential_methylation.R` | **SLURM:** `scripts/run_08_differential_methylation.sh`

Identifies Differentially Methylated Positions (DMPs) between LNCAP
(prostate cancer) and PREC (normal prostate epithelial) cell lines using
**limma** — a linear modelling framework widely adopted for methylation arrays.

For each of the 719,922 CpG sites, limma fits a linear model where:
- **Response variable (Y):** M-value at each CpG site
- **Explanatory variable (X):** Group membership (LNCAP = 1, PREC = 0)

The **empirical Bayes** moderation borrows information across all 719,922
probes to stabilize variance estimates — giving more reliable results than
a standard t-test especially with small sample sizes (n=9 per group).

A contrast matrix explicitly defines the comparison:
- Positive logFC = hypermethylated in LNCAP (cancer)
- Negative logFC = hypomethylated in LNCAP (cancer)

M-values are used for statistical testing while beta values are retained
for biological interpretation and visualization.

```bash
sbatch slurm/run_08_differential_methylation.sh
```

**Results:**

| Metric | Result |
|--------|--------|
| CpG sites tested | 719,922 |
| Significant DMPs (FDR < 0.05) | 476,854 (66%) |
| Significant DMPs (FDR < 0.05, \|M\| > 1) | 328,376 |
| Hypermethylated in LNCAP (cancer) | 196,778 (41%) |
| Hypomethylated in LNCAP (cancer) | 280,076 (59%) |

The predominance of hypomethylation in cancer is consistent with the
classical cancer epigenetics paradigm — global hypomethylation drives
genomic instability and oncogene activation, while focal hypermethylation
silences tumour suppressor genes.

---

### Step 9: Differentially Methylated Regions
**Script:** `scripts/09_dmr_analysis.R` | **SLURM:** `scripts/run_09_dmr_analysis.sh`

Identifies Differentially Methylated Regions (DMRs) using **DMRcate**.
Unlike DMPs which test individual CpG sites independently, DMR analysis
groups nearby CpGs that show coordinated methylation differences into
genomic regions — providing a more biologically interpretable and
statistically robust view of epigenetic changes.

DMRcate works in two steps:
1. **`cpg.annotate()`** — runs the same limma model as Step 8 but
   attaches statistical information to each CpG's genomic coordinates
2. **`dmrcate()`** — applies a Gaussian kernel smoother to identify
   regions where nearby significant CpGs cluster together

Key parameters:
- `lambda = 1000`: smoothing bandwidth in base pairs
- `C = 2`: scaling factor for the bandwidth
- `pcutoff = 0.05`: only regions with p < 0.05 are reported

**Note:** DMRs require coordinated changes across multiple
nearby sites — making false positives far less likely than with
single-probe testing.

```bash
sbatch scripts/run_09_dmr_analysis.sh
```

**Results:**

| Metric | Result |
|--------|--------|
| Total DMRs identified | 102,937 |
| Median DMR width | 734 bp |
| Median CpGs per DMR | 3 |
| High confidence DMRs (≥ 5 CpGs) | 29,105 |
| Hypermethylated DMRs in cancer | 42,352 (41%) |
| Hypomethylated DMRs in cancer | 60,585 (59%) |

**Top biological finding:** The HOXA gene cluster (chr7) was the most
significantly differentially methylated region, showing strong
hypermethylation across 138 CpGs spanning 25kb. HOXA genes are
well-established tumour suppressors in prostate cancer — their silencing
through promoter hypermethylation is a classic cancer epigenetics finding
independently rediscovered by this pipeline.

---

### Step 10: Functional Enrichment Analysis
**Script:** `scripts/10_functional_enrichment.R` | **SLURM:** `scripts/run_10_functional_enrichment.sh`

Performs functional enrichment analysis to identify which biological
processes and pathways are affected by the methylation differences,
using two complementary approaches:

**CpG-based enrichment (`gometh`)** uses the `missMethyl` package which
is specifically designed for methylation array data. Standard enrichment
tools assume each gene is equally likely to be selected, but methylation
arrays have more probes in gene-rich regions and near CpG islands —
creating a bias. `gometh` corrects for this probe-count bias, making
results more reliable than a standard hypergeometric test.

**DMR-based gene extraction** takes a complementary region-based approach
by extracting all unique genes physically overlapping with DMRs. These
are the genes most likely to have their expression affected by the
methylation changes.

```bash
sbatch scripts/run_10_functional_enrichment.sh
```

**Results:**

| Analysis | Result |
|----------|--------|
| Significant KEGG pathways (FDR < 0.05) | 25 |
| Top pathway | Cadherin signaling |
| Other key pathways | Ras signaling, Pathways in cancer, Proteoglycans in cancer |
| Unique genes overlapping DMRs | 29,654 |
| Genes in hypermethylated DMRs | 20,248 |
| Genes in hypomethylated DMRs | 21,561 |
| High confidence DMR genes (≥ 5 CpGs) | 23,055 |

---

### Step 11: Enrichment Visualization
**Script:** `scripts/11_enrichment_visualization.R` | **SLURM:** `scripts/run_11_enrichment_visualization.sh`

Generates publication-quality visualizations of the functional enrichment
results using **clusterProfiler** with two complementary approaches:

**Gene Set Enrichment Analysis (GSEA)** uses the full ranked list of genes
ordered by methylation logFC — testing whether genes in a pathway are
systematically enriched at the top (hypermethylated) or bottom
(hypomethylated) of the ranking. More sensitive than ORA because it uses
the entire gene list rather than an arbitrary significance cutoff.

**Over-Representation Analysis (ORA)** separately tests hypermethylated
and hypomethylated gene sets against GO terms — revealing whether silenced
genes and activated genes affect different biological processes. This
directional separation is a key biological question in cancer epigenetics.

Both GO (Biological Process, Molecular Function, Cellular Component) and
KEGG pathway databases are tested. Gene symbols are extracted from the
`UCSC_RefGene_Name` annotation column and converted to Entrez IDs for
KEGG analysis.

```bash
sbatch scripts/run_11_enrichment_visualization.sh
```

**Output plots:**

| Plot | Description |
|------|-------------|
| `go_gsea_dotplot.pdf` | Top 15 GO terms by GSEA, split by direction |
| `kegg_gsea_dotplot.pdf` | Top 15 KEGG pathways by GSEA, split by direction |
| `gsea_top_kegg_pathway.pdf` | Detailed GSEA enrichment plot for top pathway |
| `go_ora_hyper_dotplot.pdf` | GO ORA for hypermethylated genes in cancer |
| `go_ora_hypo_dotplot.pdf` | GO ORA for hypomethylated genes in cancer |
