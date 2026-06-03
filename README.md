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
**Scripts:** `environment/conda_setup.sh` → `slurm/run_00_setup.sh`

The EPIC v2.0 analysis requires a specific set of R and Bioconductor packages.
Because the Amarel cluster only provides R 3.4.1 (too old), I use conda to
install R 4.5.3 in a personal environment. Core Bioconductor infrastructure
packages (minfi, GenomicRanges, etc.) are installed via conda to handle
system-level compilation dependencies, while higher-level packages are installed
via BiocManager on top.

```bash
# Run once to create the conda environment
bash environment/conda_setup.sh

# Then install R packages
sbatch slurm/run_00_setup.sh
```

---

### Step 1: Download Raw Data
**Script:** `slurm/01_download_data.sh`

Downloads the raw IDAT files for GSE240469 directly from NCBI GEO (~580MB).
IDAT is Illumina's proprietary binary format that stores the raw fluorescence
intensity signals from the array — one file for the Green channel (unmethylated)
and one for the Red channel (methylated) per sample, giving 80 files total for
40 samples

```bash
sbatch slurm/01_download_data.sh
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
sbatch slurm/run_02_create_samplesheet.sh
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

**Why this matters:** The `RGChannelSetExtended` is the starting point for all
downstream analysis. Loading extended data gives a richer QC information that
helps identify unreliable probes before any biological analysis begins.

```bash
sbatch slurm/run_03_load_data.sh
```

**Result:** 1,105,209 probes × 40 samples loaded successfully.

---

### Step 4: Quality Control
**Script:** `scripts/04_quality_control.R` | **SLURM:** `slurm/run_04_quality_control.sh`

Performs quality control using detection p-values — a statistical measure of
whether each probe's signal is distinguishable from background noise (i.e., empty probe).
A probe "fails" if its p-value exceeds 0.01, meaning we cannot confidently distinguish
its signal from background. Samples with more than 5% failed probes are removed,
and probes failing in more than 10% of samples are flagged for removal.


```bash
sbatch slurm/run_04_quality_control.sh
```

**QC Results:**
| Metric | Result |
|--------|--------|
| Samples removed | 1 (GSM7698439_LNCAP_250_1 — poor hybridization) |
| Probes flagged | 8,331 out of 936,990 (<1%) |
| Samples remaining | 39 |

## Requirements
- Conda (tested with conda 25.3.1)
- R 4.5.3 via conda-forge
- Bioconductor 3.21
- See `environment/conda_setup.sh` for full dependency list

## Outputs
All plots and reports are generated by running the scripts in order.
Final publication-quality figures are tracked in `plots/`.
Intermediate outputs and large data files are excluded from Git
but will appear locally after running the pipeline.

## Notes
- Raw IDAT files and processed `.rds` objects are not tracked in Git
  (too large). Run `slurm/01_download_data.sh` to download the raw data.
- All scripts are designed to run non-interactively via SLURM on the
  Amarel HPC cluster at Rutgers University.
- Set the `PROJECT_DIR` environment variable to change the project root path.
