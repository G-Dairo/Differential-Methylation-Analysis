#!/bin/bash
# Run this ONCE before submitting any SLURM jobs
# Sets up conda environment with R 4.5.3 and core Bioconductor packages

conda create -n methylation_env -c conda-forge r-base=4.5.3

conda install -n methylation_env -c conda-forge -c bioconda \
    bioconductor-genomicranges \
    bioconductor-summarizedexperiment \
    bioconductor-annotationdbi \
    bioconductor-biostrings \
    bioconductor-xvector \
    bioconductor-rsamtools \
    bioconductor-genomicalignments \
    bioconductor-genomicfeatures \
    bioconductor-rtracklayer \
    bioconductor-delayedarray \
    bioconductor-hdf5array \
    bioconductor-delayedmatrixstats \
    bioconductor-geoquery \
    bioconductor-bumphunter \
    bioconductor-genefilter \
    bioconductor-annotate \
    bioconductor-minfi \
    r-xml \
    bioconductor-dmrcate
