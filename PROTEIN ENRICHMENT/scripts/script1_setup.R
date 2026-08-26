
---

## Script 01 — `01_setup.R`

```r
# ============================================================
# 01_setup.R
# Project setup and package installation
# ============================================================

# CRAN packages
cran_packages <- c(
    "dplyr",
    "ggplot2",
    "igraph"
)

for (pkg in cran_packages) {

    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }

}

# Bioconductor
if (!requireNamespace(
    "BiocManager",
    quietly = TRUE
)) {
    install.packages("BiocManager")
}

bioc_packages <- c(
    "STRINGdb",
    "ComplexHeatmap"
)

for (pkg in bioc_packages) {

    if (!requireNamespace(
        pkg,
        quietly = TRUE
    )) {

        BiocManager::install(
            pkg,
            ask = FALSE,
            update = FALSE
        )

    }

}

# Load packages

library(STRINGdb)
library(dplyr)
library(igraph)
library(ggplot2)
library(ComplexHeatmap)

# Create directories

dir.create(
    "protein_enrichment",
    showWarnings = FALSE
)

dir.create(
    "protein_enrichment/input",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "protein_enrichment/results",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "protein_enrichment/figures",
    recursive = TRUE,
    showWarnings = FALSE
)

cat(
    "\nProtein enrichment environment initialized.\n"
)
