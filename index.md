
<!-- index.md is generated from index.Rmd. Please edit that file -->

# switchgrassGWAS

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/Alice-MacQueen/switchgrassGWAS.svg?branch=master)](https://travis-ci.org/Alice-MacQueen/switchgrassGWAS)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/Alice-MacQueen/switchgrassGWAS?branch=master&svg=true)](https://ci.appveyor.com/project/Alice-MacQueen/switchgrassGWAS)
<!-- badges: end -->

The R package <b>switchgrassGWAS</b> provides functions for genome-wide
association analysis on the *<b>P</b>anicum <b>v</b>irgatum*
<b>div</b>ersity (<b>pvdiv</b>) panel.

## Installation

You can install the development version from
[GitHub](https://github.com/Alice-MacQueen/switchgrassGWAS) from within
R:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("Alice-MacQueen/switchgrassGWAS")
```

This will give you access to the [package
functions](https://alice-macqueen.github.io/switchgrassGWAS/docs/reference/index.html),
example and previously published
[phenotypes](https://alice-macqueen.github.io/switchgrassGWAS/reference/phenotypes.html),
and the currently available information about the
[genotypes](https://alice-macqueen.github.io/switchgrassGWAS/reference/pvdiv_metadata.html)
in the switchgrass diversity panel.

### Installations for additional functionality

Some <b>switchgrassGWAS</b> functions require the installation of
additional packages.

These packages can be installed from within R with:

``` r
install.packages("bigsnpr")
install.packages("mashr")
devtools::install_github("lcolladotor/dots")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("curl")
BiocManager::install(c("GenomicFeatures", "VariantAnnotation"))
```

## Background

The switchgrass (*Panicum virgatum*) diversity panel, the <b>pvdiv</b>
panel, is being grown at many common gardens across the United States
and Mexico. Many researchers are measuring phenotypes on this panel to
understand the genes and genetic regions affecting these phenotypes.

This package provides the code for fast, less memory intensive
genome-wide association (GWAS) using
[bigsnpr](https://privefl.github.io/bigsnpr/). It also provides
functions to link diversity panel phenotypic data with publicly
available [SNP data](https://doi.org/10.18738/T8/J377KE), to prepare
GWAS results plots using ggplot, and to prepare multiple univariate GWAS
results for use in the downstream application mash.

## Usage

For a start, have a look at the code examples provided for
[`pvdiv_standard_gwas`](https://alice-macqueen.github.io/switchgrassGWAS/reference/pvdiv_standard_gwas.html).

Download the SNP data [here](https://doi.org/10.18738/T8/J377KE).

Look at the
[metadata](https://alice-macqueen.github.io/switchgrassGWAS/reference/pvdiv_metadata.html)
for genotypes in the diversity panel, and the publicly available
[phenotypes](https://alice-macqueen.github.io/switchgrassGWAS/reference/phenotypes.html).

## Documentation

The HTML documentation of the development version is available on
[Github](https://alice-macqueen.github.io/switchgrassGWAS/).
