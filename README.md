
<!-- README.md is generated from README.Rmd. Please edit that file -->
switchgrassGWAS
===============

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/Alice-MacQueen/switchgrassGWAS.svg?branch=master)](https://travis-ci.org/Alice-MacQueen/switchgrassGWAS) <!-- badges: end -->

The goal of switchgrassGWAS is to allow fast, powerful genome-wide association analysis on the Panicum virgatum diversity panel.

The switchgrass (Panicum virgatum) diversity panel is now being grown at many locations across the United States and Mexico. Many researchers are measuring phenotypes on this panel to understand the genes and genetic regions affecting these phenotypes. This package provides the code for fast, less memory intensive genome-wide association (GWAS) using bigsnpr. It also provides functions to link diversity panel phenotypic data with SNP data, to prepare basic plots in ggplot for further customization, and to prepare multiple GWAS results for use in the downstream application mash.

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Alice-MacQueen/switchgrassGWAS")
```

This will give you access to the code, example phenotypes, and the currently available information about the genotypes in the switchgrass diversity panel.

At this moment, not all of the genomic information is publicly available. If you would like access to this information pre-publication, please contact Tom Juenger at tjuenger AT utexas DOT edu.

To use this package prior to publication of the switchgrass genomic data, you'll need to request access to the genomic information, then download this information and put it in your working directory. Then, you can load the information into your R environment using the following commands:

``` r
library(bigsnpr)
library(AnnotationDbi)

snp <- snp_attach("Pvirgatum_4x_784g_imp_phased_maf0.02_QC.rds")
load("svd0.rda")
txdb <- loadDb(file = "Pvirgatum_516_v5.1.gene.txdb.sqlite")
```

Genome-Wide Association
-----------------------

You can use `pvdiv_gwas` to run linear or logistic univariate regression on 21 million SNPs with a minor allele frequency of 2% or higher. The following example demonstrates running a genome-wide association on a continuous trait via linear regression. The trait is tiller count at the end of the 2018 season in Brookings, South Dakota.

``` r
NCORES <- nb_cores()

one_phenotype <- data(phenotypes) %>%
  dplyr::select(PLANT_ID, BRKG_TC_EOS_2018)

gwas <- pvdiv_gwas(df = one_phenotype, type = "linear",
                     snp = snp, covar = svd0, ncores = NCORES)
```

You can then plot the results of this GWAS using built in functions in bigsnpr.

Note that while the GWAS will run quickly, the plotting functions take about ten times as long to run, because of the large number of datapoints to plot (&gt;21 million).

``` r
snp_manhattan(gwas, infos.chr = CHRN$CHRN, infos.pos = POS)
snp_qq(gwas)
```
