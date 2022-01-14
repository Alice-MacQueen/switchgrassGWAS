## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/gwas-",
  fig.width = 6,
  fig.asp = 0.618,
  out.width = "70%",
  fig.align = "center"
)

## -----------------------------------------------------------------------------
# get the example bedfile from the package switchgrassGWAS
bedfile <- system.file("extdata", "example.bed", package = "switchgrassGWAS")

## ----setup--------------------------------------------------------------------
# Load packages bigsnpr and switchgrassGWAS
library(switchgrassGWAS)
library(bigsnpr)

# Read from bed/bim/fam to create the new files that bigsnpr uses.
# Let's put them in an temporary directory for this demo.
tmpfile <- tempfile()
snp_readBed(bedfile, backingfile = tmpfile)

# Attach the "bigSNP" object to the R session.
snp_example <- snp_attach(paste0(tmpfile, ".rds"))
# What does the bigSNP object look like?
str(snp_example, max.level = 2, strict.width = "cut")

# Load the pvdiv phenotypes into the R session.
data(pvdiv_phenotypes) 
# Make an example dataframe of one phenotype where the first column is PLANT_ID.
# This "phenotype", 'GWAS_CT', is the number of times a plant successfully 
# clonally replicated to plant in the common gardens in 2018.
one_phenotype <- pvdiv_phenotypes %>%
  dplyr::select(PLANT_ID, GWAS_CT)

## -----------------------------------------------------------------------------
# Save the output to a temporary directory for this demo.
tempdir <- tempdir()

pvdiv_standard_gwas(snp = snp_example, df = one_phenotype, 
                    type = "linear", outputdir = tempdir, 
                    savegwas = FALSE, saveplots = TRUE,
                    saveannos = FALSE, ncores = 1)

