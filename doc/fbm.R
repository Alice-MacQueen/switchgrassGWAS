## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/fbm-",
  fig.width = 6,
  fig.asp = 0.618,
  out.width = "70%",
  fig.align = "center"
)
library(knitr)

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
four_phenotype <- pvdiv_phenotypes %>%
  dplyr::select(PLANT_ID, GWAS_CT, matches("BIOMASS"))

## -----------------------------------------------------------------------------
# Save the output to a temporary directory for this demo.
tempdir <- tempdir()

pvdiv_standard_gwas(snp = snp_example, df = four_phenotype, 
                    type = "linear", outputdir = tempdir, 
                    savegwas = TRUE, savetype = "fbm", 
                    saveplots = FALSE,
                    saveannos = FALSE, ncores = 1)

## -----------------------------------------------------------------------------
gwas <- big_attach(file.path(tempdir, "gwas_effects_0M_SNPs.rds"))
gwas_metadata <- read.csv(file.path(tempdir, "gwas_effects_0M_SNPs_associated_metadata.csv"))
gwas_colnames <- read.csv(file.path(tempdir, "gwas_effects_0M_SNPs_column_names.csv"))

## -----------------------------------------------------------------------------
kable(gwas_metadata)

## -----------------------------------------------------------------------------
kable(gwas_colnames)

## -----------------------------------------------------------------------------
pvdiv_fbm_qq(effects = gwas, metadata = gwas_metadata, e_row = 2, o_row = 3:4, 
             outputdir = tempdir, thr = 0)

## -----------------------------------------------------------------------------
upset_df <- pvdiv_fbm_upset_df(effects = gwas, snp = snp_example, 
                               metadata = gwas_metadata, thr = 2)
upset_df

## -----------------------------------------------------------------------------
library(dplyr)

upset_df_1Mb <- upset_df %>%
  mutate(POS_binned = floor(POS/1000000)*1000000) %>%
  select(-POS) %>%
  select(CHR, POS_binned, everything()) %>%
  group_by(CHR, POS_binned) %>%
  summarise(across(everything(), ~ max(.x)))

## -----------------------------------------------------------------------------
sel_reg_biomass <- upset_df_1Mb %>%
  filter(across(CLMB_BIOMASS:PKLE_BIOMASS, ~ .x == 1)) %>%
  mutate(Sel_region = paste(CHR, POS_binned, sep = "_"))

sel_snp_biomass <- upset_df %>%
  mutate(POS_binned = floor(POS/1000000)*1000000,
         Sel_region = paste(CHR, POS_binned, sep = "_")) %>%
  filter(Sel_region %in% sel_reg_biomass$Sel_region)

kable(sel_snp_biomass)

## ---- include = FALSE---------------------------------------------------------
anno_file <- system.file("extdata", "Annotation_tables_3_BIOMASS.rds", 
                         package = "switchgrassGWAS")
anno_tables <- readRDS(anno_file)

## -----------------------------------------------------------------------------

kable(anno_tables)

