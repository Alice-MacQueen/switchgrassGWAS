# code to create phenotypes.rda goes here

# These phenotypes were merged and then carefully checked for errors and
# outliers in the `pvdiv-phenotypes` git repository.

# Thereafter, specific phenotypes were published as part of the switchgrass
# genome paper - seven environmental phenotypes, overwinter survival at three
# northern sites, and biomass in 2019 at the three core common gardens.

# Additionally, we include GWAS_CT, the number of times that genotype was
# clonally propagated and planted at these common gardens.

library(tidyverse)

phe_ex <- readRDS("~/Github/pvdiv-phenotypes/data/pre_replant/Phenotypes_all_pre_2019_replant.rds")
exampleGWAS <- phe_ex %>%
  filter(!(PLANT_ID %in% c("AP13", "UNK"))) %>%
  group_by(PLANT_ID, PLOT_GL) %>%
  tally() %>% tally(name = "GWAS_CT")

envGWAS <- readRDS("~/Github/pvdiv-fitness-2019/data/Seven_climate_gwas.rds")
fitnessGWAS <- readRDS("~/Github/pvdiv-fitness-2019/data/Phenotypes_fitness_linear.rds") %>%
  dplyr::select(PLANT_ID, FRAC_SRV_THREE, CLMB_BIOMASS, KBSM_BIOMASS, PKLE_BIOMASS)

pvdiv_phenotypes <- exampleGWAS %>%
  full_join(envGWAS) %>%
  full_join(fitnessGWAS)

usethis::use_data(pvdiv_phenotypes, compress = "gzip", overwrite = TRUE)
