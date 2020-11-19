## code to prepare `pvdiv_metadata` dataset goes here
##
pvdiv_metadata <- readRDS("~/Github/pvdiv-phenotypes/data/metadata.rds")

usethis::use_data(pvdiv_metadata, overwrite = TRUE)
