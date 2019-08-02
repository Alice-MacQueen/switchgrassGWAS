## code to prepare `gff_metadata` dataset goes here
library(readr)

gff_metadata <- read_delim(file = file.path("C:", "Users", "ahm543",
                                            "OneDrive", "Juenger Projects",
                                            "SwitchgrassGWAS", "Package",
                                            "Public", "switchgrassGWAS",
                                            "data-raw",
                                            "Pvirgatum_516_v5.1.gene.gff3.gz"),
                      col_names = FALSE, n_max = 3, comment = "@", delim = ";")

usethis::use_data(gff_metadata)
