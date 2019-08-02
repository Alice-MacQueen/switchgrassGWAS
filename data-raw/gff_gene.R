## code to prepare `gff_gene` dataset goes here
library(tidyverse)

gff_gene <- read_delim(file = file.path("C:", "Users", "ahm543", "OneDrive",
                                    "Juenger Projects", "SwitchgrassGWAS",
                                    "Package", "Public", "switchgrassGWAS",
                                    "data-raw",
                                    "Pvirgatum_516_v5.1.gene.gff3.gz"),
                   col_names = FALSE, delim = "\t", skip = 3)

gff_gene <- gff_gene %>%
  rename(seqid = X1, source = X2, type = X3, start = X4, end = X5, score = X6,
         strand = X7, phase = X8, attributes = X9)

usethis::use_data(gff_gene, compress = "gzip")
