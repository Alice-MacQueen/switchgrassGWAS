## code to prepare `anno_info` dataset goes here
library(tidyverse)
all_anno_info <-
  read_delim(file = file.path("~", "Github", "pvdiv-genome",
                              "Pvirgatum_516_v5.1.annotation_info.txt"),
             col_names = TRUE, delim = "\t")
anno_info <- tbl_df(all_anno_info) %>%
  distinct(locusName, .keep_all = TRUE)

usethis::use_data(anno_info)
