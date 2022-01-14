## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/annotations-",
  fig.width = 6,
  fig.asp = 0.618,
  out.width = "70%",
  fig.align = "center"
)
library(knitr)

## -----------------------------------------------------------------------------
library(bigsnpr)

snpfile <- system.file("extdata", "example_bigsnp.rds", package = "switchgrassGWAS")
gwasfile <- system.file("extdata", "GWAS_datatable_example_GWAS_CT.rds", 
                        package = "switchgrassGWAS")
                        
gwas <- readRDS(gwasfile)
snp <- snp_attach(snpfile)

## ---- include = FALSE---------------------------------------------------------
#library(AnnotationDbi)
#library(VariantAnnotation)
#library(switchgrassGWAS)
#txdb <- loadDb(file = "../../pvdiv-genome/Pvirgatum_516_v5.1.gene.txdb.sqlite")

anno_file <- system.file("extdata", "Annotation_tables_GWAS_CT.rds", 
                         package = "switchgrassGWAS")
anno_tables <- readRDS(anno_file)

## -----------------------------------------------------------------------------
kable(anno_tables[[2]])

