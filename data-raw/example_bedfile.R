## code to prepare the example bed/bim/fam files goes here.
library(bigsnpr)
snp <- snp_attach("~/Github/pvdiv-genome/tensite_twoyear/Pvirgatum_V5_GWAS_630g_33M_tensite_twoyear.rds")

example_set <- sample(snp$genotypes$ncol, size = 18*200) # ~200 SNPs per chromosome
example_set <- sort(example_set)
snp2 <- snp_subset(snp, ind.col = example_set)
snpex <- snp_attach(snp2)

snp_writeBed(snpex, bedfile = "inst/extdata/example.bed")
