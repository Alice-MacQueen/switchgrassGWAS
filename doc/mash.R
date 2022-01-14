## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.asp = 0.618,
  out.width = "70%",
  fig.align = "center"
)

## ---- include = FALSE---------------------------------------------------------
library(switchgrassGWAS)

mashfile <- system.file("extdata", "Strong_Effects5000SNPs.rds", 
                        package = "switchgrassGWAS")
pairwisefilepath <- system.file("extdata",
                            "Pairwise_sharing_Strong_Effects_5000SNPs.rds", 
                            package = "switchgrassGWAS")
mash_output <- readRDS(mashfile)
pairwise_output <- readRDS(pairwisefilepath)

## -----------------------------------------------------------------------------
nbycond <- mash_plot_sig_by_condition(m = mash_output)
nbycond$ggobject

## ---- fig.width = 9-----------------------------------------------------------
mashhattan <- mash_plot_manhattan_by_condition(m = mash_output)
mashhattan$ggmanobject

## -----------------------------------------------------------------------------
pairwise_plot <- mash_plot_pairwise_sharing(effectRDS = pairwisefilepath,
                                            reorder = TRUE) 
pairwise_plot$gg_corr

## -----------------------------------------------------------------------------
effects <- mash_plot_effects(m = mash_output, n = 1)
effects$ggobject 

