#' Wrapper for bigsnpr for GWAS on Panicum virgatum.
#'
#' @description Given a dataframe of phenotypes associated with PLANT_IDs, this
#'     function is a wrapper around bigsnpr functions to conduct linear or
#'     logistic regression on Panicum virgatum. The main advantages of this
#'     function over just using the bigsnpr functions is that it automatically
#'     removes individual genotypes with missing phenotypic data, that it
#'     converts switchgrass chromosome names to the format bigsnpr requires,
#'     and that it can run GWAS on multiple phenotypes sequentially.
#'
#' @param df Dataframe of phenotypes where the first column is PLANT_ID.
#' @param type Character string. Type of univarate regression to run for GWAS.
#'     Options are "linear" or "logistic".
#' @param snp Genomic information to include for Panicum virgatum. Contact
#'     tjuenger <at> utexas <dot> edu to obtain this information pre-publication.
#' @param covar Optional covariance matrix to include in the regression. You
#'     can generate these using \code{bigsnpr::snp_autoSVD()}.
#' @param ncores Number of cores to use. Default is one.
#'
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom dplyr mutate rename case_when
#' @importFrom purrr as_vector
#' @importFrom tibble as_tibble enframe
#' @importFrom rlang .data
#'
#' @return The gwas results for the last phenotype in the dataframe. That
#'     phenotype, as well as the remaining phenotypes, are saved as RDS objects
#'     in the working directory.
#'
#' @export
pvdiv_gwas <- function(df, type = c("linear", "logistic"), snp,
                       covar = NA, ncores = 1){

  G <- snp$genotypes
  CHR <- snp$map$chromosome
  POS <- snp$map$physical.pos
  #NCORES <- nb_cores()

  # Make the switchgrass chromosome names numeric, which bigsnpr requires.
  # Scaffolds that remain will be called 19, but note for some analyses that
  # they need to be ordered (so two scaffolds can't have the same number)
  CHRN <- enframe(CHR, name = NULL) %>%
    dplyr::rename(CHR = .data$value) %>%
    mutate(CHRN = case_when(.data$CHR == "Chr01K" ~ 1,
                            .data$CHR == "Chr01N" ~ 2,
                            .data$CHR == "Chr02K" ~ 3,
                            .data$CHR == "Chr02N" ~ 4,
                            .data$CHR == "Chr03K" ~ 5,
                            .data$CHR == "Chr03N" ~ 6,
                            .data$CHR == "Chr04K" ~ 7,
                            .data$CHR == "Chr04N" ~ 8,
                            .data$CHR == "Chr05K" ~ 9,
                            .data$CHR == "Chr05N" ~ 10,
                            .data$CHR == "Chr06K" ~ 11,
                            .data$CHR == "Chr06N" ~ 12,
                            .data$CHR == "Chr07K" ~ 13,
                            .data$CHR == "Chr07N" ~ 14,
                            .data$CHR == "Chr08K" ~ 15,
                            .data$CHR == "Chr08N" ~ 16,
                            .data$CHR == "Chr09K" ~ 17,
                            .data$CHR == "Chr09N" ~ 18,
                            TRUE ~ 19
    ))

  for(i in seq_along(names(df))[-1]){
    y1 <- as_vector(df[!is.na(df[,i]), i])
    ind_y <- which(!is.na(df[,i]))
    if(!is.na(covar[1])){
      ind_u <- covar$u[!is.na(df[,i]),] # 4 PC's for phenotyped individuals.
      gwaspc <- big_univLinReg(G, y.train = y1, covar.train = ind_u,
                               ind.train = ind_y, ncores = ncores)
    } else {
      gwaspc <- big_univLinReg(G, y.train = y1, ind.train = ind_y,
                               ncores = ncores)
    }

    saveRDS(gwaspc, file = paste0("GWAS_object_", names(df)[i], ".rds"))

  }
  return(gwaspc)
}
