#' Return lambda_GC for different numbers of PCs for GWAS on Panicum virgatum.
#'
#' @description Given a dataframe of phenotypes associated with PLANT_IDs and
#'     output from a PCA to control for population structure, this function will
#'     return a .csv file of the lambda_GC values for the GWAS upon inclusion
#'     of different numbers of PCs. This allows the user to choose a number of
#'     PCs that returns a lambda_GC close to 1, and thus ensure that they have
#'     done adequate correction for population structure.
#'
#' @param df Dataframe of phenotypes where the first column is PLANT_ID and each
#'     PLANT_ID occurs only once in the dataframe.
#' @param type Character string. Type of univarate regression to run for GWAS.
#'     Options are "linear" or "logistic".
#' @param snp Genomic information to include for Panicum virgatum. Contact
#'     tjuenger <at> utexas <dot> edu to obtain this information pre-publication.
#' @param covar Covariance matrix to include in the regression. You
#'     can generate these using \code{bigsnpr::snp_autoSVD()}.
#' @param ncores Number of cores to use. Default is one.
#' @param npcs Integer vector of principle components to use.
#'     Defaults to c(0:10).
#' @param saveoutput Logical. Should output be saved as a csv to the working directory?
#'
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom dplyr mutate rename case_when
#' @importFrom purrr as_vector
#' @importFrom tibble as_tibble enframe
#' @importFrom rlang .data
#' @importFrom readr write_csv
#' @importFrom utils tail
#'
#' @return A dataframe containing the lambda_GC values for each number of PCs
#'     specified. This is also saved as a .csv file in the working directory.
#'
#' @export
pvdiv_lambda_GC <- function(df, type = c("linear", "logistic"), snp,
                       covar = NA, ncores = 1, npcs = c(0:10),
                       saveoutput = FALSE){
  if(colnames(df)[1] != "PLANT_ID"){
    stop("First column of phenotype dataframe (df) must be 'PLANT_ID'.")
    }
  if(length(covar) == 1){
    stop(paste0("Need to specify covariance matrix (covar) and a vector of",
                " PC #'s to test (npcs)."))
  }
  if(saveoutput == FALSE){
    message("saveoutput is FALSE, so lambda_GC values won't be saved to a csv.")
  }

  G <- snp$genotypes
  CHR <- snp$map$chromosome
  POS <- snp$map$physical.pos

  # Make the switchgrass chromosome names numeric, which bigsnpr requires.
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

  LambdaGC <- as_tibble(matrix(data =
                                 c(npcs, rep(NA, (ncol(df) - 1)*length(npcs))),
                               nrow = length(npcs), ncol = ncol(df),
                               dimnames = list(npcs, colnames(df))))
  LambdaGC <- LambdaGC %>%
    dplyr::rename("NumPCs" = .data$PLANT_ID)

  for(i in seq_along(names(df))[-1]){
    y1 <- as_vector(df[!is.na(df[,i]), i])
    ind_y <- which(!is.na(df[,i]))
    for(k in c(1:length(npcs))){
      if(npcs[k] == 0){
        gwaspc <- big_univLinReg(G, y.train = y1, ind.train = ind_y,
                                 ncores = ncores)
      } else {
        ind_u <- matrix(covar$u[which(!is.na(df[,i])),1:npcs[k]], ncol = npcs[k])
        gwaspc <- big_univLinReg(G, y.train = y1, covar.train = ind_u,
                                 ind.train = ind_y, ncores = ncores)
      }
      gwas2 <- gwaspc[which(!is.na(gwaspc$score)),]
      LambdaGC[k,i] <- bigsnpr:::getLambdaGC(gwas = gwas2)
      message(paste0("Finished Lambda_GC calculation for ", names(df)[i], " using ", npcs[k], " PCs."))
    }
    if(saveoutput == TRUE){
      write_csv(LambdaGC, path = paste0("Lambda_GC_", names(df)[i], ".csv"))
    }
    message(paste0("Finished phenotype ", i-1, ": ", names(df)[i]))
  }
  if(saveoutput == TRUE){
    write_csv(LambdaGC, path = paste0("Lambda_GC_", names(df)[2], "_to_",
                                      tail(names(df), n = 1), "_Phenotypes_",
                                      npcs[1], "_to_", tail(npcs, n = 1),
                                      "_PCs.csv"))
    best_LambdaGC <- pvdiv_best_PC_df(df = LambdaGC)
    write_csv(best_LambdaGC, path = paste0("Best_Lambda_GC_", names(df)[2], "_to_",
                                           tail(names(df), n = 1), "_Phenotypes_",
                                           npcs[1], "_to_", tail(npcs, n = 1),
                                           "_PCs.csv"))
  }
  return(LambdaGC)
}

#' Return best number of PCs in terms of lambda_GC for Panicum virgatum.
#'
#' @description Given a dataframe created using pvdiv_lambda_GC, this function
#'     returns the first lambda_GC less than 1.05, or the smallest lambda_GC,
#'     for each column in the dataframe.
#'
#' @param df Dataframe of phenotypes where the first column is NumPCs and
#'     subsequent column contains lambda_GC values for some phenotype.
#'
#' @importFrom dplyr filter top_n select full_join arrange
#' @importFrom tidyr gather
#' @importFrom rlang .data sym !!
#'
#' @return A dataframe containing the best lambda_GC value and number of PCs
#'     for each phenotype in the data frame.
pvdiv_best_PC_df <- function(df){
  column <- names(df)[ncol(df)]
  bestPCs <- df %>%
    filter(!! sym(column) < 1.05) %>%
    top_n(n = -1, wt = .data$NumPCs) %>%
    select(.data$NumPCs, column)

  for(i in c((ncol(df)-2):1)){
    column <- names(df)[i+1]

    bestPCs <- df %>%
      filter(!! sym(column) < 1.05 | !! sym(column) == min(!! sym(column))) %>%
      top_n(n = -1, wt = .data$NumPCs) %>%
      select(.data$NumPCs, column) %>%
      full_join(bestPCs, by = "NumPCs")
  }

  bestPCdf <- bestPCs %>%
    arrange(.data$NumPCs) %>%
    gather(key = "trait", value = "lambda_GC", 2:ncol(bestPCs)) %>%
    filter(!is.na(.data$lambda_GC))

  return(bestPCdf)
}

#' Return best number of PCs in terms of lambda_GC following Cattrell's rule.
#'
#' @description Given a dataframe created using pvdiv_lambda_GC, this function
#'     returns the lambda_GC that is closest to 1 for each column in the
#'     dataframe.
#'
#' @param df Dataframe of phenotypes where the first column is NumPCs and
#'     subsequent column contains lambda_GC values for some phenotype.
#'
#' @importFrom dplyr filter top_n select full_join arrange
#' @importFrom tidyr gather
#' @importFrom rlang .data sym !!
#'
#' @return A dataframe containing the best lambda_GC value and number of PCs
#'     for each phenotype in the data frame.
asv_best_PC_df <- function(df){
  column <- names(df)[ncol(df)]
  bestPCs <- df %>%
    filter(abs(!! sym(column)-1) == min(abs(!! sym(column)-1))) %>%
    top_n(n = -1, wt = .data$NumPCs) %>%
    select(.data$NumPCs, column)

  for(i in c((ncol(df)-2):1)){
    column <- names(df)[i+1]

    bestPCs <- df %>%
      filter(abs(!! sym(column)-1) == min(abs(!! sym(column)-1))) %>%
      top_n(n = -1, wt = .data$NumPCs) %>%
      select(.data$NumPCs, column) %>%
      full_join(bestPCs, by = "NumPCs")
  }

  bestPCdf <- bestPCs %>%
    arrange(.data$NumPCs) %>%
    gather(key = "trait", value = "lambda_GC", 2:ncol(bestPCs)) %>%
    filter(!is.na(.data$lambda_GC))

  return(bestPCdf)
}

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
#' @param npcs Number of principle components to use. Default is 10.
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
                       covar = NA, ncores = 1, npcs = 10){

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
      ind_u <- matrix(covar$u[which(!is.na(df[,i])),1:npcs], ncol = npcs)
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



#' Return a number rounded to some number of digits
#'
#' @description Given some x, return the number rounded to some number of digits.
#'
#' @param x A number or vector of numbers
#' @param at Numeric. Rounding factor or size of the bin to round to.
#'
#' @return A number or vector of numbers
round2 <- function(x, at) ceiling(x / at) * at

#' Return a dataframe binned into 2-d bins by some x and y.
#'
#' @description Given a dataframe of x and y values (with some optional
#'     confidence intervals surrounding the y values), return only the unique
#'     values of x and y in some set of 2-d bins.
#'
#' @param x Numeric vector. The first vector for binning.
#' @param y Numeric vector. the second vector for binning
#' @param cl Numeric vector. Optional confidence interval for the y vector,
#'     lower bound.
#' @param cu Numeric vector. Optional confidence interval for the y vector,
#'     upper bound.
#' @param roundby Numeric. The amount to round the x and y vectors by for 2d
#'     binning.
#'
#' @return A dataframe containing the 2-d binned values for x and y, and their
#'     confidence intervals.
round_xy <- function(x, y, cl = NA, cu = NA, roundby = 0.001){
  expected <- round2(x, at = roundby)
  observed <- round2(y, at = roundby)
  if(!is.na(cl[1]) & !is.na(cu[1])){
    clower <- round2(cl, at = roundby)
    cupper <- round2(cu, at = roundby)
    tp <- cbind(expected, observed, clower, cupper)
    return(tp[!duplicated(tp),])
  } else {
    tp <- cbind(expected, observed)
    return(tp[!duplicated(tp),])
  }
}

#' Create a quantile-quantile plot with ggplot2.
#'
#' @description Assumptions for this quantile quantile plot:
#'     Expected P values are uniformly distributed.
#'     Confidence intervals assume independence between tests.
#'     We expect deviations past the confidence intervals if the tests are
#'     not independent.
#'     For example, in a genome-wide association study, the genotype at any
#'     position is correlated to nearby positions. Tests of nearby genotypes
#'     will result in similar test statistics.
#'
#' @param ps Numeric vector of p-values.
#' @param ci Numeric. Size of the confidence interval, 0.95 by default.
#' @param lambdaGC Logical. Add the Genomic Control coefficient as subtitle to
#'     the plot?
#' @param tol Numeric. Tolerance for optional Genomic Control coefficient.
#'
#' @import ggplot2
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom stats qbeta median uniroot ppoints
#'
#' @return A ggplot2 plot.
#'
#' @export
pvdiv_qqplot <- function(ps, ci = 0.95, lambdaGC = FALSE, tol = 1e-8) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  df_round <- round_xy(df$expected, df$observed, cl = df$clower, cu = df$cupper)
  log10Pe <- expression(paste("Expected -log"[10], plain("("), italic(p-value),
                              plain(")")))
  log10Po <- expression(paste("Observed -log"[10], plain("("), italic(p-value),
                              plain(")")))
  p1 <- ggplot(as_tibble(df_round)) +
    geom_point(aes(.data$expected, .data$observed), shape = 1, size = 1) +
    geom_abline(intercept = 0, slope = 1, size = 1.5, color = "red") +
    geom_line(aes(.data$expected, .data$cupper), linetype = 2) +
    geom_line(aes(.data$expected, .data$clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po)

  if (lambdaGC) {
    xtr <- log10(ps)
    MEDIAN <- log10(0.5)
    f.opt <- function(x) (x - MEDIAN)
    xtr_p <- stats::median(xtr) / stats::uniroot(f.opt, interval = range(xtr),
                                                 check.conv = TRUE,
                                                 tol = tol)$root
    lamGC <- signif(xtr_p)
    expr <- substitute(expression(lambda[GC] == l), list(l = lamGC))
    p1 + labs(subtitle = eval(expr))
  } else {
    p1
  }
}

