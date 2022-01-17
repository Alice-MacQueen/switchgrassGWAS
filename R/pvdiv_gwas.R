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
#' @param snp Genomic information to include for Panicum virgatum. SNP data
#'     is available at doi:10.18738/T8/ET9UAU
#' @param covar Covariance matrix to include in the regression. You
#'     can generate these using \code{bigsnpr::snp_autoSVD()}.
#' @param ncores Number of cores to use. Default is one.
#' @param npcs Integer vector of principle components to use.
#'     Defaults to c(0:10).
#' @param saveoutput Logical. Should output be saved as a csv to the
#'     working directory?
#'
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom dplyr mutate rename case_when mutate_if
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
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
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

  LambdaGC <- as_tibble(matrix(data =
                                 c(npcs, rep(NA, (ncol(df) - 1)*length(npcs))),
                               nrow = length(npcs), ncol = ncol(df),
                               dimnames = list(npcs, colnames(df))))
  LambdaGC <- LambdaGC %>%
    dplyr::rename("NumPCs" = .data$PLANT_ID) %>%
    mutate_if(is.integer, as.numeric)

  for(i in seq_along(names(df))[-1]){

    for(k in c(1:length(npcs))){

      if(type == "linear"){

        y1 <- as_vector(df[which(!is.na(df[,i])), i])
        ind_y <- which(!is.na(df[,i]))

        if(npcs[k] == 0){

          gwaspc <- big_univLinReg(G, y.train = y1, ind.train = ind_y,
                                   ncores = ncores)
        } else {

          ind_u <- matrix(covar$u[which(!is.na(df[,i])),1:npcs[k]],
                          ncol = npcs[k])
          gwaspc <- big_univLinReg(G, y.train = y1, covar.train = ind_u,
                                   ind.train = ind_y, ncores = ncores)
        }
    } else if(type == "logistic"){
      message(paste0("For logistic models, if convergence is not reached by ",
      "the main algorithm for some SNPs, the corresponding `niter` element ",
      "is set to NA, and glm is used instead. If glm can't ",
      "converge either, those SNP estimations are set to NA."))
      y1 <- as_vector(df[which(!is.na(df[,i])), i])
      ind_y <- which(!is.na(df[,i]))
        if(npcs[k] == 0){
          gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                    ind.train = ind_y,
                                                    ncores = ncores))
        } else {
          ind_u <- matrix(covar$u[which(!is.na(df[,i])),1:npcs[k]],
                          ncol = npcs[k])
          gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                    covar.train = ind_u,
                                                    ind.train = ind_y,
                                                    ncores = ncores))
        }
    }
      ps <- predict(gwaspc, log10 = FALSE)
      LambdaGC[k,i] <- get_lambdagc(ps = ps)
      message(paste0("Finished Lambda_GC calculation for ", names(df)[i],
                     " using ", npcs[k], " PCs."))
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
    best_LambdaGC <- get_best_PC_df(df = LambdaGC)
    write_csv(best_LambdaGC, path = paste0("Best_Lambda_GC_", names(df)[2],
                                           "_to_", tail(names(df), n = 1),
                                           "_Phenotypes_", npcs[1], "_to_",
                                           tail(npcs, n = 1), "_PCs.csv"))
  }
  return(LambdaGC)
}

#' Return best number of PCs in terms of lambda_GC for Panicum virgatum.
#' Return best number of PCs in terms of lambda_GC for the CDBN.
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
#' @importFrom tidyselect all_of
#'
#' @return A dataframe containing the best lambda_GC value and number of PCs
#'     for each phenotype in the data frame.
get_best_PC_df <- function(df){
  column <- names(df)[ncol(df)]
  bestPCs <- df %>%
    filter(!! sym(column) < 1.05| !! sym(column) == min(!! sym(column))) %>%
    top_n(n = -1, wt = .data$NumPCs) %>%
    select(.data$NumPCs, all_of(column))

  if(ncol(df) > 2){
    for(i in c((ncol(df)-2):1)){
      column <- names(df)[i+1]

      bestPCs <- df %>%
        filter(!! sym(column) < 1.05 | !! sym(column) == min(!! sym(column))) %>%
        top_n(n = -1, wt = .data$NumPCs) %>%
        select(.data$NumPCs, all_of(column)) %>%
        full_join(bestPCs, by = c("NumPCs", (column)))
    }
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
#' @param snp Genomic information to include for Panicum virgatum. SNP data
#'     is available at doi:10.18738/T8/ET9UAU#'
#' @param covar Optional covariance matrix to include in the regression. You
#'     can generate these using \code{bigsnpr::snp_autoSVD()}.
#' @param ncores Number of cores to use. Default is one.
#' @param npcs Number of principle components to use. Default is 10.
#' @param saveoutput Logical. Should output be saved as a rds to the
#'     working directory?
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
                       covar = NA, ncores = 1, npcs = 10, saveoutput = FALSE){
  stopifnot(type %in% c("linear", "logistic"))
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(colnames(df)[1] != "PLANT_ID"){
    stop("First column of phenotype dataframe (df) must be 'PLANT_ID'.")
  }
  G <- snp$genotypes

  for(i in seq_along(names(df))[-1]){
    y1 <- as_vector(df[which(!is.na(df[,i])), i])
    ind_y <- which(!is.na(df[,i]))

    if(type == "linear"){
      if(!is.na(covar[1])){
        ind_u <- matrix(covar$u[which(!is.na(df[,i])),1:npcs], ncol = npcs)
        gwaspc <- big_univLinReg(G, y.train = y1, covar.train = ind_u,
                                 ind.train = ind_y, ncores = ncores)
      } else {
        gwaspc <- big_univLinReg(G, y.train = y1, ind.train = ind_y,
                                 ncores = ncores)
      }
    } else if(type == "logistic"){
      message(paste0("For logistic models, if convergence is not reached by ",
        "the main algorithm for any SNP, the corresponding `niter` element ",
        "is set to NA, and glm is used instead. If glm can't ",
        "converge either, those SNP estimations are set to NA."))
      if(!is.na(covar[1])){
        ind_u <- matrix(covar$u[which(!is.na(df[,i])),1:npcs], ncol = npcs)
        gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                  covar.train = ind_u,
                                                  ind.train = ind_y,
                                                  ncores = ncores))
      } else {
        gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                  ind.train = ind_y,
                                                  ncores = ncores))
      }
    } else {
      stop(paste0("Type of GWAS not recognized: please choose one of 'linear'",
                  " or 'logistic'"))
    }

    if(saveoutput){
      saveRDS(gwaspc, file = paste0("GWAS_object_", names(df)[i], ".rds"))
    } else {
        print("saveoutput is FALSE so GWAS object will not be saved to disk.")
      }

  }
  return(gwaspc)
}



#' Return a number rounded to some number of digits
#'
#' @description Given some x, return the number rounded to some number of
#'     digits.
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
#' @importFrom dplyr distinct
#'
#' @return A dataframe containing the 2-d binned values for x and y, and their
#'     confidence intervals.
round_xy <- function(x, y, cl = NA, cu = NA, roundby = 0.001){
  expected <- round2(x, at = roundby)
  observed <- round2(y, at = roundby)
  if(!is.na(cl[1]) & !is.na(cu[1])){
    clower <- round2(cl, at = roundby)
    cupper <- round2(cu, at = roundby)
    tp <- tibble(expected = expected, observed = observed, clower = clower,
                 cupper = cupper)
    tp <- tp %>% dplyr::distinct()
    return(tp)
  } else {
    tp <- tibble(expected = expected, observed = observed)
    tp <- tp %>% dplyr::distinct()
    return(tp)
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
#' @param effects a gwas effects FBM object created using 'pvdiv_standard_gwas'.
#'     Saved under the name "gwas_effects_{suffix}.rds" and can be loaded into
#'     R using the bigstatsr function "big_attach".
#' @param ind If effects is a FBM object, this should be the row number of the
#'     phenotype from the associated metadata for the FBM object.
#' @param ci Numeric. Size of the confidence interval, 0.95 by default.
#' @param lambdaGC Logical. Add the Genomic Control coefficient as subtitle to
#'     the plot?
#'
#' @import ggplot2
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom stats qbeta ppoints
#' @param tol Numeric. Tolerance for optional Genomic Control coefficient.
#'
#' @return A ggplot2 plot.
#'
#' @export
pvdiv_qqplot <- function(ps, effects = NULL, ind = NULL, ci = 0.95,
                       lambdaGC = FALSE, tol = 1e-8) {
  if(!is.null(effects) & !is.null(ind)){
    ind <- ind*3
    roundFBM <- function(X, ind, at) ceiling(X[, ind] / at) * at
    observed <- big_apply(effects, ind = ind, a.FUN = roundFBM, at = 0.01,
                          a.combine = 'plus', ncores = 1)
    ps <- 10^-observed
  }
  ps <- ps[which(!is.na(ps))]
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
    ylab(log10Po) +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.line.x = element_line(size = 0.35, colour = 'grey50'),
          axis.line.y = element_line(size = 0.35, colour = 'grey50'),
          axis.ticks = element_line(size = 0.25, colour = 'grey50'),
          legend.justification = c(1, 0.75), legend.position = c(1, 0.9),
          legend.key.size = unit(0.35, 'cm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 9),
          legend.text.align = 0, legend.background = element_blank(),
          plot.subtitle = element_text(size = 10, vjust = 0),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0.5, size = 10 ,vjust = 0),
          strip.placement = 'outside', panel.spacing.x = unit(-0.4, 'cm'))

  if (lambdaGC) {
    lamGC <- get_lambdagc(ps = ps, tol = tol)
    expr <- substitute(expression(lambda[GC] == l), list(l = lamGC))
    p1 + labs(subtitle = eval(expr))
  } else {
    p1
  }
}

#' Create a Manhattan plot with ggplot2.
#'
#' @description Create a Manhattan plot using ggplot2 on either RDS or FBM
#'     object GWAS results.
#'
#' @param effects Either a gwas effects RDS or FBM object created using
#'     'pvdiv_standard_gwas' (with savetype = "rds" or savetype = "fbm"). If a
#'     fbm, this file is saved under the name "gwas_effects_{suffix}.rds" and
#'     can be loaded into R using the bigstatsr function "big_attach".
#' @param ind If effects is a FBM object, this should be the row number of the
#'     phenotype from the associated metadata for the FBM object.
#' @param snp If effects is a FBM object, you must also supply a "bigSNP" object;
#'     load into R with \code{bigsnpr::snp_attach()}.
#' @param thr Numeric. Significance threshold plotted as a horizontal line.
#'     Default is Bonferroni.
#' @param ncores Integer. Number of cores to use for parallelization.
#'
#' @import ggplot2
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom dplyr distinct mutate
#' @importFrom stats qbeta ppoints
#'
#' @return A ggplot2 plot.
#'
#' @export
pvdiv_manhattan <- function(effects, ind = NULL, snp = NULL, thr = NULL,
                            ncores = 1){
  if(!is.null(ind) & !is.null(snp)){
    if (attr(snp, "class") != "bigSNP") {
      stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
    }
    ind <- ind*3
    roundFBM <- function(X, ind, at) ceiling(X[, ind] / at) * at
    observed <- big_apply(effects, ind = ind, a.FUN = roundFBM, at = 0.01,
                          a.combine = 'plus', ncores = ncores)

    plot_data <- tibble(CHR = snp$map$chromosome, POS = snp$map$physical.pos,
                        observed = observed)
  } else if(is.null(ind) & !is.null(snp)) {
    stop("must specify both ind and snp if effects is a fbm")
  } else if(!is.null(ind) & is.null(snp)) {
    stop("must specify both ind and snp if effects is a fbm")
  } else {
    plot_data <- tibble(CHR = effects$CHR, POS = effects$POS,
                        observed = effects$log10p) %>%
      mutate(observed = round2(.data$observed, at = 0.01))
  }
  if(is.null(thr)) {
    thr <- -log10(.05/nrow(plot_data))
  }
  # If more than half a million markers, round data slightly to reduce the
  # number of markers plotted.
  if (nrow(plot_data) >= 500000) {
    plot_data <- plot_data %>%
      mutate(POS = round2(.data$POS, at = 250000))
  }
  plot_data <- plot_data %>% dplyr::distinct()

  nchr <- length(unique(plot_data$CHR))
  log10P <- expression(paste("-log"[10], plain("("), italic(p-value),
                             plain(")")))
  p1 <- plot_data %>%
    ggplot(aes(x = .data$POS, y = .data$observed)) +
    geom_point(aes(color = .data$CHR, fill = .data$CHR)) +
    geom_hline(yintercept = thr, color = "black", linetype = 2,
               size = 1) +
    facet_wrap(~ .data$CHR, nrow = 1, scales = "free_x",
               strip.position = "bottom") +
    scale_color_manual(values = rep(c("#1B0C42FF", "#48347dFF",
                                      "#95919eFF"), ceiling(nchr/3)),
                       guide = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill=NA),
          legend.position = "none",
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.line.x = element_line(size = 0.35, colour = 'grey50'),
          axis.line.y = element_line(size = 0.35, colour = 'grey50'),
          axis.ticks = element_line(size = 0.25, colour = 'grey50'),
          legend.justification = c(1, 0.75),
          legend.key.size = unit(0.35, 'cm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 9),
          legend.text.align = 0, legend.background = element_blank(),
          plot.subtitle = element_text(size = 10, vjust = 0),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0.5, size = 10 ,vjust = 0),
          strip.placement = 'outside', panel.spacing.x = unit(-0.1, 'cm')) +
    labs(x = "Chromosome", y = log10P) +
    scale_x_continuous(expand = c(0.15, 0.15))
  return(p1)
}


#' Find lambda_GC value for non-NA p-values
#'
#' @description Finds the lambda GC value for some vector of p-values.
#'
#' @param ps Numeric vector of p-values. Can have NA's.
#' @param tol Numeric. Tolerance for optional Genomic Control coefficient.
#'
#' @importFrom stats median uniroot
#'
#' @return A lambda GC value (some positive number, ideally ~1)
#'
#' @export
get_lambdagc <- function(ps, tol = 1e-8){
  ps <- ps[which(!is.na(ps))]
  xtr <- log10(ps)
  MEDIAN <- log10(0.5)
  f.opt <- function(x) (x - MEDIAN)
  xtr_p <- median(xtr) / uniroot(f.opt, interval = range(xtr),
                                 check.conv = TRUE,
                                 tol = tol)$root
  lamGC <- signif(xtr_p)
  return(lamGC)
}


#' Adjusted p-values for simple multiple testing procedures
#'
#' @description This function computes adjusted p-values for simple multiple
#'     testing procedures from a vector of raw (unadjusted) p-values. The
#'     procedures include the Bonferroni, Holm (1979), Hochberg (1988), and
#'     Sidak procedures for strong control of the family-wise Type I error
#'     rate (FWER), and the Benjamini & Hochberg (1995) and Benjamini &
#'     Yekutieli (2001) procedures for (strong) control of the false discovery
#'     rate (FDR). The less conservative adaptive Benjamini & Hochberg (2000)
#'     and two-stage Benjamini & Hochberg (2006) FDR-controlling procedures are
#'     also included. This function is taken from the multtest package. It is
#'     the only function used from this package and is added to this package
#'     wholesale to reduce user installation burden.
#'
#' @usage mt.rawp2adjp(rawp, proc=c("Bonferroni", "Holm", "Hochberg", "SidakSS",
#'     "SidakSD", "BH", "BY","ABH","TSBH"), alpha = 0.05, na.rm = FALSE)
#'
#' @param rawp A vector of raw (unadjusted) p-values for each hypothesis under
#'     consideration. These could be nominal p-values, for example, from
#'     t-tables, or permutation p-values as given in mt.maxT and mt.minP. If
#'     the mt.maxT or mt.minP functions are used, raw p-values should be given
#'     in the original data order, ordered by the index of that data.
#' @param proc A vector of character strings containing the names of the
#'     multiple testing procedures for which adjusted p-values are to be
#'     computed. This vector should include any of the following: "Bonferroni",
#'     "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH", "TSBH".
#' @param alpha A nominal type I error rate, or a vector of error rates, used
#'     for estimating the number of true null hypotheses in the two-stage
#'     Benjamini & Hochberg procedure ("TSBH"). Default is 0.05.
#' @param na.rm An option for handling NA values in a list of raw p-values.
#'     If FALSE, the number of hypotheses considered is the length of the vector
#'      of raw p-values. Otherwise, if TRUE, the number of hypotheses is the
#'      number of raw p-values which were not NAs.
#'
#' @author Sandrine Dudoit, http://www.stat.berkeley.edu/~sandrine,
#' @author Yongchao Ge, yongchao.ge@mssm.edu,
#' @author Houston Gilbert, http://www.stat.berkeley.edu/~houston.
#'
#' @return A list with components: adjp, index, h0.ABH, h0.TSBH. See multtest
#'    package on Bioconductor for details.
mt.rawp2adjp <- function (rawp, proc = c("Bonferroni", "Holm", "Hochberg",
                         "SidakSS", "SidakSD", "BH", "BY",
                         "ABH", "TSBH"), alpha = 0.05, na.rm = FALSE)
{
  m <- length(rawp)
  if (na.rm) {
    mgood <- sum(!is.na(rawp))
  }
  else {
    mgood <- m
  }
  n <- length(proc)
  a <- length(alpha)
  index <- order(rawp)
  h0.ABH <- NULL
  h0.TSBH <- NULL
  spval <- rawp[index]
  adjp <- matrix(0, m, n + 1)
  dimnames(adjp) <- list(NULL, c("rawp", proc))
  adjp[, 1] <- spval
  if (is.element("TSBH", proc)) {
    TS.spot <- which(proc == "TSBH")
    TSBHs <- paste("TSBH", alpha, sep = "_")
    newprocs <- append(proc, TSBHs, after = TS.spot)
    newprocs <- newprocs[newprocs != "TSBH"]
    adjp <- matrix(0, m, n + a)
    dimnames(adjp) <- list(NULL, c("rawp", newprocs))
    adjp[, 1] <- spval
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood/i) * spval[i],
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i]))
        tmp[i] <- NA
    }
    h0.TSBH <- rep(0, length(alpha))
    names(h0.TSBH) <- paste("h0.TSBH", alpha, sep = "_")
    for (i in 1:length(alpha)) {
      h0.TSBH[i] <- mgood - sum(tmp < alpha[i]/(1 + alpha[i]),
                                na.rm = TRUE)
      adjp[, TS.spot + i] <- tmp * h0.TSBH[i]/mgood
    }
  }
  if (is.element("Bonferroni", proc)) {
    tmp <- mgood * spval
    tmp[tmp > 1] <- 1
    adjp[, "Bonferroni"] <- tmp
  }
  if (is.element("Holm", proc)) {
    tmp <- spval
    tmp[1] <- min(mgood * spval[1], 1)
    for (i in 2:m) tmp[i] <- max(tmp[i - 1], min((mgood -
                                                    i + 1) * spval[i], 1))
    adjp[, "Holm"] <- tmp
  }
  if (is.element("Hochberg", proc)) {
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood - i + 1) * spval[i],
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i]))
        tmp[i] <- NA
    }
    adjp[, "Hochberg"] <- tmp
  }
  if (is.element("SidakSS", proc))
    adjp[, "SidakSS"] <- 1 - (1 - spval)^mgood
  if (is.element("SidakSD", proc)) {
    tmp <- spval
    tmp[1] <- 1 - (1 - spval[1])^mgood
    for (i in 2:m) tmp[i] <- max(tmp[i - 1], 1 - (1 - spval[i])^(mgood -
                                                                   i + 1))
    adjp[, "SidakSD"] <- tmp
  }
  if (is.element("BH", proc)) {
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood/i) * spval[i],
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i]))
        tmp[i] <- NA
    }
    adjp[, "BH"] <- tmp
  }
  if (is.element("BY", proc)) {
    tmp <- spval
    a <- sum(1/(1:mgood))
    tmp[m] <- min(a * spval[m], 1)
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood * a/i) * spval[i],
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i]))
        tmp[i] <- NA
    }
    adjp[, "BY"] <- tmp
  }
  if (is.element("ABH", proc)) {
    tmp <- spval
    h0.m <- rep(0, mgood)
    for (k in 1:mgood) {
      h0.m[k] <- (mgood + 1 - k)/(1 - spval[k])
    }
    grab <- min(which(diff(h0.m, na.rm = TRUE) > 0), na.rm = TRUE)
    h0.ABH <- ceiling(min(h0.m[grab], mgood))
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood/i) * spval[i],
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i]))
        tmp[i] <- NA
    }
    adjp[, "ABH"] <- tmp * h0.ABH/mgood
  }
  list(adjp = adjp, index = index, h0.ABH = h0.ABH[1], h0.TSBH = h0.TSBH[1:length(alpha)])
}


#' Set an advanced argument
#'
#' From the \code{...} argument used in your function, find if a specific
#' argument was included and extract its value.
#'
#' @description This function is taken from the lcolladotor/dots package. It is
#'     the only function used from this package and is added to this package
#'     wholesale to reduce user installation burden. Please use the original
#'     package from Github if you use this in your own work.
#'
#' @param name Name of the advanced argument to look for in \code{...}
#' @param value The default value of the advanged argument. If this advanced
#' argument is used in several of your functions, we recommend using
#' \code{getOption('value')} and explaining this option in your package
#' vignette geared towards experienced users.
#' @param ... Advanced arguments. See \link[methods]{dotsMethods}.
#'
#' @details
#' Note that you can make dots() even more powerful by using \link{getOption}
#' to define \code{value}. This is particularly useful if you use the
#' same advanced argument in several functions.
#'
#' @export
#' @seealso \link[methods]{dotsMethods}
#' @aliases advanced...
#' @author L. Collado-Torres
#'
#' @examples
#'
#' ## Simple example that calculates the max between 'x' and 'y' with a
#' ## specified minimum value to return.
#'
#' minMax <- function(x, y, ...) {
#'     minValue <- dots('minValue', 0, ...)
#'     res <- max(x, y, minValue)
#'     return(res)
#' }
#' minMax(1:2, 3:4)
#' minMax(1:2, 3:4, minValue = 5)
#'
#'
#' ## Arguably these examples are simple, but the idea is that dots()
#' ## can simplify very long function calls where some parameters will be used
#' ## by a minority of the users.
#'
dots <- function(name, value, ...) {
  args <- list(...)
  if(!name %in% names(args)) {
    ## Default value
    return(value)
  } else {
    ## If the argument was defined in the ... part, return it
    return(args[[name]])
  }
}
