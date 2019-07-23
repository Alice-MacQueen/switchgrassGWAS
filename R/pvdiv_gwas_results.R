#' @title Compute Adjusted False Discovery Rate for pvdiv_gwas output.
#'
#' @description Determines the false discovery rate for p-values from a
#'     dataframe. By default, this function uses a Benjamini-Hochberg
#'     correction on a column called "p.value" to determine a new column of
#'     FDR-adjusted p-values.
#'
#' @param gwas_obj A dataframe with a \code{p.value} column.
#' @param alpha The type I error threshold. Default is 0.05.
#' @param proc The procedure (from \code{multtest::mt.rawp2adjp}) used
#' to calculate the False Discovery Rate. Default is Benjamini Hochberg, "BH".
#'
#' @return A \code{tbl_df()} of the original dataframe with a new column named
#' \code{FDR_adjusted_p}.
#'
#' @importFrom dplyr arrange bind_cols select rename mutate
#' @importFrom multtest mt.rawp2adjp
#' @importFrom magrittr %>%
#' @importFrom rlang !! .data
#' @importFrom tibble as_tibble
#'
#' @examples
#' \dontrun{
#' gwas_df <- read_csv(file = "LEN.full.GWAS.csv.gz", col_names = TRUE,
#'     col_types = "ciinnnn")
#' gwas_fdr_1 <- pvdiv_fdr(gwas_obj = gwas_df, alpha = 0.1)
#' gwas_fdr_05 <- pvdiv_fdr(gwas_obj = gwas_df)
#' gwas_fdr_bonf <- pvdiv_fdr(gwas_obj = gwas_df, proc = "Bonferroni")
#' }
#'
#' @export
pvdiv_fdr <- function(gwas_obj, alpha = 0.05, proc = "BH") {
  alpha <- as.numeric(alpha)
  stopifnot('p.value' %in% names(gwas_obj))

  if(is.null(gwas_obj)) {
    gwas_obj_adj <- NULL
  } else {
    res <- mt.rawp2adjp(gwas_obj$p.value, alpha = alpha, proc = proc)
    adj_p <- res$adjp[order(res$index), ] %>%
      as_tibble()
    gwas_obj_adj <- dplyr::bind_cols(gwas_obj, adj_p) %>%
      dplyr::select(-.data$rawp) %>%
      dplyr::arrange(.data$p.value) %>%
      dplyr::rename("FDR_adjusted_p" = !! eval(proc))
  }

  return(gwas_obj_adj)
}

#' Get number of conditions
#'
#' @param m The mash result
#'
#' @importFrom ashr get_pm
get_ncond = function(m){
  return(ncol(get_pm(m)))
}


#' Return the Bayes Factor for each effect
#'
#' @param m the mash result (from joint or 1by1 analysis); must have been computed using usepointmass=TRUE
#'
#' @return if m was fitted using usepointmass=TRUE then returns a vector of
#'     the log10(bf) values for each effect. That is, the jth element
#'     lbf_j is log10(Pr(Bj | g=ghat-nonnull)/Pr(Bj | g = 0)) where gha
#'     t-nonnull is the non-null part of ghat.  Otherwise returns NULL.
#'
#' @export
get_log10bf = function(m) {
  if(is.null(m$null_loglik)){
    return(NULL)
  } else {
    return((m$alt_loglik - m$null_loglik)/log(10))
  }
}


#' Find effects that are significant in at least one condition
#'
#' @param m the mash result (from joint or 1by1 analysis)
#' @param thresh indicates the threshold below which to call signals significant
#' @param conditions which conditions to include in check (default to all)
#' @param sig_fn the significance function used to extract significance from mash object; eg could be ashr::get_lfsr or ashr::get_lfdr. (Small values must indicate significant.)
#'
#' @return a vector containing the indices of the significant effects, by order of most significant to least
#'
#' @importFrom ashr get_lfsr
#'
#' @export
get_significant_results = function(m, thresh = 0.05, conditions = NULL,
                                   sig_fn = ashr::get_lfsr) {
  if (is.null(conditions)) {
    conditions = 1:get_ncond(m)
  }
  top = apply(sig_fn(m)[,conditions,drop=FALSE],1,min) # find top effect in each condition
  sig = which(top < thresh)
  ord = order(top[sig],decreasing=FALSE)
  sig[ord]
}


#' @title Get mash marker_df
#'
#' @description Pulls the names of the SNP markers from the mash object.
#'
#' @param m An object of type mash
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tibble enframe
#' @importFrom dplyr arrange
#'
#' @export
get_marker_df <- function(m){
  marker_df <- get_significant_results(m, thresh = 1) %>%
    enframe(name = "Marker") %>%
    arrange(.data$value)

  return(marker_df)
}





