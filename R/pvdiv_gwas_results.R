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
get_log10bf = function(m) {
  if(is.null(m$null_loglik)){
    return(NULL)
  } else {
    return((m$alt_loglik - m$null_loglik)/log(10))
  }
}


#' From a mash result, get effects that are significant in at least one condition
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
get_marker_df <- function(m){
  marker_df <- get_significant_results(m, thresh = 1) %>%
    enframe(name = "Marker") %>%
    arrange(.data$value)

  return(marker_df)
}


#' @title Make a SNP dataframe subset from annotation table or a genomic range.
#'
#' @description Given a dataframe created using pvdiv_table_topsnps(), this
#'     function creates a dataframe of SNP calls for the subset of SNPs from
#'     this annotation table.
#'
#' @note This function is a wrapper around bigsnpr functions to subset its SNP
#'     file format that may be useful if you have a small interval to look at or
#'     a small number of SNPs from an annotation table.
#'
#' @param snp A `FBM.code256` object. Genomic information for Panicum virgatum.
#'    Contact tjuenger <at> utexas <dot> edu to obtain this information
#'    pre-publication.
#' @param type One of "anno" or "range", depending on if you are using an
#'    annotation dataframe from pvdiv_table_topsnps() or a genomic interval.
#' @param anno_df One dataframe of annotations from pvdiv_table_topsnps(). This
#'    dataframe needs to contain the columns CHR and region_start. It's
#'    recommended that you set rangevector = 0 in pvdiv_table_topsnps() to get
#'    the SNP itself using this function.
#' @param chr Character string. The chromsome (e.g., "Chr01K") to get SNPs from.
#' @param pos1 Integer. The low position to start getting SNPs from.
#' @param pos2 Integer. The high position to stop getting SNPs from.
#'
#' @return A `FBM.code256` object for the subset of SNPs in the annotation data
#'    frame or in the genomic interval.
#'
#' @importFrom bigstatsr rows_along cols_along
#' @import bigsnpr
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate between
#' @importFrom tibble as_tibble
#'
#' @export
pvdiv_bigsnp_subset <- function(snp, type = c("anno", "range"), anno_df, chr,
                                pos1, pos2){
  if(type == "anno"){
    stopifnot("CHR" %in% names(anno_df) & "POS" %in% names(anno_df))
    anno_df <- anno_df %>%
      mutate(marker_ID = paste(.data$CHR, .data$POS, sep = "_"))
    pos_subset <- which(snp$map$marker.ID %in% anno_df$marker_ID)
  } else if(type == "range"){
    stopifnot(is.character(chr) & is.numeric(pos1) & is.numeric(pos2) & pos2 > pos1)
    stopifnot(chr %in% c("Chr01K", "Chr01N", "Chr02K", "Chr02N", "Chr03K",
                         "Chr03N", "Chr04K", "Chr04N", "Chr05K", "Chr05N",
                         "Chr06K", "Chr06N", "Chr07K", "Chr07N", "Chr08K",
                         "Chr08N", "Chr09K", "Chr09N"))
    pos_df <- which(between(snp$map$physical.pos, pos1, pos2))
    chr_df <- which(snp$map$chromosome %in% chr)
    pos_subset <- pos_df[which(pos_df %in% chr_df)]
  } else {
    stop("Subset type must be one of 'anno' or 'range'.")
  }
  subset_name <- subset(snp, ind.col = pos_subset)
  subset <- snp_attach(subset_name)
  return(subset)
}


#' @title Convert A `FBM.code256` SNP subset to a dataframe.
#'
#' @description This function creates a dataframe of SNP calls for a
#'    `FBM.code256` object. We recommend that this be a dataframe for a small
#'    region.
#'
#' @note bigsnpr already has functions to subset its SNP file format and return
#'     a bed file and to return its own SNP file format. This function is useful
#'     if you instead want a small(ish) dataframe of SNPs that can be
#'     manipulated using data frame tools in R.
#'
#' @param snp A `FBM.code256` object. Genomic information for Panicum virgatum.
#'    Contact tjuenger <at> utexas <dot> edu to obtain this information
#'    pre-publication.
#'
#' @return A \code{tbl_df()} of the SNP calls for all individuals for the
#'    subset of SNPs on that chromosome between the two positions specified.
#'
#' @importFrom bigstatsr rows_along cols_along
#' @import bigsnpr
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom tibble as_tibble
#'
#' @export
pvdiv_bigsnp2tibble <- function(snp){
  subset_df <- snp$genotypes[rows_along(snp$genotypes),
                             cols_along(snp$genotypes)]
  colnames(subset_df) <- snp$map$marker.ID
  subset_tibble <- data.frame(PLANT_ID = snp$fam$sample.ID, subset_df) %>%
    as_tibble() %>%
    mutate(PLANT_ID = as.character(.data$PLANT_ID))
  return(subset_tibble)
}
