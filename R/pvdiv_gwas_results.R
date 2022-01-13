#' @title Compute Adjusted False Discovery Rate for pvdiv_gwas output.
#'
#' @description Determines the false discovery rate for p-values from a
#'     dataframe. By default, this function uses a Benjamini-Hochberg
#'     correction on a column called "p.value" to determine a new column of
#'     FDR-adjusted p-values.
#'
#' @param gwas_obj A dataframe with a \code{p.value} column.
#' @param alpha The type I error threshold. Default is 0.05.
#' @param proc The procedure (used wholesale from \code{multtest::mt.rawp2adjp})
#'     used to calculate the False Discovery Rate. Default is Benjamini
#'     Hochberg, "BH".
#'
#' @return A \code{tbl_df()} of the original dataframe with a new column named
#' \code{FDR_adjusted_p}.
#'
#' @importFrom dplyr arrange bind_cols select rename mutate
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
#' @param m the mash result (from joint or 1by1 analysis); must have been
#'     computed using usepointmass=TRUE
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


#' From a mash result, get effects that are significant in at least
#'     one condition.
#'
#' @param m the mash result (from joint or 1by1 analysis)
#' @param thresh indicates the threshold below which to call signals significant
#' @param conditions which conditions to include in check (default to all)
#' @param sig_fn the significance function used to extract significance from
#'     mash object; eg could be ashr::get_lfsr or ashr::get_lfdr. (Small values
#'     must indicate significant.)
#'
#' @return a vector containing the indices of the significant effects, by
#'     order of most significant to least
#'
#' @importFrom ashr get_lfsr
#'
#' @export
get_significant_results = function(m, thresh = 0.05, conditions = NULL,
                                   sig_fn = ashr::get_lfsr) {
  if (is.null(conditions)) {
    conditions = 1:get_ncond(m)
  }
  top = apply(sig_fn(m)[,conditions,drop=FALSE],1,min) # find top effect
  # in each condition
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
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(type == "anno"){
    stopifnot("CHR" %in% names(anno_df) & "POS" %in% names(anno_df))
    anno_df <- anno_df %>%
      mutate(marker_ID = paste(.data$CHR, .data$POS, sep = "_"))
    pos_subset <- which(snp$map$marker.ID %in% anno_df$marker_ID)
  } else if(type == "range"){
    stopifnot(is.character(chr) & is.numeric(pos1) & is.numeric(pos2) & pos2 > pos1)
    if(length(which(!(chr %in% paste0("Chr0", rep(1:9, times = 2),
                                      rep(c("K", "N"), 9))))) > 0){
      stop(paste0("This function does not work on files with scaffolds; it ",
                  "needs chromosomes in the range of Chr01K to Chr09N."))
    }
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
  stopifnot(attr(snp, "class") == "bigSNP")
  subset_df <- snp$genotypes[rows_along(snp$genotypes),
                             cols_along(snp$genotypes)]
  colnames(subset_df) <- snp$map$marker.ID
  subset_tibble <- data.frame(PLANT_ID = snp$fam$sample.ID, subset_df) %>%
    as_tibble() %>%
    mutate(PLANT_ID = as.character(.data$PLANT_ID))
  return(subset_tibble)
}

#' @title Threshold a FBM
#'
#' @param X A FBM.
#' @param ind Numeric. Single column to threshold.
#' @param thr Numeric. Threshold above which row value is changed to "TRUE".
#' @param quantile Numeric. Top quantile/percentile to keep for each GWAS for
#'     comparisons.
#'
thresholdFBM <- function(X, ind, thr, quantile = NA) {
  if (!is.na(quantile)) {
    thr <- quantile(X[, ind], quantile)
  }
  case_when(X[, ind] > thr ~ TRUE,
            TRUE ~ FALSE)
}

#' @title Plot the -log10p-value distributions for two univariate GWAS against
#'     one another using ggplot.
#'
#' @description This function takes GWAS results from the dive_ functions of
#'     snpdiver that create FBM of univariate GWAS effects. It creates a
#'     dataframe from these results suitable for an Upset plot, containing only
#'     the rows/SNPs significant in at least one univariate GWAS at the
#'     -log10p threshold specified.
#'
#' @param effects fbm created using 'pvdiv_standard_gwas'.
#'     Saved under the name "gwas_effects_{suffix}.rds" and can be loaded into
#'     R using the bigstatsr function "big_attach".
#' @param snp A "bigSNP" object; load with \code{bigsnpr::snp_attach()}.
#'     Here, genomic information for Panicum virgatum.
#' @param metadata Metadata created using 'pvdiv_standard_gwas'.
#'     Saved under the name "gwas_effects_{suffix}_associated_metadata.csv"
#' @param thr Numeric. Threshold above which SNP/row is kept for comparisons.
#' @param quantile Numeric. Top quantile/percentile to keep for each GWAS for
#'     comparisons.
#'
#' @importFrom matrixStats rowMaxs
#' @importFrom tibble add_column
#'
#' @export
pvdiv_fbm_upset_df <- function(effects, snp, metadata, thr = 7, quantile = NA){
  if (attr(snp, "class") != "bigSNP") {
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if (attr(effects, "class") != "FBM") {
    stop("effects needs to be a FBM object, produced by pvdiv_standard_gwas().")
  }
  gwas_ok <- floor(effects$ncol / 3)
  if (gwas_ok != nrow(metadata)) {
    stop(paste0("metadata needs to be the dataframe saved with the FBM object",
         "produced by pvdiv_standard_gwas(). This is saved as a csv ending in",
         "associated_metadata.csv."))
  }
  ind_p <- (1:(gwas_ok))*3
  colnames_fbm <- metadata$phe
  if(effects$ncol < (sum(gwas_ok)*3 + 1))  {
    effects$add_columns(ncol_add = 1)
  }  # add a column for the threshold score if there isn't one already
  thr_df <- big_apply(effects,
                      a.FUN = function(X, ind) rowMaxs(as.matrix(effects[, ind])),
                      ind = ind_p, a.combine = 'c', block.size = 100)
  effects[,(sum(gwas_ok)*3 + 1)] <- thr_df
  thr_df <- which(thr_df > thr)
  snpfile_subset <- snp_subset(snp, ind.col = thr_df)
  snp_subset <- snp_attach(snpfile_subset)
  thr_df <- big_copy(effects, ind.row = thr_df, ind.col = ind_p)

  for (j in seq_along(1:thr_df$ncol)) {  # standardize one gwas at a time.
    thr_df[, j] <- big_apply(thr_df, a.FUN = thresholdFBM, ind = j, thr = thr,
                             a.combine = 'plus')
  }

  thr_df1 <- thr_df[1:thr_df$nrow, 1:thr_df$ncol]
  colnames(thr_df1) <- colnames_fbm
  thr_df1 <- as_tibble(thr_df1) %>%
    add_column(CHR = snp_subset$map$chromosome,
               POS = snp_subset$map$physical.pos, .before = TRUE)
  return(thr_df1)
}



#' @title Plot the -log10p-value distributions for two univariate GWAS against
#'     one another using ggplot.
#'
#' @description This function takes GWAS results from the dive_ functions of
#'     snpdiver that create FBM of univariate GWAS effects. It uses ggplot2 to
#'     plot these distributions against one another for comparison purposes.
#'
#' @param effects GWAS output saved as an FBM with an .rds and .bk file generated
#'     by dive_phe2effects or dive_phe2mash functions of snpdiver. Load this
#'     into the R environment using bigstatsr::big_attach.
#' @param metadata The metadata associated with GWAS output generated by
#'     dive_phe2effects or dive_phe2mash functions of snpdiver. Eventually
#'     this and gwas should be rolled into a new object type for R, but not yet.
#' @param e_row Integer. The row number of gwas_meta that corresponds to the
#'     expected GWAS. This will be plotted on the x-axis for all comparisons.
#' @param o_row Integer vector. The rownumbers of gwas_meta that are the
#'     observed GWAS. These will be compared to the expected GWAS in e_row.
#' @param thr Numeric. -log10p threshold. Only SNPs with an expected -log10p
#'     value above this threshold will be plotted.
#' @param suffix Optional character vector to give saved files a unique search
#'     string/name.
#' @param outputdir String or file.path() to the output directory. Default is
#'     the working directory.
#'
#' @return Plots saved to disk in a "analysis/gwas_comps" folder comparing the
#'     expected distribution, e_row, to all observed gwas distributions, o_row.
#'
#' @importFrom dplyr filter
#' @import ggplot2
#' @importFrom cowplot save_plot
#' @importFrom rlang .data
#' @import hexbin
#'
#' @export
pvdiv_fbm_qq <- function(effects, metadata, e_row = 1,
                             o_row = 2:nrow(metadata), thr = 5,
                             suffix = NA, outputdir = ".") {
  requireNamespace("hexbin")
  e_row <- e_row*3
  o_row <- o_row*3
  if (!dir.exists(outputdir)) {
    stop(paste0("Output directory (outputdir) specified does not exist. ",
                "Please specify an existing output directory."))
  }
  if (!grepl("^_", suffix) & !is.na(suffix)){
    suffix <- paste0("_", suffix)
  } else if(is.na(suffix)) {
    suffix <- ""
  }
  for (i in seq_along(o_row)) {
    plot_qq <- filter(data.frame(effects[,c(e_row, o_row[i])]),
                      .data$X1 >= thr) %>%
      ggplot() +
      geom_hex(aes(x = .data$X1, y = .data$X2)) +
      scale_fill_gradient(trans = "log") +
      geom_smooth(aes(x = .data$X1, y = .data$X2), ) +
      geom_abline(slope = 1, linetype = 2, color = "red") +
      theme(legend.position = "right") +
      labs(x = metadata[e_row/3,1], y = metadata[o_row[i]/3,1])

    save_plot(filename = file.path(outputdir,
                                    paste0("GWAS_significance_thr_log10p_", thr,
                                           "_", metadata[e_row/3,1], "_vs_",
                                           metadata[o_row[i]/3,1], suffix, ".png")),
              plot = plot_qq, base_height = 5, base_asp = 1.61)
  }
}
