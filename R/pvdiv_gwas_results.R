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
#' @import rlang
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
#'
#' @export
get_marker_df <- function(m){
  marker_df <- get_significant_results(m, thresh = 1) %>%
    enframe(name = "Marker") %>%
    arrange(.data$value)

  return(marker_df)
}


#' @title Create a Table of Annotated Top SNPs for Panicum virgatum
#'
#' @description Loads in one set of GWAS results from a bigsnp GWAS, a mash
#'    object, RQTL2 output, or a dataframe with 'CHR', 'start', and 'end'
#'    columns for the genome of Panicum virgatum. Then, it constructs SNP
#'    tables meeting different criteria using these genomic intervals.
#'
#' @param df a data frame or tbl_df. Can be bigsnpr output, mashr output
#'    (loaded into R), r/qtl2 output (specify the path to a saved .csv file),
#'    or, for Panicum virgatum intervals in another format, a
#'    data frame containing columns 'CHR', 'start', and 'end'.
#' @param type Type of Panicum virgatum genomic marker input specified by the
#'    df parameter. Options are "bigsnp", "mash", "rqtl2", and "df". Defaults
#'    to 'bigsnp'.
#' @param n An integer or integer vector The numberof most significant SNPs to
#'     select (by p-value). Set to NA to omit this table. Default is 10.
#' @param FDRalpha The false discovery rate. Numeric, a number or vector of
#'     numbers between 0 and 1. Set to NA to omit this table. Default is 0.1.
#' @param rangevector How far from the significant SNP should annotations be
#'    pulled? Can be an integer or a vector of integers. Default is 0 (the SNP
#'    itself) and a 10 kbp window around the SNP.
#' @param markers For data frames of type "mash" or "bigsnp", the same set of
#'    markers (with CHR and POS columns) as your df object.
#' @param anno_info Gene information from
#'    Pvirgatum_516_v5.1.annotation_info.txt
#' @param txdb Annotation information from Pvirgatum_516_v5.1.gene.txdb.sqlite
#'
#' @return A list containing dataframes of SNPs. If more than one dataframe is
#'     returned, they are named using the criteria used to select the SNPs in the dataframe.
#'
#' @importFrom AnnotationDbi loadDb
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom ashr get_lfsr
#' @importFrom dplyr arrange select rename mutate top_n left_join filter
#' @importFrom GenomicFeatures genes
#' @importFrom magrittr %>%
#' @importFrom multtest mt.rawp2adjp
#' @importFrom readr read_csv
#' @import rlang
#' @importFrom rlist list.append
#' @importFrom stats predict
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom tidyr unite separate
#' @importFrom VariantAnnotation locateVariants AllVariants PromoterVariants IntergenicVariants
#'
#' @export
pvdiv_table_topsnps <- function(df, type = c("bigsnp", "mash", "rqtl2", "df"),
                                n = 10, FDRalpha = 0.1,
                                rangevector = c(0, 10000), markers = NULL,
                                anno_info = NULL, txdb = NULL){
  n <- as.integer(n)
  FDRalpha <- as.numeric(FDRalpha)
  rangevector <- as.integer(rangevector)
  stopifnot(type %in% c("bigsnp", "mash", "rqtl2", "df"), !is_null(anno_info),
            !is_null(txdb))
  if(type %in% c("bigsnp", "mash") & is.na(n) & is.na(FDRalpha)){
    stop(paste0("For 'mash' and 'bigsnp' types, need to specify at least one",
                "of n (as an integer) or FDR (between 0 and 1)."))
  }
  topsnp_inputlist <- list()
  ## Prepare a dataframe for each type for further analysis.
  if(type == "bigsnp"){
    stopifnot(!is_null(markers))

    input_df <- df %>%
      mutate(p.value = predict(df, log10 = FALSE))
    res <- multtest::mt.rawp2adjp(input_df$p.value, alpha = FDRalpha,
                        proc = "BH")
    adj_p <- res$adjp[order(res$index), ]
    input_df <- cbind(markers, input_df, adj_p) %>%
      as_tibble() %>%
      dplyr::rename(`p value` = .data$p.value,
                    `FDR Adjusted p value` = .data$`BH`,
                    `SNP Effect` = .data$estim,
                    `SNP standard error` = .data$std.err)

  }
  if(type == "mash"){
    stopifnot(!is_null(markers))

    input_df <- get_log10bf(m = df) %>%
      as.data.frame() %>%
      rownames_to_column(var = "value") %>%
      mutate(value = as.integer(.data$value)) %>%
      as_tibble() %>%
      left_join(get_marker_df(m = df), by = "value") %>%
      dplyr::rename(log10BayesFactor = .data$V1) %>%
      dplyr::select(-.data$value) %>%
      separate(.data$Marker, into = c("CHR", "POS"), sep = "_") %>%
      mutate(POS = as.integer(.data$POS))
  }

  if(type == "rqtl2"){

    input <- read_csv(file = df)
    topsnp_inputlist[[1]] <- input %>%
      separate(.data$marker, into = c("CHR", "marker_pos"), sep = "_") %>%
      separate(.data$flank_lo, into = c("cl", "flank_lo_pos"), sep = "_") %>%
      separate(.data$flank_hi, into = c("ch", "flank_hi_pos"), sep = "_") %>%
      dplyr::select(-.data$chr, -.data$cl, -.data$ch) %>%
      mutate(POS = as.numeric(.data$marker_pos)*1E6,
             start = as.numeric(.data$flank_lo_pos)*1E6,
             end = as.numeric(.data$flank_hi_pos)*1E6,
             POS = as.integer(.data$POS),
             start = as.integer(.data$start),
             end = as.integer(.data$end)
      ) %>%
      dplyr::select(-(.data$marker_pos:.data$flank_hi_pos))
  }
  if(type == "df"){
    stopifnot(c("CHR", "start", "end") %in% names(df))
    topsnp_inputlist[[1]] <- df %>%
      mutate(start = as.integer(.data$start),
             end = as.integer(.data$end))
  }
  if(type %in% c("bigsnp", "mash")){
    ## Now make multiple dataframes for different categories of top SNPs.

    ### Get as_tibble for each SNP criterion specified (both a number of top SNPs
    ### and FDR of some thresholds are supported)

    if(!is.na(n[1]) & !is.na(FDRalpha[1])){
      for(i in seq_along(n)){
        if(type == "bigsnp"){
          topsnp_inputlist[[i]] <- input_df %>%
            top_n( -n[i], .data$`p value`)
        } else {
          topsnp_inputlist[[i]] <- input_df %>%
            top_n( n[i], .data$log10BayesFactor)
        }
      }
      for(i in seq_along(FDRalpha)){
        if(type == "bigsnp"){
          FDR1 <- input_df %>%
            filter(.data$`FDR Adjusted p value` <= FDRalpha[i])
        } else {
          BFalpha <- -log10(FDRalpha[i])+1
          FDR1 <- input_df %>%
            filter(.data$log10BayesFactor >= BFalpha)
        }
        topsnp_inputlist <- rlist::list.append(topsnp_inputlist, FDR1)
      }
      names1 <- c(paste0("top", n, "SNPs_"), paste0("FDR", FDRalpha, "_"))
    } else if(!is.na(n[1]) & is.na(FDRalpha[1])){
      for(i in seq_along(n)){
        if(type == "bigsnp"){
          topsnp_inputlist[[i]] <- input_df %>%
            top_n( -n[i], .data$`p value`)
        } else {
          topsnp_inputlist[[i]] <- input_df %>%
            top_n( n[i], .data$log10BayesFactor)
        }
      }
      names1 <- c(paste0("top", n, "SNPs_"))
    } else if (is.na(n[1]) & !is.na(FDRalpha[1])){
      for(i in seq_along(FDRalpha)){
        if(type == "bigsnp"){
          FDR1 <- input_df %>%
            filter(.data$`FDR Adjusted p value` <= FDRalpha[i])
        } else {
          BFalpha <- -log10(FDRalpha[i])+1
          FDR1 <- input_df %>%
            filter(.data$log10BayesFactor >= BFalpha)
        }
        topsnp_inputlist <- rlist::list.append(topsnp_inputlist, FDR1)
      }
      names1 <- c(paste0("FDR", FDRalpha, "_"))
    } else {
      stop(paste0("Need to specify at least one of n (as an integer) or FDR",
                  " (between 0 and 1)."))
    }
  }

  #### Find annotations for each tbl_df found for each SNP criterion specified.
  for(k in seq_along(rangevector)){
    range <- rangevector[k]
    for(i in seq_along(topsnp_inputlist)){
      loop_input <- topsnp_inputlist[[i]]
      if(nrow(loop_input) > 0){

        if(type %in% c("bigsnp", "mash")){
          # Prepare input dataframe
          input <- loop_input %>%
            mutate(start = .data$POS - (range/2),
                   end = .data$POS + (range/2))
        } else {
          input <- loop_input # Other types have start and end already
        }

        ## Make input into a GRanges object for querying with locateVariants
        target <- with(input, GRanges(seqnames = Rle(.data$CHR),
                                      ranges = IRanges(start = .data$start,
                                                       end = .data$end,
                                                       names = NULL),
                                      strand = Rle(strand("*"))))

        ### Find genes that overlap input SNPs with VariantAnnotation
        loc <-
          locateVariants(target, txdb,
                         AllVariants(promoter =
                                       PromoterVariants(upstream = 2000L,
                                                        downstream = 2000L),
                                     intergenic =
                                       IntergenicVariants(upstream = 20000L,
                                                          downstream = 20000L,
                                                          idType = "gene")))

        #### Convert a GRanges object to an output data frame
        out <- as_tibble(loc)

        filter1 <- list(gene_id = out$GENEID)
        genedf <- as_tibble(genes(txdb, filter = filter1)) %>%
          dplyr::rename(gene_width = .data$width, gene_start = .data$start, gene_end = .data$end,
                        gene_strand = .data$strand) %>%
          mutate(seqnames = as.character(.data$seqnames))

        topsnp_outputlist[[(i + (k-1)*length(topsnp_inputlist))]] <- out %>%
          mutate(seqnames = as.character(.data$seqnames)) %>%
          dplyr::select(.data$seqnames:.data$TXID, .data$GENEID) %>%
          left_join(input, by = c("seqnames" = "CHR", "start", "end")) %>%
          left_join(genedf, by = c("seqnames", "GENEID" = "gene_id")) %>%
          left_join(anno_info, by = c("GENEID" = "locusName")) %>%
          dplyr::rename(region_start = .data$start,
                        region_end = .data$end,
                        region_strand = .data$strand) %>%
          mutate(d2_start = .data$POS - .data$gene_start,
                 d2_end = .data$POS - .data$gene_end,
                 distance = ifelse(abs(.data$d2_start) < abs(.data$d2_end),
                                   abs(.data$d2_start),
                                   abs(.data$d2_end))) %>%
          dplyr::rename(CHR = .data$seqnames,
                        `Annotation` = .data$LOCATION,
                        `Gene ID` = .data$GENEID,
                        `Arabidopsis thaliana homolog` = .data$`Best-hit-arabi-name`,
                        `A. thaliana gene name` = .data$`arabi-symbol`,
                        `A. thaliana gene annotation` = .data$`arabi-defline`,
                        `Rice homolog` = .data$`Best-hit-rice-name`,
                        `Rice gene name` = .data$`rice-symbol`,
                        `Rice gene annotation` = .data$`rice-defline`,
                        `Panther Categories` = .data$Panther,
                        `distance from gene` = .data$distance) %>%
          dplyr::select(-.data$LOCSTART, -.data$LOCEND)
      } else {
        topsnp_outputlist[[(i + (k-1)*length(topsnp_inputlist))]] <-
          as_tibble(NA)
      }
    }
  }

  if(type %in% c("mash", "bigsnp")){
    names(topsnp_outputlist) <- paste0(rep(names1, length(rangevector)),
                                       "within",
                                       rep(rangevector, each = length(names1)),
                                       "bp")
  }

  return(topsnp_outputlist)
}





