#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @importFrom multtest mt.rawp2adjp
#' @importFrom stats predict
bigsnp2anno <- function(df, markers, FDRalpha){
  input_df <- df %>%
    mutate(p.value = predict(df, log10 = FALSE))
  if(!is.na(FDRalpha)[1]){
    res <- mt.rawp2adjp(input_df$p.value, alpha = FDRalpha,
                        proc = "BH")
    adj_p <- res$adjp[order(res$index), ]
    input_df <- cbind(markers, input_df, adj_p) %>%
      as_tibble() %>%
      dplyr::rename(`p value` = .data$p.value,
                    `FDR Adjusted p value` = .data$`BH`,
                    `SNP Effect` = .data$estim,
                    `SNP standard error` = .data$std.err)
  } else {
    input_df <- cbind(markers, input_df) %>%
      as_tibble() %>%
      dplyr::rename(`p value` = .data$p.value,
                    `SNP Effect` = .data$estim,
                    `SNP standard error` = .data$std.err)
  }
  return(input_df)
}

#' @importFrom ashr get_lfsr
#' @importFrom dplyr rename select left_join mutate
#' @importFrom magrittr %>%
#' @importFrom tidyr separate
#' @importFrom tibble as_tibble rownames_to_column
mash2anno <- function(df, markers){
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
  return(input_df)
}

#' @importFrom dplyr select mutate
#' @importFrom magrittr %>%
#' @importFrom readr read_csv
#' @importFrom tidyr separate
rqtl2anno <- function(df){
  input <- read_csv(file = df)
  input_df <- input %>%
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
  return(input_df)
}

#' @importFrom dplyr filter top_n
#' @importFrom rlist list.append
#' @importFrom magrittr %>%
get_top_snps <- function(input_df, n, FDRalpha, type = c("bigsnp", "mash")){
  topsnp_inputlist <- list()
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
      topsnp_inputlist <- list.append(topsnp_inputlist, FDR1)
    }
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
      topsnp_inputlist <- list.append(topsnp_inputlist, FDR1)
    }
  } else {
    stop(paste0("Need to specify at least one of n (as an integer) or FDR",
                " (between 0 and 1)."))
  }
  return(topsnp_inputlist)
}

#' @importFrom dplyr select mutate left_join rename
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
get_tidy_annos <- function(df, input, anno_info, txdb){
  out <- as_tibble(df)
  requireNamespace("GenomicFeatures")

  filter1 <- list(gene_id = out$GENEID)
  genedf <- as_tibble(GenomicFeatures::genes(txdb, filter = filter1)) %>%
    dplyr::rename(gene_width = .data$width, gene_start = .data$start,
                  gene_end = .data$end, gene_strand = .data$strand) %>%
    mutate(seqnames = as.character(.data$seqnames))

  outdf <- out %>%
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
   return(outdf)
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
#'    df parameter. Options are "bigsnp", "mash", "rqtl2", and "table". Defaults
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
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom multtest mt.rawp2adjp
#' @importFrom rlang is_null
#' @importFrom tibble as_tibble
#'
#' @export
pvdiv_table_topsnps <- function(df, type = c("bigsnp", "mash", "rqtl2", "table"),
                                n = 10, FDRalpha = 0.1,
                                rangevector = c(0, 10000), markers = NULL,
                                anno_info = NULL, txdb = NULL){
  requireNamespace("VariantAnnotation")

  n <- as.integer(n)
  FDRalpha <- as.numeric(FDRalpha)
  rangevector <- as.integer(rangevector)
  stopifnot(type %in% c("bigsnp", "mash", "rqtl2", "table"), !is_null(anno_info),
            !is_null(txdb))
  if(type == "table" & !("CHR" %in% names(df)) & !("start" %in% names(df)) &
     !("end" %in% names(df))){
    stop(paste0("For 'table' type, need to have columns 'CHR', 'start', and ",
                "'end' in your data frame."))
  }
  if(type %in% c("bigsnp", "mash") & is.na(n)[1] & is.na(FDRalpha)[1]){
    stop(paste0("For 'mash' and 'bigsnp' types, need to specify at least one",
                "of n (as an integer) or FDR (between 0 and 1)."))
  }
  if(type %in% c("bigsnp", "mash")){
    stopifnot(!is_null(markers[1]))
  }
  topsnp_inputlist <- list()
  topsnp_outputlist <- list()
  ## Prepare a dataframe for each type for further analysis.
  if(type == "bigsnp"){
    input_df <- bigsnp2anno(df = df, markers = markers, FDRalpha = FDRalpha)
  }
  if(type == "mash"){
    input_df <- mash2anno(df = df, markers = markers)
  }
  if(type == "rqtl2"){
    topsnp_inputlist[[1]] <- rqtl2anno(df = df)
  }
  if(type == "table"){
    topsnp_inputlist[[1]] <- df %>%
      mutate(start = as.integer(.data$start),
             end = as.integer(.data$end),
             POS = (.data$end + .data$start)/2)
  }
  if(type %in% c("bigsnp", "mash")){
    ## Now make multiple dataframes for different categories of top SNPs.
    topsnp_inputlist <- get_top_snps(input_df = input_df, n = n,
                                     FDRalpha = FDRalpha, type = type)
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
        target <- with(input, GRanges(seqnames = Rle(CHR),
                                      ranges = IRanges(start = start,
                                                       end = end,
                                                       names = NULL),
                                      strand = Rle(strand("*"))))

        ### Find genes that overlap input SNPs with VariantAnnotation
        loc <-
          VariantAnnotation::locateVariants(target, txdb,
                                            VariantAnnotation::AllVariants(promoter =
                                                                             VariantAnnotation::PromoterVariants(upstream = 2000L,
                                                        downstream = 2000L),
                                     intergenic =
                                       VariantAnnotation::IntergenicVariants(upstream = 20000L,
                                                          downstream = 20000L,
                                                          idType = "gene")))

        #### Convert a GRanges object to an output data frame
        topsnp_outputlist[[(i + (k-1)*length(topsnp_inputlist))]] <-
          get_tidy_annos(df = loc, input = input, anno_info = anno_info,
                         txdb = txdb)
      } else {
        topsnp_outputlist[[(i + (k-1)*length(topsnp_inputlist))]] <-
          as_tibble(NA)
      }
    }
  }
  if(type %in% c("bigsnp", "mash")){
    if(!is.na(n)[1] & !is.na(FDRalpha)[1]){
    names1 <- c(paste0("top", n, "SNPs_"), paste0("FDR", FDRalpha, "_"))
    } else if(!is.na(n)[1]){
    names1 <- c(paste0("top", n, "SNPs_"))
    } else {
    names1 <- c(paste0("FDR", FDRalpha, "_"))
    }
  names(topsnp_outputlist) <- paste0(rep(names1, length(rangevector)),
                                     "within",
                                     rep(rangevector, each = length(names1)),
                                     "bp")
  }
  return(topsnp_outputlist)
}
