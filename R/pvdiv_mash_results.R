#' @title Get column names from a mash object
#'
#' @description This function extracts the column names from the local false
#' sign rate table of a mash object's results. This can tell you the condition
#' names or phenotype names used in the mash object. That can be useful for
#' looking at a subset of these columns, say.
#'
#' @param m An object of type mash
#'
#' @return A vector of phenotype names
#'
#' @examples
#'     \dontrun{get_colnames(m = mash_obj)}
#'
#' @export
get_colnames <- function(m){
  column_names <- colnames(m$result$lfsr)
  return(column_names)
}

#' Count number of conditions each effect is significant in
#'
#' @param m the mash result (from joint or 1by1 analysis)
#' @param thresh indicates the threshold below which to call signals significant
#' @param conditions which conditions to include in check (default to all)
#' @param sig_fn the significance function used to extract significance from mash object; eg could be ashr::get_lfsr or ashr::get_lfdr
#'
#' @return a vector containing the number of significant conditions
#'
#' @export
get_n_significant_conditions = function(m, thresh = 0.05, conditions = NULL,
                                        sig_fn = get_lfsr){
  if (is.null(conditions)) {
    conditions = 1:get_ncond(m)
  }
  return(apply(sig_fn(m)[,conditions,drop=FALSE] < thresh, 1, sum))
}



#' @title Reorder correlation matrix
#'
#' @description  Reorder correlation coefficients from a matrix of things
#'     (including NA's) and hierarchically cluster them
#'
#' @param cormat A correlation matrix
#'
#' @importFrom cluster daisy
#' @importFrom stats hclust
#'
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- daisy(cormat, metric = "gower")
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}


# --- Plot & Save Plots ---------

#' @title Significant SNPs per number of conditions
#'
#' @description For some number of columns in a mash object that correspond to
#'     conditions, find the number of SNPs that are significant for that number
#'     of conditions.
#'
#' @param m An object of type mash
#' @param conditions A vector of conditions
#' @param saveoutput Logical. Save plot output to a file? Default is FALSE.
#' @param thresh What is the threshold to call an effect significant? Default is
#'     0.05.
#'
#' @return A list containing a dataframe of the number of SNPs significant per
#'     number of conditions, and a ggplot object using that dataframe.
#'
#' @import ggplot2
#' @importFrom tibble enframe
#' @importFrom dplyr rename summarise filter group_by
#'
#' @examples
#'   \dontrun{mash_plot_sig_by_condition(m = mash_obj, saveoutput = TRUE)}
#'
#' @export
mash_plot_sig_by_condition <- function(m, conditions = NA, saveoutput = FALSE,
                                       thresh = 0.05){

  thresh <- as.numeric(thresh)
  num_sig_in_cond <- c()

  if(is.na(conditions)[1]){
    cond <- get_colnames(m = m)
  }

  SigHist <- get_n_significant_conditions(m = m, thresh = thresh,
                                          conditions = cond) %>%
    enframe(name = "Marker") %>%
    rename(Number_of_Conditions = .data$value) %>%
    group_by(.data$Number_of_Conditions) %>%
    summarise(Significant_SNPs = n()) %>%
    filter(.data$Number_of_Conditions != 0)

  vis <- ggplot(SigHist, aes(x = .data$Number_of_Conditions, y = .data$Significant_SNPs)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    xlab(label = "Number of Conditions") +
    ylab(label = "Number of Significant SNPs")

  if(saveoutput == TRUE){
    ggsave(paste0("SNPs with significant effects in n conditions ",
                  str_replace_all(Sys.time(), ":", "."),
                  ".bmp"), width = 5, height = 3, units = "in", dpi = 400)
  }

  return(list(sighist = SigHist, ggobject = vis))
}


#' @title Manhattan plot in ggplot colored by significant conditions
#'
#' @description Takes a mash object and, for some vector of phenotypes, returns
#'     a Manhattan plot ggplot object (and its dataframe). Each SNP in the plot
#'     is colored by the number of phenotypes it is significant for. Even and
#'     odd chromosomes have different shapes for their SNPs, so that
#'     chromosome identity can be determined.
#'
#' @param m A mash object (outputted by mash).
#' @param cond A vector of phenotypes. Defaults to the names of each
#'     column in the mash object.
#' @param saveoutput Logical. Should the output be saved to the path?
#' @param thresh Numeric. The threshold used for the local false sign rate to
#'     call significance in a condition.
#'
#' @return A \code{tbl_df()} of the data used to make the Manhattan plot, and a
#'     ggplot object containing the Manhattan.
#'
#' @importFrom cowplot save_plot
#' @importFrom dplyr rename select arrange mutate left_join
#' @import ggplot2
#' @importFrom tibble as_tibble rownames_to_column enframe
#' @importFrom tidyr separate
#' @import viridis
#' @importFrom stringr str_replace_all
#'
#' @examples
#' \dontrun{manhattan_out <- mash_ggman_by_condition(m = m, saveoutput = TRUE)}
#'
#' @export
mash_plot_manhattan_by_condition <- function(m, cond = NA,
                                             saveoutput = FALSE, thresh = 0.05){
  num_sig_in_cond <- c()

  if(is.na(cond)[1]){
    cond <- get_colnames(m = m)
  }

  log10bf_df <- get_log10bf(m = m) %>%
    as.data.frame() %>%
    rownames_to_column(var = "value") %>%
    mutate(value = as.integer(.data$value)) %>%
    as_tibble() %>%
    left_join(get_marker_df(m = m)) %>%
    dplyr::rename(log10BayesFactor = .data$V1) %>%
    dplyr::select(-.data$value)

  ggman_df <- get_n_significant_conditions(m = m, thresh = thresh,
                                           conditions = cond) %>%
    enframe(name = "Marker") %>%
    rename(Num_Sig_Conditions = .data$value) %>%
    separate(.data$Marker, into = c("Chr", "Pos"), remove = FALSE, sep = "_",
             extra = "merge") %>%
    mutate(Pos = as.numeric(.data$Pos)) %>%
    left_join(log10bf_df, by = "Marker") %>%
    arrange(.data$Chr, .data$Pos)

  log10BF <- expression(paste("log"[10], plain("(Bayes Factor)")))

  ggmanobject <- ggplot(data = ggman_df, aes(x = .data$Pos, y = .data$log10BayesFactor)) +
    geom_point(aes(color = .data$Num_Sig_Conditions, fill = .data$Num_Sig_Conditions,
                   shape = as.factor(.data$Chr))) +
    facet_wrap(~ .data$Chr, nrow = 1, scales = "free_x", strip.position = "bottom") +
    scale_color_viridis(option = "B") + scale_fill_viridis(option = "B") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill=NA)) +
    labs(x = "Chromosome", y = log10BF) +
    scale_x_continuous(expand = c(0.3, 0.3)) +
    scale_shape_manual(values = rep(c(21,22),9), guide = FALSE)

  if(saveoutput == TRUE){
    save_plot(paste0("Manhattan_mash_", str_replace_all(Sys.time(), ":", "."),
                     ".png"), plot = ggmanobject, base_aspect_ratio = 2.5,
              base_height = 4)
  }

  return(list(ggman_df = ggman_df, ggmanobject = ggmanobject))
}


