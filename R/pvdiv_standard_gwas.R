#' @title Get current date-time in a filename-appropriate format.
#'
#' @description Converts the current \code{Sys.time()} system time to a format
#'     that is acceptable to include in a filename. Changes punctuation that
#'     won't work in a filename.
#'
#' @return A string containing the current date-time with spaces and colons
#'     replaced with underscores and periods, respectively.
#'
#' @importFrom stringr str_replace_all
get_date_filename <- function(){
  str_replace_all(str_replace_all(Sys.time(), ":", "."), " ", "_")
}

#' @title Wrapper for the snp_autoSVD function for switchgrass.
#'
#' @description This is a wrapper to determine population structure for GWAS
#'     for a SNP file with the switchgrass chromosomes, which are not numeric.
#'     Arguments that are recognized by bigsnpr::snp_autoSVD can also be
#'     specified in this function.
#'
#' @param snp A "bigSNP" object; load with bigsnpr::snp_attach().
#' @param k Integer. The number of principal components to find. Default is 10.
#' @param ncores Integer. Number of cores to use. Default is one.
#' @param saveoutput Logical. Should the output be saved to the working
#'     directory?
#' @param ... Other arguments to \code{\link{snp_autoSVD}}.
#'
#' @return A big_SVD object.
#'
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when
#' @importFrom tibble enframe
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' svd20 <- pvdiv_autoSVD(snp = snp, k = 20, saveoutput = TRUE)
#' }
#'
#' @export
pvdiv_autoSVD <- function(snp, k = 10, ncores = 1, saveoutput = FALSE, ...){
  requireNamespace("dots")
  fun.scaling <- dots::dots(name = 'fun.scaling', value = snp_scaleBinom(), ...)
  thr.r2 <- dots::dots(name = 'thr.r2', value = 0.2, ...)
  size <- dots::dots(name = 'size', value = 100/thr.r2, ...)
  roll.size <- dots::dots(name = 'roll.size', value = 50, ...)
  int.min.size <- dots::dots(name = 'int.min.size', value = 20, ...)
  #alpha.tukey <- dots::dots(name = 'alpha.tukey', value = 0.05, ...)
  #min.mac <- dots::dots(name = 'min.mac', value = 10, ...)
  #max.iter <- dots::dots(name = 'max.iter', value = 5, ...)
  verbose <- dots::dots(name = 'verbose', value = FALSE, ...)

  # Argument error checks; Set up numeric chromosome data frame.
  stopifnot(attr(snp, "class") == "bigSNP")
  G <- snp$genotypes
  CHR <- snp$map$chromosome
  POS <- snp$map$physical.pos
  plants <- snp$fam$sample.ID
  if(length(which(!(CHR %in% paste0("Chr0", rep(1:9, times = 2),
                                    rep(c("K", "N"), 9))))) > 0){
    stop(paste0("This function does not work on files with scaffolds; it ",
    "needs chromosomes in the range of Chr01K to Chr09N."))
  }
  if(saveoutput == FALSE){
    message("'saveoutput' is FALSE, so the svd will not be saved to the working directory.")
  }
  CHRN <- enframe(CHR, name = NULL, value = "CHR") %>%
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

  # Determine population structure
  svd <- snp_autoSVD(G = G, infos.chr = CHRN$CHRN, infos.pos = POS,
                     ncores = ncores, k = k, fun.scaling = fun.scaling,
                     thr.r2 = thr.r2, size = size, roll.size = roll.size,
                     int.min.size = int.min.size, #alpha.tukey = alpha.tukey,
                     #min.mac = min.mac, max.iter = max.iter,
                     verbose = verbose)
  if(saveoutput){
  saveRDS(svd, file = paste0("SVD_", length(plants), "g_", k, "PCs.rds"))
  }
  return(svd)
}


#' @title Find a kinship matrix using the van Raden method.
#'
#' @description Calculate the kinship matrix using code from GAPIT and bigsnpr
#'     and the methods of VanRaden (2009, J. Dairy Sci. 91:4414???C4423). Note
#'     that this matrix cannot currently be used with the GWAS methods in
#'     bigsnpr; however, this matrix could be used with other GWAS programs.
#'
#' @param snp A "bigSNP" object; load with bigsnpr::snp_attach().
#' @param ind.row An integer vector of the rows (individuals) to find a
#'     kinship matrix for. Defaults to all rows.
#' @param saveoutput Logical. Should the output be saved to the working
#'     directory?
#'
#' @return A kinship matrix with labeled rows and columns.
#'
#' @import bigsnpr
#' @import bigstatsr
#'
#' @examples
#' \dontrun{
#' K <- pvdiv_kinship(snp = snp, saveoutput = TRUE)
#' }
#'
#' @export
pvdiv_kinship <- function(snp, ind.row = NA, saveoutput = FALSE){
  stopifnot(attr(snp, "class") == "bigSNP")
  if(saveoutput == FALSE){
    message(paste0("'saveoutput' is FALSE, so the kinship matrix will not ",
                   "be saved to the working directory."))
  }

  G <- snp$genotypes
  if(!is.na(ind.row[1])){
    nInd <- length(ind.row)
  } else {
    nInd <- snp$genotypes$nrow
  }
  # Centered (mean of rows subtracted) transposed crossproduct of snp file.
  K <- big_tcrossprodSelf(G, ind.row = ind.row,
                          fun.scaling = big_scale(center = TRUE,
                                                  scale = FALSE))

  #Extract diagonals
  i = 1:nInd
  j = (i - 1)*nInd
  index = i + j
  d = K[index]
  DL = min(d)
  DU = max(d)
  floor = min(K[i, i])

  K = (K[i, i] - floor)/(DL - floor)
  MD = (DU - floor)/(DL - floor)

  if(!is.na(ind.row[1])){
    rownames(K) <- snp$fam$sample.ID[ind.row]
    colnames(K) <- snp$fam$sample.ID[ind.row]
  } else {
    rownames(K) <- snp$fam$sample.ID
    colnames(K) <- snp$fam$sample.ID
  }

  if(MD > 2){
    K[index] <- K[index]/(MD-1)+1
  }

  if(saveoutput){
    saveRDS(K, paste0("Kinship_van_Raden_method_", nInd, "_individuals_",
    ".rds"))
  }
  return(K)
}

#' @title Wrapper for Juenger lab standard GWAS functions for switchgrass.
#'
#' @description This is a wrapper to
#'
#' @param snp A "bigSNP" object; load with bigsnpr::snp_attach(). Here, genomic
#'     information for Panicum virgatum. Contact tjuenger <at> utexas <dot> edu
#'     to obtain this information pre-publication.
#' @param df Dataframe of phenotypes where the first column is PLANT_ID.
#' @param type Character string. Type of univarate regression to run for GWAS.
#'     Options are "linear" or "logistic".
#' @param covar Optional covariance matrix to include in the regression. You
#'     can generate these using \code{pvdiv_autoSVD()}.
#' @param ncores Number of cores to use. Default is one.
#' @param lambdagc Default is TRUE - should lambda_GC be used to find the best
#'     population structure correction? Alternatively, you can provide a data
#'     frame containing "NumPCs" and the phenotype names containing lambda_GC
#'     values. This is saved to the output directory by pvdiv_standard_gwas and
#'     saved or generated by pvdiv_lambda_GC.
#' @param outputdir String or file.path() to the output directory. Default is
#'     the working directory.
#' @param savegwas Logical. Should the gwas output be saved as a rds to the
#'     working directory? These files are typically quite large. Default is
#'     FALSE.
#' @param saveplots Logical. Should Manhattan and QQ-plots be generated and
#'     saved to the working directory? Default is TRUE.
#' @param saveannos Logical. Should annotation tables for top SNPs be generated
#'     and saved to the working directory? Default is FALSE. Can take
#'     additional arguments; requires a txdb.sqlite object used in
#'     AnnotationDbi.
#' @param txdb A txdb object such as 'Pvirgatum_516_v5.1.gene.txdb.sqlite'.
#'     Load this into your environment with AnnotationDbi::loadDb.
#' @param ... Other arguments to \code{\link{pvdiv_lambda_GC}} or
#'     \code{\link{pvdiv_table_topsnps}}.
#'
#' @return A big_SVD object.
#'
#' @import bigsnpr
#' @import bigstatsr
#' @import ggplot2
#' @importFrom data.table := data.table
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when
#' @importFrom tibble enframe as_tibble tibble
#' @importFrom rlang .data
#' @importFrom readr write_rds
#' @importFrom stats p.adjust
#' @importFrom tidyselect all_of
#'
#' @examples
#' \dontrun{
#' pvdiv_standard_gwas(snp, df = phenotypes, type = "linear", covar = svd,
#'     ncores = NCORES, lambdaGC = TRUE, savegwas = TRUE, saveplots = TRUE,
#'     saveannos = TRUE, txdb = txdb)
#' }
#'
#' @export
pvdiv_standard_gwas <- function(snp, df = switchgrassGWAS::phenotypes,
                                type = c("linear", "logistic"), covar = NULL,
                                ncores = 1, lambdagc = TRUE, outputdir = ".",
                                savegwas = FALSE, saveplots = TRUE,
                                saveannos = FALSE, txdb = NULL, ...){

  stopifnot(attr(snp, "class") == "bigSNP")
  if(colnames(df)[1] != "PLANT_ID"){
    stop("First column of phenotype dataframe (df) must be 'PLANT_ID'.")
  }
  stopifnot(type %in% c("linear", "logistic"))
  if(is.null(covar)){
    stop(paste0("Need to specify covariance matrix (covar) - you can generate this using pvdiv_autoSVD()."))
  }
  if(saveannos == TRUE & is.null(txdb)){
    stop(paste0("Need to specify a txdb object created using AnnotationDbi ",
                "in order to generate data frames containing annotated top ",
                "SNPs. If you don't have this, set saveannos = FALSE."))
  }
  if(lambdagc == TRUE){
    message("'lambdagc' is TRUE, so lambda_GC will be used to find the best population structure correction using the covariance matrix.")
  }
  if(savegwas == FALSE){
    message("'savegwas' is FALSE, so the gwas results will not be saved to disk.")
  }

  plants <- snp$fam$sample.ID
  bonferroni <- -log10(0.05/length(snp$map$physical.pos))
  markers <- tibble(CHR = snp$map$chromosome, POS = snp$map$physical.pos)

  for(i in 2:nrow(df)){

# ------- Choose best pop structure correction --------------------------

    df1 <- df %>%
      dplyr::select(.data$PLANT_ID, all_of(i))

    lambdagc <- pvdiv_lambda_GC(df = df1, type = "linear", snp = snp,
                                covar = svd, ncores = ncores,
                                npcs = c(0:20), saveoutput = TRUE)


    # ------ Run GWAS with best pop structure correction -----
    # Use a new function asr_best_PC_df(df) to find the NumPCs where lambda_GC
    # is closest to one.

    PCdf <- pvdiv_best_PC_df(lambdagc) # asv_best_PC_df(lambdagc)
    PCdf1 <- PCdf[1,]
    # PCdf1 <- data.frame(NumPCs = 20)

    if(PCdf1$NumPCs == 0){
      gwas <- pvdiv_gwas(df = df1, type = "linear", snp = snp, ncores = ncores)
    } else {
      gwas <- pvdiv_gwas(df = df1, type = "linear", snp = snp, covar = svd,
                         ncores = ncores, npcs = PCdf1$NumPCs)
    }

    # Save a data.table object with the GWAS results
    gwas_data <- data.table(CHR = markers$CHR, POS = markers$POS,
                            estim = gwas$estim, std_err = gwas$std.err,
                            bigsnpscore = gwas$score)
    gwas_data[,.data$pvalue := predict(gwas, log10 = FALSE)]
    gwas_data[,.data$log10p := -log10(.data$pvalue)]
    gwas_data[,.data$FDR_adj := p.adjust(.data$pvalue, method = "BH")]
    if(savegwas == TRUE){
      write_rds(gwas_data, path = paste0("GWAS_datatable_", names(df1)[2], "_",
                                         PCdf1$NumPCs, "_PCs", "_.rds"),
                compress = "gz")
    }

    if(saveplots == TRUE){
      # Find 10% FDR threshold
      FDRthreshhi <- gwas_data %>%
        as_tibble() %>%
        filter(between(.data$FDR_adj, 0.10001, 1)) %>%  # 10% FDR threshold currently.
        summarise(thresh = max(.data$log10p))
      FDRthreshlo <- gwas_data %>%
        as_tibble() %>%
        filter(between(.data$FDR_adj, 0, 0.09999)) %>%  # 10% FDR threshold currently.
        summarise(thresh = min(.data$log10p))
      if(FDRthreshhi$thresh[1] > 0 & FDRthreshlo$thresh[1] > 0){
        FDRthreshold = (FDRthreshhi$thresh[1] + FDRthreshlo$thresh[1])/2
      } else if(FDRthreshhi$thresh[1] > 0){
        FDRthreshold = FDRthreshhi$thresh[1]
      } else if(FDRthreshlo$thresh[1] > 0){
        FDRthreshold = FDRthreshlo$thresh[1]
      } else {
        FDRthreshold = NA
      }
      # Save a Manhattan plot with 10% FDR
      ggmanobject1 <- gwas_data %>%
        filter(.data$log10p > 1) %>%
        ggplot(aes(x = .data$POS, y = .data$log10p)) +
        geom_hline(yintercept = c(5, 10), color = "lightgrey") +
        geom_point(aes(color = .data$CHR, fill = .data$CHR)) +
        geom_hline(yintercept = FDRthreshold, color = "black", linetype = 2,
                   size = 1) +
        facet_wrap(~ .data$CHR, nrow = 1, scales = "free_x",
                   strip.position = "bottom") +
        scale_color_viridis(option = "B", end = 0.8, discrete = TRUE) +
        scale_fill_viridis(option = "B", end = 0.8, discrete = TRUE) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.background = element_rect(fill=NA), legend.position = "none") +
        labs(x = "Chromosome", y = "-log10(p value)") +
        scale_x_continuous(expand = c(0.18, 0.18))

      save_plot(paste0("Manhattan_", names(df1)[2], "_", PCdf1$NumPCs, "_PCs_",
                       "10percent_FDR_", get_date_filename(),
                       ".png"), plot = ggmanobject1, base_asp = 4, base_height = 4)

      # Save a QQplot
      ggqqplot <- pvdiv_qqplot(ps = gwas_data$pvalue, lambdaGC = TRUE)
      save_plot(paste0("QQplot_", names(df1)[2], "_", PCdf1$NumPCs, "_PCs_FDR_",
                       get_date_filename(), ".png"),
                plot = ggqqplot, base_asp = 1, base_height = 4)

      # Save a Manhattan plot with Bonferroni
      ggmanobject2 <- gwas_data %>%
        filter(.data$log10p > 1) %>%
        ggplot(aes(x = .data$POS, y = .data$log10p)) +
        geom_hline(yintercept = c(5, 10), color = "lightgrey") +
        geom_point(aes(color = .data$CHR, fill = .data$CHR)) +
        geom_hline(yintercept = bonferroni, color = "black", linetype = 2,
                   size = 1) +
        facet_wrap(~ .data$CHR, nrow = 1, scales = "free_x",
                   strip.position = "bottom") +
        scale_color_viridis(option = "B", end = 0.8, discrete = TRUE) +
        scale_fill_viridis(option = "B", end = 0.8, discrete = TRUE) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.background = element_rect(fill=NA), legend.position = "none") +
        labs(x = "Chromosome", y = "-log10(p value)") +
        scale_x_continuous(expand = c(0.18, 0.18))

      save_plot(paste0("Manhattan_", names(df1)[2], "_", PCdf1$NumPCs,
                       "_PCs_Bonferroni_", get_date_filename(),
                       ".png"), plot = ggmanobject2, base_asp = 4, base_height = 4)
    }

    if(saveannos == TRUE){
      ## Save annotation tables for the top associations
      anno_tables <- pvdiv_table_topsnps(df = gwas, type = "bigsnp", n = c(10,100),
                                         FDRalpha = 0.1,
                                         rangevector = c(0, 20000, 100000),
                                         markers = markers,
                                         anno_info = switchgrassGWAS::gff_gene,
                                         txdb = txdb)
      saveRDS(anno_tables, file.path("/", "home", "alice", "Github",
                                     "pvdiv-phenology-gxe", "analysis",
                                     "all-phenotypes-two-years",
                                     paste0("Annotation_tables_", names(df1)[2],
                                            "_", PCdf1$NumPCs, "_PCs", ".rds")))
    }
  }
}
