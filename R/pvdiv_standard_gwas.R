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
#' snpfile <- system.file("extdata", "example_bigsnp.rds", package = "switchgrassGWAS")
#' library(bigsnpr)
#' snp <- snp_attach(snpfile)
#' svd5 <- pvdiv_autoSVD(snp = snp, k = 5, saveoutput = FALSE)
#'
#' @export
pvdiv_autoSVD <- function(snp, k = 10, ncores = 1, saveoutput = FALSE, ...){
  requireNamespace("dots")
  fun.scaling <- dots::dots(name = 'fun.scaling', value = snp_scaleBinom(),
                            ...)
  thr.r2 <- dots::dots(name = 'thr.r2', value = 0.2, ...)
  size <- dots::dots(name = 'size', value = 100/thr.r2, ...)
  roll.size <- dots::dots(name = 'roll.size', value = 50, ...)
  int.min.size <- dots::dots(name = 'int.min.size', value = 20, ...)
  #alpha.tukey <- dots::dots(name = 'alpha.tukey', value = 0.05, ...)
  #min.mac <- dots::dots(name = 'min.mac', value = 10, ...)
  #max.iter <- dots::dots(name = 'max.iter', value = 5, ...)
  verbose <- dots::dots(name = 'verbose', value = FALSE, ...)

  # Argument error checks; Set up numeric chromosome data frame.
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
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
    message(paste0("'saveoutput' is FALSE, so the svd will not be saved to ",
                   "the working directory."))
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
#' @param hasInbred Logical. Does the SNP file contain inbred individuals or
#'     closely related individuals, like siblings?
#' @param saveoutput Logical. Should the output be saved to the working
#'     directory?
#'
#' @return A kinship matrix with labeled rows and columns.
#'
#' @import bigsnpr
#' @import bigstatsr
#'
#' @examples
#' snpfile <- system.file("extdata", "example_bigsnp.rds", package = "switchgrassGWAS")
#' library(bigsnpr)
#' snp <- snp_attach(snpfile)
#' K <- pvdiv_kinship(snp = snp, saveoutput = FALSE)
#'
#' @export
pvdiv_kinship <- function(snp, ind.row = NA, hasInbred = TRUE,
                          saveoutput = FALSE){
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(saveoutput == FALSE){
    message(paste0("'saveoutput' is FALSE, so the kinship matrix will not ",
                   "be saved to the working directory."))
  }

  G <- snp$genotypes
  if(!is.na(ind.row[1])){
    nInd <- length(ind.row)
  } else {
    nInd <- snp$genotypes$nrow
    ind.row <- 1:nInd
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
  #Handler of inbred
  if(MD < 2 & hasInbred){
    K = 2*K / ((DU - floor) / (DL - floor))
  }

  if(saveoutput){
    saveRDS(K, paste0("Kinship_van_Raden_method_", nInd, "_individuals_",
    ".rds"))
  }
  return(K)
}

#' @title Juenger lab standard GWAS function.
#'
#' @description This function is a wrapper around the standard GWAS procedures
#'     in the Juenger lab. Singular value decomposition of the SNPs is done to
#'     get principal components for population structure correction; the 'best'
#'     number of PCs is chosen as the one that makes lambda_GC, the Genomic
#'     Control coefficient, closest to 1. (See the `lambdagc` parameter to set
#'     this yourself.) Next, genome-wide association is conducted, and the GWAS
#'     output can be saved, as well as Manhattan plots, QQ-plots, and annotation
#'     information for the top SNPs for each phenotype.
#'
#' @param snp A "bigSNP" object; load with \code{bigsnpr::snp_attach()}.
#'     Here, genomic information for Panicum virgatum. Contact tjuenger <at>
#'     utexas <dot> edu to obtain this information pre-publication.
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
#' @param minphe Integer. What's the minimum number of phenotyped individuals
#'     to conduct a GWAS on? Default is 200. Use lower values with caution.
#' @param ... Other arguments to \code{\link{pvdiv_lambda_GC}} or
#'     \code{\link{pvdiv_table_topsnps}}.
#'
#' @return A big_SVD object.
#'
#' @import bigsnpr
#' @import bigstatsr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when
#' @importFrom tibble enframe as_tibble tibble
#' @importFrom rlang .data
#' @importFrom readr write_rds
#' @importFrom stats p.adjust
#' @importFrom tidyselect all_of
#' @importFrom cowplot save_plot
#'
#' @examples
#' snpfile <- system.file("extdata", "example_bigsnp.rds", package = "switchgrassGWAS")
#' library(bigsnpr)
#' snp <- snp_attach(snpfile)
#' pvdiv_standard_gwas(snp, df = phenotypes, type = "linear", savegwas = FALSE,
#'     saveplots = FALSE, ncores = 1)
#'
#' \dontrun{
#' # Here we specify that we do want to generate and save the gwas dataframes,
#' # the Manhattan and QQ-plots, and the annotation tables.
#' pvdiv_standard_gwas(snp, df = phenotypes, type = "linear", covar = svd,
#'     ncores = nb_cores(), lambdagc = TRUE, savegwas = TRUE, saveplots = TRUE,
#'     saveannos = TRUE, txdb = txdb)
#' }
#'
#' @export
pvdiv_standard_gwas <- function(snp, df = switchgrassGWAS::phenotypes,
                                type = c("linear", "logistic"),
                                ncores = nb_cores(),
                                outputdir = ".", covar = NULL, lambdagc = TRUE,
                                savegwas = FALSE, saveplots = TRUE,
                                saveannos = FALSE, txdb = NULL, minphe = 200,
                                ...){

  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
    }
  if(colnames(df)[1] != "PLANT_ID"){
    stop("First column of phenotype dataframe (df) must be 'PLANT_ID'.")
    }
  stopifnot(type %in% c("linear", "logistic"))
  if(saveannos == TRUE & is.null(txdb)){
    stop(paste0("Need to specify a txdb object created using AnnotationDbi ",
                "in order to generate data frames containing annotated top ",
                "SNPs. If you don't have this, set saveannos = FALSE."))
    }
  if(is.null(colnames(lambdagc))){
    message(paste0("'lambdagc' is TRUE, so lambda_GC will be used to find ",
                   "the best population structure correction using the ",
                   "covariance matrix."))
    } else if(!(colnames(lambdagc)[1] == "NumPCs" &
                colnames(lambdagc)[2] %in% colnames(df))){
    stop(paste0("If lambdagc is a dataframe, the column names must include",
                " NumPCs and the names of the phenotypes to run GWAS on ",
                "(the names of 'df'). You can generate this data frame with ",
                "pvdiv_lambda_GC()."))
      }
  if(savegwas == FALSE){
    message(paste0("'savegwas' is FALSE, so the gwas results will not be ",
                   "saved to disk."))
    }

  nSNP_M <- round(snp$genotypes$ncol/1000000, digits = 1)
  nInd <- snp$genotypes$nrow
  if(is.null(covar)){
    message(paste0("Covariance matrix (covar) was not supplied - this will be",
                   " generated using pvdiv_autoSVD()."))
    requireNamespace("dots")
    k <- dots::dots(name = 'k', value = 15, ...)
    covar <- pvdiv_autoSVD(snp, k = k, ncores = ncores, saveoutput = FALSE)
    if(savegwas == TRUE){
      saveRDS(covar, file = file.path(outputdir, paste0("SVD_", nInd, "g_",
                                                        nSNP_M, "M_SNPs_",
                                                        k, "PCs.rds")))
      }
    } else {
      stopifnot(attr(covar, "class") == "big_SVD")
      }

  plants <- snp$fam$sample.ID
  bonferroni <- -log10(0.05/length(snp$map$physical.pos))
  markers <- tibble(CHR = snp$map$chromosome, POS = snp$map$physical.pos)
  df <- plants %>%
    enframe(name = NULL, value = "PLANT_ID") %>%
    left_join(df, by = "PLANT_ID")

  for(i in 2:ncol(df)){

    df1 <- df %>%
      dplyr::select(.data$PLANT_ID, all_of(i))
    phename <- names(df1)[2]
    nPhe <- length(which(!is.na(df1[,2])))
    nLev <- nrow(unique(df1[which(!is.na(df1[,2])),2]))
    # Checks for correct combinations of phenotypes and GWAS types.
    if(nPhe < minphe){
      message(paste0("The phenotype ", phename, " does not have the minimum ",
                     "number of phenotyped PLANT_ID's, (", minphe, ") and so ",
                     "will not be used for GWAS."))
      next
    } else if(nLev < 2){
      message(paste0("The phenotype ", phename, " does not have two or more ",
                     "distinct non-NA values and will not be used for GWAS."))
      next
    } else if(nLev > 2 & type == "logistic"){
      message(paste0("The phenotype ", phename, " has more than two distinct ",
                     "non-NA values and will not be used for GWAS with 'type=",
                     "logistic'."))
      next
    } else if(!(unique(df1[which(!is.na(df1[,2])),2])[1,1] %in% c(0,1)) &
              !(unique(df1[which(!is.na(df1[,2])),2])[2,1] %in% c(0,1)) & type == "logistic"){
      message(paste0("The phenotype ", phename, " has non-NA values that are ",
                     "not 0 or 1 and will not be used for GWAS with 'type=",
                     "logistic'."))
      next
    } else {
    message(paste0("Now starting GWAS pipeline for ", phename, "."))

      if(nLev == 2 & type == "linear"){
        message(paste0("The phenotype ", phename, " has only two distinct non-",
                       "NA values; consider using a logistic model instead.",
                       "(Set type = 'logistic')."))
      }

    if(is.null(colnames(lambdagc))){
      pc_max = ncol(covar$u)
      message(paste0("Now determining lambda_GC for GWAS models with ",
                     pc_max+1, " sets of PCs. This will take some time."))
      lambdagc_df <- pvdiv_lambda_GC(df = df1, type = type, snp = snp,
                                  covar = covar, ncores = ncores,
                                  npcs = c(0:pc_max), saveoutput = FALSE)
      if(saveplots == TRUE){
        write_csv(lambdagc_df, path = file.path(outputdir,
                                                paste0("Lambda_GCs_by_PC_Num_",
                                                       phename, "_", type,
                                                       "_model_", nPhe, "g_",
                                                       nSNP_M, "M_SNPs",
                                                       ".csv")))
      }
      PCdf <- get_best_PC_df(lambdagc_df) # asv_best_PC_df(lambdagc_df)
    } else {
      PC1 <- lambdagc %>%
        dplyr::select(.data$NumPCs, phename)
      PCdf <- get_best_PC_df(PC1) # asv_best_PC_df(PC1)
    }
    PCdf1 <- PCdf[1,]

  # ------ Run GWAS with best pop structure correction -----
  message("Now running GWAS with the best population structure correction.")
    if(PCdf1$NumPCs == 0){
      gwas <- pvdiv_gwas(df = df1, type = type, snp = snp, ncores = ncores)
      } else {
      gwas <- pvdiv_gwas(df = df1, type = type, snp = snp, covar = covar,
                         ncores = ncores, npcs = PCdf1$NumPCs)
      }

    gwas_data <- tibble(CHR = markers$CHR, POS = markers$POS,
                        estim = gwas$estim, std_err = gwas$std.err,
                        bigsnpscore = gwas$score,
                        pvalue = predict(gwas, log10 = FALSE))
    gwas_data <- gwas_data %>%
      mutate(log10p = -log10(.data$pvalue))
    gwas_data$FDR_adj <- p.adjust(gwas_data$pvalue, method = "BH")
    if(savegwas == TRUE){
      # Save a data.table object with the GWAS results
      write_rds(gwas_data, file = file.path(outputdir,
                                            paste0("GWAS_datatable_", phename,
                                                   "_", type, "_model_", nPhe,
                                                   "g_", nSNP_M, "M_SNPs_",
                                                   PCdf1$NumPCs, "_PCs_",
                                                   "_.rds")), compress = "gz")
      }
    if(saveplots == TRUE){
      message("Now generating and saving Manhattan and QQ plots.")
      # ggplot settings for Manhattans and QQ plots.
      theme_oeco <- theme_classic() +
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
              strip.placement = 'outside', panel.spacing.x = unit(-0.5, 'cm'))
      # Find 10% FDR threshold
      FDRthreshhi <- gwas_data %>%
        as_tibble() %>%
        filter(between(.data$FDR_adj, 0.10001, 1)) %>%  # 10% FDR threshold
        summarise(thresh = max(.data$log10p))
      FDRthreshlo <- gwas_data %>%
        as_tibble() %>%
        filter(between(.data$FDR_adj, 0, 0.09999)) %>%  # 10% FDR threshold
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
        theme_oeco +
        geom_hline(yintercept = c(5, 10), color = "lightgrey") +
        geom_point(aes(color = .data$CHR, fill = .data$CHR)) +
        geom_hline(yintercept = FDRthreshold, color = "black", linetype = 2,
                   size = 1) +
        facet_wrap(~ .data$CHR, nrow = 1, scales = "free_x",
                   strip.position = "bottom") +
        scale_color_manual(values = rep(c("#1B0C42FF", "grey"), 9),
                           guide = FALSE) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.background = element_rect(fill=NA),
              legend.position = "none") +
        labs(x = "Chromosome", y = "-log10(p value)") +
        scale_x_continuous(expand = c(0.18, 0.18))

      save_plot(filename = file.path(outputdir,
                                     paste0("Manhattan_", phename, "_", type,
                                            "_model_", nPhe, "g_", nSNP_M,
                                            "M_SNPs_", PCdf1$NumPCs,
                                            "_PCs_10percent_FDR_",
                                            get_date_filename(), ".png")),
                plot = ggmanobject1, base_asp = 4, base_height = 4)

      # Save a QQplot
      ggqqplot <- pvdiv_qqplot(ps = gwas_data$pvalue, lambdaGC = TRUE)
      save_plot(filename = file.path(outputdir,
                                     paste0("QQplot_", phename, "_", type,
                                            "_model_", nPhe, "g_", nSNP_M,
                                            "M_SNPs_", PCdf1$NumPCs,
                                            "_PCs_FDR_",
                                            get_date_filename(), ".png")),
                plot = ggqqplot + theme_oeco, base_asp = 1, base_height = 4)

      # Save a Manhattan plot with Bonferroni
      ggmanobject2 <- gwas_data %>%
        filter(.data$log10p > 1) %>%
        ggplot(aes(x = .data$POS, y = .data$log10p)) +
        theme_oeco +
        geom_hline(yintercept = c(5, 10), color = "lightgrey") +
        geom_point(aes(color = .data$CHR, fill = .data$CHR)) +
        geom_hline(yintercept = bonferroni, color = "black", linetype = 2,
                   size = 1) +
        facet_wrap(~ .data$CHR, nrow = 1, scales = "free_x",
                   strip.position = "bottom") +
        scale_color_manual(values = rep(c("#1B0C42FF", "grey"), 9),
                           guide = FALSE) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.background = element_rect(fill=NA),
              legend.position = "none") +
        labs(x = "Chromosome", y = "-log10(p value)") +
        scale_x_continuous(expand = c(0.18, 0.18))

      save_plot(filename = file.path(outputdir,
                                     paste0("Manhattan_", phename, "_", type,
                                            "_model_", nPhe, "g_", nSNP_M,
                                            "M_SNPs_", PCdf1$NumPCs,
                                            "_PCs_Bonferroni_",
                                            get_date_filename(), ".png")),
                plot = ggmanobject2, base_asp = 4, base_height = 4)
      }

    if(saveannos == TRUE){
      message(paste0("Now creating annotation data frames for the top 10 & ",
                     "top 500 SNPs by p-value, and for SNPs above a 10% FDR."))
      ## Save annotation tables for the top associations
      requireNamespace("dots")
      n <- dots::dots(name = 'n', value = c(10, 500), ...)
      FDRalpha <- dots::dots(name = 'FDRalpha', value = 0.1, ...)
      rangevector <- dots::dots(name = 'rangevector', value = c(0, 50000), ...)
      anno_tables <- pvdiv_table_topsnps(df = gwas, type = "bigsnp",
                                         n = n, FDRalpha = FDRalpha,
                                         rangevector = rangevector,
                                         snp = snp, txdb = txdb)
      saveRDS(anno_tables, file.path(outputdir,
                                     paste0("Annotation_tables_", phename,
                                            "_", type, "_model_", nPhe, "g_",
                                            nSNP_M, "M_SNPs_",PCdf1$NumPCs,
                                            "_PCs", ".rds")))
      }
    }
    }
  }
