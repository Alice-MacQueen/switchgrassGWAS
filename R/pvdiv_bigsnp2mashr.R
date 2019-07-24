#' Identify phenotype names from bigsnpr results in a folder.
#'
#' @description Creates a vector of phenotype names from bigsnpr results.
#'
#' @param path File path to the files from bigsnpr, a character string.
#'     Defaults to the current working directory.
#' @param pattern Pattern within the filename to match. Default is "*.rds".
#'
#' @return A vector of phenotype names.
#'
#' @examples
#' \dontrun{gwas_results_in_folder(path = system.file("inst/extdata",
#' package = "switchgrassGWAS"))}
#' \dontrun{gwas_results_in_folder(path = "path/to/gwas/results")}
#'
#' @export
pvdiv_results_in_folder <- function(path = ".", pattern = "*.rds"){
  result_files <- list.files(path = path, pattern = pattern)
  return(result_files)
}


#' Step One of bigsnp2mashr
#'
#' @param path Path
#' @param gwas_rds RDS file with gwas results
#' @param phenotype Character vector. Single phenotype name
#' @param numSNPs Integer. Number of top SNPs to choose.
#' @param model One of linear or logistic. Type of GWAS model.
#' @param keep_score Logical. Keep p-values for clumping?
#' @param markers Marker names for the GWAS you ran
#'
#' @import bigsnpr
#' @importFrom stats predict
#' @importFrom dplyr left_join mutate top_n select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
pvdiv_top_effects_log10p <- function(path, gwas_rds, phenotype, numSNPs,
                                     model = c("linear", "logistic"),
                                     keep_score = FALSE, markers){
  gwas_obj <- readRDS(file = file.path(path, gwas_rds))
  log10p <- predict(gwas_obj)
  # Use bigsnpr predict function to generate p-values
  pre_mash_1 <- as_tibble(cbind(gwas_obj, markers, log10p))
  # Use SNP object and gwas_obj generated using SNP object to make a full
  # pre_mash object for one phenotype.
  if(keep_score == FALSE){
    top_set <- pre_mash_1 %>%
      top_n(-numSNPs, .data$log10p) %>%
      dplyr::select(.data$CHR, .data$POS)
  } else {
    top_set <- pre_mash_1 %>%
      top_n(-numSNPs, .data$log10p) %>%
      mutate(score_log10p = -log10p) %>%
      dplyr::select(.data$CHR, .data$POS, .data$score_log10p)
    names(top_set)[3] <- paste0(phenotype, "_score_log10p")
  }
  rm(pre_mash_1)
  rm(gwas_obj)
  rm(log10p)
  # optionally save this pre_mash_1 for reuse in step 2, if it's faster to
  # read/write this than to generate a new pre_mash_1 in both steps. Both are
  # time consuming so I can compare their speeds during testing.

  # name top_set columns appropriately so that the joined file will have
  # informative column names.
  return(top_set)
}

#' Step Two of bigsnp2mashr
#'
#' @param path Path
#' @param gwas_rds RDS file with gwas results
#' @param phenotype Character vector. Single phenotype name
#' @param top_set Top markers chosen
#' @param random_sample Numeric vector. Random sample of SNPs
#' @param markers Marker names for the GWAS you ran
#' @param markers2 Markers to include if SNPs are LD clumped.
#' @param model One of linear or logistic. Type of GWAS model.
#' @param clumped Logical. Clump SNPs?
#'
#' @import bigsnpr
#' @importFrom stats predict
#' @importFrom dplyr left_join mutate
s_hat_bigsnp <- function(path, gwas_rds, phenotype, top_set, random_sample,
                         markers = switchgrassGWAS::markers,
                         markers2 = switchgrassGWAS::markers2,
                         model = c("linear", "logistic"), clumped = FALSE){
  gwas_obj <- readRDS(file = file.path(path, gwas_rds))
  log10p <- predict(gwas_obj)
  # Use bigsnpr predict function to generate p-values
  pre_mash_1 <- as_tibble(cbind(gwas_obj, markers, log10p))
  standardization <- max(abs(pre_mash_1$estim), na.rm = TRUE)

  pre_mash_strong <- top_set %>%
    dplyr::left_join(pre_mash_1, by = c("CHR", "POS")) %>%
    dplyr::mutate(stderr_d = .data$std.err / standardization,
                  effect_d = .data$estim / standardization,
                  stderr_d = ifelse(is.na(.data$stderr_d),
                                    10,
                                    .data$stderr_d),
                  effect_d = ifelse(is.na(.data$effect_d),
                                    0,
                                    .data$effect_d))
  if(clumped == FALSE){
    pre_mash_random <- pre_mash_1[random_sample,] %>%
      dplyr::mutate(stderr_d = .data$std.err / standardization,
                    effect_d = .data$estim / standardization,
                    stderr_d = ifelse(is.na(.data$stderr_d),
                                      10,
                                      .data$stderr_d),
                    effect_d = ifelse(is.na(.data$effect_d),
                                      0,
                                      .data$effect_d))
  } else {
    pre_mash_random <- markers2 %>%
      left_join(pre_mash_1, by = c("CHR", "POS"))
    pre_mash_random <- pre_mash_random[random_sample,] %>%
      dplyr::mutate(stderr_d = .data$std.err / standardization,
                    effect_d = .data$estim / standardization,
                    stderr_d = ifelse(is.na(.data$stderr_d),
                                      10,
                                      .data$stderr_d),
                    effect_d = ifelse(is.na(.data$effect_d),
                                      0,
                                      .data$effect_d))
  }

  if(model == "linear"){
    # CHR POS max_score_log10p CHRN estim std.err score log10p stderr_d effect_d
    # CHR POS estim std.err tscore log10p stderr_d effect_d

    names(pre_mash_strong)[3] <- paste0(phenotype, "_maxscore")
    names(pre_mash_strong)[4] <- paste0("Bhat_", phenotype)
    names(pre_mash_strong)[5] <- paste0("Shat_",phenotype)
    names(pre_mash_strong)[6] <- paste0(phenotype, "_tscore")
    names(pre_mash_strong)[7] <- paste0(phenotype, "_log10p")
    names(pre_mash_strong)[8] <- paste0("Stand_Shat", phenotype)
    names(pre_mash_strong)[9] <- paste0("Stand_Bhat", phenotype)

    names(pre_mash_random)[3] <- paste0("Bhat_", phenotype)
    names(pre_mash_random)[4] <- paste0("Shat_",phenotype)
    names(pre_mash_random)[5] <- paste0(phenotype, "_tscore")
    names(pre_mash_random)[6] <- paste0(phenotype, "_log10p")
    names(pre_mash_random)[7] <- paste0("Stand_Shat", phenotype)
    names(pre_mash_random)[8] <- paste0("Stand_Bhat", phenotype)
  }
  if(model == "logistic"){
    # CHR POS max_score_log10p CHRN estim std.err score log10p stderr_d effect_d
    # CHR POS estim std.err niter zscore log10p stderr_d effect_d
    names(pre_mash_strong)[3] <- paste0(phenotype, "_maxscore")
    names(pre_mash_strong)[4] <- paste0("Bhat_", phenotype)
    names(pre_mash_strong)[5] <- paste0("Shat_",phenotype)
    names(pre_mash_strong)[6] <- paste0(phenotype, "_niter")
    names(pre_mash_strong)[7] <- paste0(phenotype, "_zscore")
    names(pre_mash_strong)[8] <- paste0(phenotype, "_log10p")
    names(pre_mash_strong)[9] <- paste0("Stand_Shat_",phenotype)
    names(pre_mash_strong)[10] <- paste0("Stand_Bhat_",phenotype)

    names(pre_mash_random)[3] <- paste0("Bhat_", phenotype)
    names(pre_mash_random)[4] <- paste0("Shat_",phenotype)
    names(pre_mash_random)[5] <- paste0(phenotype, "_niter")
    names(pre_mash_random)[6] <- paste0(phenotype, "_zscore")
    names(pre_mash_random)[7] <- paste0(phenotype, "_log10p")
    names(pre_mash_random)[8] <- paste0("Stand_Shat_",phenotype)
    names(pre_mash_random)[9] <- paste0("Stand_Bhat_",phenotype)
  }
  # name top_set columns appropriately so that the joined file will have
  # informative column names.

  return(list(strong_df_1 = pre_mash_strong, random_df_1 = pre_mash_random))
}


#' Convert bigsnpr output to mashr input dataframes.
#'
#' @description This function converts bigsnpr output, saved as rds files to
#'    the specified path, to four dataframes used in the R package mashr. It
#'    can clump SNPs based on LD and the maximum -log10(p-value) across all
#'    included GWAS. It can also set the random effect data frames to come from
#'    a subsample of SNPs clumped by MAF and LD.
#'
#' @param path File path to the rds files saved from bigsnpr, a character
#'     string. Defaults to the working directory.
#' @param gwas_rds A character vector of saved gwas rds objects from bigsnpr. If NA, all *.rds files in the path will be used.
#' @param phenotypes A character vector of phenotype names for the GWAS RDS
#'    objects. Must be the same length as gwas_rds, or NA. If NA, these will be
#'    the rds file names.
#' @param clumping Logical. Should SNPs be clumped by LD & p-value to
#'    standardize signal strength across different LD blocks? Default is TRUE.
#' @param standardization Logical. Should marker effects in each condition be
#'    standardized to fall between -1 and 1? Default is TRUE.
#' @param G snp$genotypes from a FBS object read in by \code{snp_attach().}
#' @param markers A dataframe containing two columns: CHR & POS. POS should be
#'    an integer vector with marker positions. CHR contains marker chromosomes.
#' @param markers2 A dataframe of SNPs clumped by LD & MAF containing two
#'    columns, as for markers. Used when clumping = TRUE.
#' @param numSNPs The number of most significant SNPs selected from each GWAS.
#'     Ideally this will give 1 million or fewer total cells in the resultant
#'     mash dataframes. Defaults to 1000.
#' @param model Regression used in bigstatsr. One of "logistic" or "linear".
#'     Default is "linear".
#' @param saveoutput Logical. Should the function's output also be saved to RDS
#' files? Default is FALSE.
#'
#' @return A list containing five data frames: the SNPs selected, the B_hat
#'    and S_hat matrices for the strong SNP set and for a random SNP set that
#'    is twice the size.
#'
#' @note To create a vector of phenotype names, use the
#'     \code{\link{pvdiv_results_in_folder}} function.
#'
#' @examples
#' \dontrun{bigsnp2mashr(path = system.file("inst/extdata"), numSNPs = 20,
#'     model = "linear")}
#' \dontrun{bigsnp2mashr(numSNPs = 10000, model = "logistic")}
#' \dontrun{bisgnp2mashr(numSNPs = 20000, model = "linear", saveoutput = TRUE)}
#' \dontrun{phenotype_vector <- pvdiv_results_in_folder(path = system.file(
#'     "inst/extdata"))
#'     numSNPs <- 1000000 / length(phenotype_vector)^2
#'     bigsnp2mashr(phenotypes = phenotype_vector, numSNPs = numSNPs,
#' model = "linear", saveoutput = TRUE)}
#'
#' @import bigsnpr
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stringr str_sub
#' @import tidyr
#' @importFrom tidyselect starts_with
#'
#' @export
pvdiv_bigsnp2mashr <- function(path = ".", gwas_rds = NA, phenotypes = NA,
                               clumping = TRUE, standardization = TRUE, G,
                               markers = NULL,
                               markers2 = NULL,
                               numSNPs = 1000, model = c("linear", "logistic"),
                               saveoutput = FALSE){
  match.arg(model, c("linear", "logistic"))
  if(is.na(gwas_rds)[1]){
    phe_col <- pvdiv_results_in_folder(path = path)
  } else {
    phe_col <- gwas_rds
  }
  if(is.null(phe_col)){
    stop("Can't find any bigsnp rds files in this path.")
  }
  if(is.na(phe_col[1])){
    stop("Can't find any bigsnp rds files in this path.")
  }
  if(is.na(phenotypes)[1]){
    phenotypes <- str_sub(phe_col, end = -5)
  }
  numSNPs <- as.integer(numSNPs)
  #load("data/markers.rda")
  message(paste0("Starting part one: Making a data frame of all SNPs that are",
                 " in the top ", numSNPs, " SNPs
                 by maximum -log10(p-values) for at least one phenotype."))

  top_set <- pvdiv_top_effects_log10p(path = path, gwas_rds = phe_col[1],
                                      phenotype = phenotypes[1],
                                      numSNPs = numSNPs, model = model,
                                      keep_score = TRUE, markers = markers)

  for(i in seq_along(phe_col)[-1]){
    message(paste0("Picking top SNPs for GWAS for ", phenotypes[i]," (", i,
                   " of ", length(phe_col), ")."))
    top_set_new <- pvdiv_top_effects_log10p(path = path, gwas_rds = phe_col[i],
                                            phenotype = phenotypes[i],
                                            numSNPs = numSNPs, model = model,
                                            keep_score = TRUE, markers = markers)
    top_set <- top_set %>%
      full_join(top_set_new, by = c("CHR", "POS"))
  }
  top_set <- top_set %>%
    dplyr::arrange(.data$CHR, .data$POS) %>%
    gather(key = "Condition", value = "score_log10p", -(1:2)) %>%
    group_by(.data$CHR, .data$POS) %>%
    filter(!is.na(.data$score_log10p)) %>%
    mutate(max_score_log10p = max(.data$score_log10p)*n()) %>%
    slice(which.max(.data$max_score_log10p)) %>%
    dplyr::select(-.data$Condition, -.data$score_log10p) %>%
    arrange(.data$CHR, .data$POS) %>%
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
  if(clumping == TRUE){
    top_clumps <- snp_clumping(G, infos.chr = top_set$CHRN,
                               S = top_set$max_score_log10p,
                               infos.pos = top_set$POS) # not sure if this will work
    top_set <- top_set[top_clumps,] %>%
      dplyr::select(.data$CHR, .data$POS, .data$max_score_log10p)
  } else {
    top_set <- top_set %>%
      dplyr::select(.data$CHR, .data$POS)
  }
  if(saveoutput == TRUE){
    saveRDS(top_set, file = file.path(path, paste0("Part-One-Output_",
                                                   "Top-Effects-", numSNPs,
                                                   "-SNPs.rds")))
  }

  message(paste0("Part One: data frame of SNPs to keep complete.
                 ----------------"))
  message(paste0("Starting Part Two: Creating strong and random dataframes of ",
                 "B_hat and S_hat
                 values for use in mashr."))

  set.seed(1234) # Makes the random data frames reproducible.
  if(clumping == TRUE){
    random_sample <- sample(1:nrow(markers2), length(top_clumps)*10) %>%
      sort()
  } else {
    random_sample <- sample(1:nrow(markers), nrow(top_set)) %>%
      sort()
  }
  mash_list_1 <- s_hat_bigsnp(path = path, gwas_rds = phe_col[1],
                              phenotype = phenotypes[1],
                              top_set = top_set,
                              random_sample = random_sample, model = model,
                              markers = markers, markers2 = markers2,
                              clumped = clumping)

  mash_df_strong <- mash_list_1$strong_df_1
  mash_df_random <- mash_list_1$random_df_1

  for(i in seq_along(phe_col)[-1]){
    message(paste0("Determining standardized B_hat and S_hat values for ",
                   nrow(top_clumps), " SNPs in GWAS for ", phenotypes[i]," (",
                   i, " of ", length(phe_col), ")."))
    mash_list_1 <- s_hat_bigsnp(path = path, gwas_rds = phe_col[i],
                                phenotype = phenotypes[i],
                                top_set = top_set,
                                random_sample = random_sample, model = model,
                                markers = markers, markers2 = markers2,
                                clumped = clumping)

    mash_df_strong <- mash_df_strong %>%
      dplyr::left_join(mash_list_1$strong_df_1, by = c("CHR", "POS"))
    mash_df_random <- mash_df_random %>%
      dplyr::left_join(mash_list_1$random_df_1, by = c("CHR", "POS"))
  }
  if(standardization == TRUE){
    bhat_strong <- mash_df_strong %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Stand_Bhat"))
    shat_strong <- mash_df_strong %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Stand_Shat"))
    bhat_random <- mash_df_random %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Stand_Bhat"))
    shat_random <- mash_df_random %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Stand_Shat"))
  } else {
    bhat_strong <- mash_df_strong %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Bhat"))
    shat_strong <- mash_df_strong %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Shat"))
    bhat_random <- mash_df_random %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Bhat"))
    shat_random <- mash_df_random %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Shat"))
    # make a SNP column by uniting CHR and POS
    # make Bhat only and Shat only tables for both strong and random sets.
  }
  B_hat_random <- data.frame(bhat_random, row.names = "SNP")
  S_hat_random <- data.frame(shat_random, row.names = "SNP")
  B_hat_strong <- data.frame(bhat_strong, row.names = "SNP")
  S_hat_strong <- data.frame(shat_strong, row.names = "SNP")

  if(saveoutput == TRUE){
    saveRDS(B_hat_strong, file = file.path(path, paste0("B_hat_strong_df_",
                                                        numSNPs, "topSNPs.rds")))
    saveRDS(S_hat_strong, file = file.path(path, paste0("S_hat_strong_df_",
                                                        numSNPs, "topSNPs.rds")))
    saveRDS(B_hat_random, file = file.path(path, paste0("B_hat_random_df_",
                                                        numSNPs, "topSNPs.rds")))
    saveRDS(S_hat_random, file = file.path(path, paste0("S_hat_random_df_",
                                                        numSNPs, "topSNPs.rds")))
  }
  return(list(top_set = top_set,
              random_sample = random_sample,
              B_hat_strong = B_hat_strong,
              S_hat_strong = S_hat_strong,
              B_hat_random = B_hat_random,
              S_hat_random = S_hat_random))
}
