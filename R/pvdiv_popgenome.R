#' Format the VCF subset
#'
#' @param path file.path() or character string. Path to the VCF file subset.
#' @param filename Character string. VCF subset filename
#'
#' @importFrom readr read_delim write_delim
#' @importFrom dplyr mutate case_when
#' @importFrom stringr str_sub
#' @importFrom magrittr %>%
#'
pvdiv_vcf_popgenome <- function(path, filename){
  vcf <- read_delim(file = file.path(path, filename), col_names = FALSE,
                    skip = 5, delim = "\t")
  vcf <- vcf %>%
    mutate(X1 = case_when(.data$X1 == "Chr01K" ~ 1,
                          .data$X1 == "Chr01N" ~ 2,
                          .data$X1 == "Chr02K" ~ 3,
                          .data$X1 == "Chr02N" ~ 4,
                          .data$X1 == "Chr03K" ~ 5,
                          .data$X1 == "Chr03N" ~ 6,
                          .data$X1 == "Chr04K" ~ 7,
                          .data$X1 == "Chr04N" ~ 8,
                          .data$X1 == "Chr05K" ~ 9,
                          .data$X1 == "Chr05N" ~ 10,
                          .data$X1 == "Chr06K" ~ 11,
                          .data$X1 == "Chr06N" ~ 12,
                          .data$X1 == "Chr07K" ~ 13,
                          .data$X1 == "Chr07N" ~ 14,
                          .data$X1 == "Chr08K" ~ 15,
                          .data$X1 == "Chr08N" ~ 16,
                          .data$X1 == "Chr09K" ~ 17,
                          .data$X1 == "Chr09N" ~ 18,
                          TRUE ~ 19
    ))

  vcfmeta <- read_delim(file = file.path(path, filename), n_max = 5,
                        comment = "@", col_names = FALSE, delim = ";")

  outputfile <- str_sub(filename, end = -5) %>%
    paste0(.data, "_num.vcf")
  dir.create(path = file.path(path, "VCF"), showWarnings = FALSE)

  write_delim(vcfmeta, path = file.path(path, "VCF", outputfile), delim = ";",
              col_names = FALSE)
  write_delim(vcf, path = file.path(path, "VCF", outputfile), delim = "\t",
              col_names = FALSE, append = TRUE)
  return(outputfile)
}

#' Get GFF subsets
#'
#' @description Prepare a subset of the P. virgatum VCF file to be used in the
#'     PopGenome R package. Given a VCF subset (created using VCFtools or some
#'     other commandline program), this function will make a GFF and VCF folder
#'     in the path, if they do not already exist, and place the VCF subset in
#'     the correct format in the VCF folder, and the GFF subset in the correct
#'     format in the GFF folder. These files can then be read into the R
#'     package PopGenome.
#'
#' @param path file.path() or character string. Path to the VCF file subset.
#' @param filename Character string. VCF subset filename
#' @param seqid Character string. name of the chromosome or scaffold; here, of Panicum virgatum.
#' @param start Integer. Start position of the feature, with sequence numbering starting at 1.
#' @param end Integer. End position of the feature, with sequence numbering starting at 1.
#'
#' @importFrom readr write_delim
#' @importFrom dplyr mutate case_when filter
#' @importFrom stringr str_sub
#' @importFrom magrittr %>%
#'
#' @export
pvdiv_subset_popgenome <- function(path, filename, seqid, start, end){

  outputfile <- pvdiv_vcf_popgenome(path = path, filename = filename)
  dir.create(path = file.path(path, "GFF"), showWarnings = FALSE)

  gff_sub <- switchgrassGWAS::gff_gene %>%
    filter(.data$seqid == seqid & .data$start >= start & .data$end >= start &
             .data$start <= end & .data$end <= end) %>%
    mutate(seqid = case_when(.data$seqid == "Chr01K" ~ 1,
                             .data$seqid == "Chr01N" ~ 2,
                             .data$seqid == "Chr02K" ~ 3,
                             .data$seqid == "Chr02N" ~ 4,
                             .data$seqid == "Chr03K" ~ 5,
                             .data$seqid == "Chr03N" ~ 6,
                             .data$seqid == "Chr04K" ~ 7,
                             .data$seqid == "Chr04N" ~ 8,
                             .data$seqid == "Chr05K" ~ 9,
                             .data$seqid == "Chr05N" ~ 10,
                             .data$seqid == "Chr06K" ~ 11,
                             .data$seqid == "Chr06N" ~ 12,
                             .data$seqid == "Chr07K" ~ 13,
                             .data$seqid == "Chr07N" ~ 14,
                             .data$seqid == "Chr08K" ~ 15,
                             .data$seqid == "Chr08N" ~ 16,
                             .data$seqid == "Chr09K" ~ 17,
                             .data$seqid == "Chr09N" ~ 18,
                             TRUE ~ 19
    ))

  write_delim(switchgrassGWAS::gff_metadata, path = file.path(path, "GFF", outputfile),
              delim = "\t", col_names = FALSE)
  write_delim(gff_sub, path = file.path(path, "GFF", outputfile), delim = "\t",
              col_names = FALSE, append = TRUE)
}
