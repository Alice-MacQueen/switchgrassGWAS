#' Phenotypes available for the Panicum virgatum diversity panel.
#'
#' A dataset containing phenotypes for the Panicum virgatum diversity panel.
#'     The variables are as follows:
#'
#' \itemize{
#'   \item PLANT_ID. ID of the Panicum virgatum genotype sequenced in the diversity panel.
#'   \item GWAS_CT. The number of times this genotype was cloned and planted in the ten 2018 sites for the Panicum virgatum diversity panel.
#'   \item BRKG_TC_EOS_2018.
#'   \item BRKG_DEAD_2018.
#' }
#'
#' @name phenotypes
#' @docType data
#' @author Alice MacQueen \email{alice.macqueen@@utexas.edu}
#' @references \url{data_blah.com}
#' @keywords data
#' @usage data(phenotypes)
#' @format A data frame with 785 rows and 4 variables
NULL

#' GFF3 Gene information
#'
#' A dataset containing the gene information in GFF3 format.
#'     The variables are as follows:
#'
#' \itemize{
#'   \item seqid. name of the chromosome or scaffold; here, of Panicum virgatum.
#'   \item source. name of the program that generated this feature.
#'   \item type. type of feature. Term or accession from the SOFA sequence ontology.
#'   \item start. Start position of the feature, with sequence numbering starting at 1.
#'   \item end. End position of the feature, with sequence numbering starting at 1.
#'   \item score. A floating point value.
#'   \item strand. defined as + (forward) or - (reverse).
#'   \item phase. One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on.
#'   \item attributes. A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent.
#' }
#'
#' @name gff_gene
#' @docType data
#' @author Alice MacQueen \email{alice.macqueen@@utexas.edu}
#' @references \url{data_blah.com}
#' @keywords data
#' @usage data(gff_gene)
#' @format A data frame with 1236881 rows and 9 variables
NULL

#' Metadata lines for the gene GFF3
#'
#' A dataset containing the three metadata lines for the gff3 file gff_gene.
#'
#' @name gff_metadata
#' @docType data
#' @author Alice MacQueen \email{alice.macqueen@@utexas.edu}
#' @references \url{data_blah.com}
#' @keywords data
#' @usage data(gff_metadata)
#' @format A data frame with 3 rows and 1 variable.
NULL
