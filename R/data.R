#' Phenotypes available for the Panicum virgatum diversity panel.
#'
#' A dataset containing phenotypes for the Panicum virgatum diversity panel.
#'     The variables are as follows:
#'
#' \itemize{
#'   \item PLANT_ID. ID of the Panicum virgatum genotype sequenced in the
#'      diversity panel.
#'   \item GWAS_CT. The number of times that genotype was clonally propagated
#'       and planted at these common gardens.
#'   \item MAT. Mean annual temperature at the genotype's location of origin.
#'   \item bio17. Precipitation in the driest quarter at the genotype's location
#'       of origin.
#'   \item bio4. Temperature seasonality at the genotype's location of origin.
#'   \item bio16. Precipitation in the wettest quarter at the genotype's
#'      location of origin.
#'   \item AHM. Annual heat moisture index, calculated as (MAT+10)/(MAP/1000),
#'      at the genotype's location of origin.
#'   \item bio2. Mean diurnal temperature range at the genotype's location of
#'      origin.
#'   \item bio5. Max temperature of the warmest month at the genotype's location
#'      of origin.
#'   \item FRAC_SRV_THREE. Fraction of individuals surviving winter 2018/2019
#'      at three northern gardens (BRKG, LINC, CLMB).
#'   \item CLMB_BIOMASS. Biomass at the end of the 2019 season in Columbia,
#'      Missouri.
#'   \item KBSM_BIOMASS. Biomass at the end of the 2019 season in Hickory
#'      Corners, Michigan (Kellog Biological Station).
#'   \item PKLE_BIOMASS. Biomass at the end of the 2019 season in Austin, Texas
#'       (Pickle Research Campus).
#' }
#'
#' @name pvdiv_phenotypes
#' @docType data
#' @author Alice MacQueen \email{alice.macqueen@@utexas.edu}
#' @references \url{switchgrassgenomepaper.link}
#' @keywords data
#' @usage data(pvdiv_phenotypes)
#' @format A data frame with 1114 rows and 13 variables
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
#' @keywords data
#' @usage data(gff_metadata)
#' @format A data frame with 3 rows and 1 variable.
NULL

#' Annotation information from Pvirgatum_516_v5.1.annotation_info
#'
#' A dataset containing the v5.1 annotation information for Panicum virgatum.
#'     The variables are as follows:
#'
#' \itemize{
#'   \item pacID ID of the Panicum virgatum gene from the annotation file.
#'   \item locusName The gene name from Panicum virgatum. Only distinct locus
#'       names are included in this dataset.
#'   \item transcriptName The transcript name from P. virgatum. Only the first
#'       transcript variant is included in this dataset.
#'   \item peptideName The protein name for the variant. Only the first protein
#'       variant is included in this dataset.
#'   \item Pfam
#'   \item Panther Panther categories associated with this locus.
#'   \item KOG assignments from the euKaryotic Orthologous Groups (KOG) database
#'   \item ec EC numbers in Enzyme Nomenclature.
#'   \item KO KEGG Orthology (KO) numbers associated with this locusName.
#'   \item GO Gene ontology entries associated with this locusName.
#'   \item Best-hit-arabi-name Best Arabidopsis thaliana homolog locus name.
#'   \item arabi-symbol Best A. thaliana homolog gene name.
#'   \item arabi-defline Best A. thaliana homolog gene annotation.
#'   \item Best-hit-rice-name Best rice homolog locus name.
#'   \item rice-symbol Best rice homolog gene name.
#'   \item rice-defline Best rice homolog gene annotation.
#' }
#'
#' @name anno_info
#' @docType data
#' @author Alice MacQueen \email{alice.macqueen@@utexas.edu}
#' @references \url{switchgrassgenomepaper.link}
#' @keywords data
#' @usage data(anno_info)
#' @format A data frame with 80278 rows and 16 variables.
NULL

#' Oeco ggplot theme
#'
#' A list containing the oeco ggplot2 theme for use in cdbn_standard_gwas and
#'     mash_plot_ ... functions.
#'
#' @name theme_oeco
#' @docType data
#' @author Alice MacQueen \email{alice.macqueen@@utexas.edu}
#' @references \url{switchgrassgenomepaper.link}
#' @keywords data
#' @usage data(theme_oeco)
#' @format A list containing ggplot2 specifications
NULL

#' Metadata table for the Panicum virgatum diversity panel
#'
#' A dataset containing the collected metadata information for individuals in
#'     the Panicum virgatum diversity panel.
#'     The variables are as follows:
#'
#' \itemize{
#'   \item PLANT_ID Individual plant name. Unique to each plant. The column we
#'       should be using for analyses.
#'   \item SUBPOP Using clustering from an additive kinship matrix, which
#'       genetic subpopulation does the PLANT_ID belong to? 4X and 8X
#'       individuals have been assigned a ploidy by flow cytometry, but not a
#'       subpopulation by sequencing, as of the end of 2020.
#'   \item LIBRARY Four letter library code used for sequencing. IF NA, the
#'       plant has not been sequenced as of the end of 2020.
#'   \item GENETIC_SUBPOP_KINSHIP Using clustering from an additive kinship
#'       matrix, which genetic subpopulation does the PLANT_ID belong to? NA
#'       individuals were not included in the kinship matrix.
#'   \item KINSHIP_ORDER Order of individuals in an additive kinship matrix.
#'       Roughly speaking, closer numbers represent more similar individuals.
#'   \item GENETIC_SUBPOP_50PER From hierarchical k=3 analysis done early 2020.
#'       Which genetic subpopulation does the PLANT_ID belong to, at 50%
#'       assignment or more. Admixed individuals are unassigned PLANT_IDs
#'       included in the analysis; NA individuals were not included.
#'   \item GENETIC_SUBPOP_95PER From hierarchical k=3 analysis done early 2020.
#'       Which genetic subpopulation does the PLANT_ID belong to, at 95%
#'       assignment or more. Admixed individuals are unassigned PLANT_IDs
#'       included in the analysis; NA individuals were not included.
#'   \item GULF_Q From hierarchical k=3 analysis done early 2020. Structure Q
#'       value for assignment to the Gulf genetic subpopulation. NA individuals
#'       were not included in the Structure analysis.
#'   \item ATLANTIC_Q From hierarchical k=3 analysis done winter 2020. Structure
#'       Q value for assignment to the Atlantic subpopulation. NA individuals
#'       were not included in the Structure analysis.
#'   \item MIDWEST_DA From hierarchical k=3 analysis done winter 2020. DAPC
#'       value for assignment to the Midwest subpopulation. NA individuals were
#'       not included in the DAPC analysis.
#'   \item ECOTYPE From John's phenotypic analysis of Ecotypes. Which phenotypic
#'       group does the PLANT_ID belong to, at 50% assignment or more. NA
#'       individuals did not have phenotypic data for this classification.
#'   \item ECOTYPE_NNET From John's neural network extending ecotype
#'       classification to individuals without phenotypic data.
#'   \item ECOTYPE_PC1 From John's phenotypic analysis of ecotypes done summer
#'       2019. The first principle component separating individuals into
#'       ecotypes.
#'   \item UPLAND_ASSIGNMENT From John's phenotypic analysis of Ecotypes done
#'       summer 2019. Assignment to the Upland ecotype.
#'   \item LOWLAND_ASSIGNMENT From John's phenotypic analysis of Ecotypes done
#'       summer 2019. Assignment to the Texas ecotype.
#'   \item COASTAL_ASSIGNMENT: From John's phenotypic analysis of Ecotypes done
#'       summer 2019. Assignment to the Coastal ecotype.
#'   \item LIB_GROWN For tetraploid plants with sequencing libraries, was the
#'       PLANT_ID grown at one or more of the 10 common garden sites in 2018?
#'       "Y" = yes, "N" = no. Does not include plants where LIB_PHENO is "N".
#'   \item LIB_BIOCLIM Do we have bioclim variables for this PLANT_ID?
#'       "Y" = yes, "N" = no. Does not include plants where LIB_CLIMATE is "N".
#'   \item LIB_CLIMATE John & Jason's decision on whether this PLANT_ID should
#'       be used for climate analyses done summer 2019. "Y" = yes, "N" = no.
#'   \item LIB_PHENO John & Jason's decision on whether this PLANT_ID should be
#'       used for phenotypic analyses done summer 2019. "Y" = yes, "N" = no.
#'   \item LIB_ACC John & Jason's decision on whether this PLANT_ID should be
#'       used for analyses of a single plant per location done summer 2019 .
#'       "Y" = yes, "N" = no.
#'   \item LIB_NOTES Notes on John & Jason's decisions for the above 3 columns.
#'       This analysis was done summer 2019.
#'   \item PLOIDY Ploidy as determined by both sequencing & flow cytometry
#'       analyses done in 2019. Typically 4X (tetraploid) or 8X (octoploid).
#'   \item COLLECTION_TYPE Original plant collection metadata. What kind of
#'       collection was this? Breeding selection, Natural collection, Cultivar.
#'   \item COLLECTION_METHOD How was the plant collected?
#'   \item COLLECTOR Who collected this plant?
#'   \item STATE State of the United States that the plant was collected from.
#'   \item COUNTY County that the plant was collected from.
#'   \item LATITUDE Latitude of origin of the collected plant.
#'   \item LONGITUDE Longitude of origin of the collected plant.
#'   \item ELEVATION Elevation of origin of the collected plant.
#'   \item NOTE_LATLONG Notes on exactness of latitude & longitude of collection
#'   \item LOCALITY Additional details on plant collection location.
#'   \item TAXON Genus & species of the collected plant.
#'   \item HABITAT Details on habitat of plant collection location, if known.
#'   \item COLL_DATE Information on plant collection date, if known.
#'   \item PLANTED_2018 All plants clonally replicated and transplanted into
#'       one or more common garden in the summer of 2018.
#' }
#'
#' @name pvdiv_metadata
#' @docType data
#' @author Alice MacQueen \email{alice.macqueen@@utexas.edu}
#' @references \url{switchgrassgenomepaper.link}
#' @keywords data
#' @usage data(pvdiv_metadata)
#' @format A data frame with 1347 rows and 37 variables.
NULL
