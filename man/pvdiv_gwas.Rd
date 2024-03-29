% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvdiv_gwas.R
\name{pvdiv_gwas}
\alias{pvdiv_gwas}
\title{Wrapper for bigsnpr for GWAS on Panicum virgatum.}
\usage{
pvdiv_gwas(
  df,
  type = c("linear", "logistic"),
  snp,
  covar = NA,
  ncores = 1,
  npcs = 10,
  saveoutput = FALSE
)
}
\arguments{
\item{df}{Dataframe of phenotypes where the first column is PLANT_ID.}

\item{type}{Character string. Type of univarate regression to run for GWAS.
Options are "linear" or "logistic".}

\item{snp}{Genomic information to include for Panicum virgatum. SNP data
is available at doi:10.18738/T8/ET9UAU#'}

\item{covar}{Optional covariance matrix to include in the regression. You
can generate these using \code{bigsnpr::snp_autoSVD()}.}

\item{ncores}{Number of cores to use. Default is one.}

\item{npcs}{Number of principle components to use. Default is 10.}

\item{saveoutput}{Logical. Should output be saved as a rds to the
working directory?}
}
\value{
The gwas results for the last phenotype in the dataframe. That
phenotype, as well as the remaining phenotypes, are saved as RDS objects
in the working directory.
}
\description{
Given a dataframe of phenotypes associated with PLANT_IDs, this
function is a wrapper around bigsnpr functions to conduct linear or
logistic regression on Panicum virgatum. The main advantages of this
function over just using the bigsnpr functions is that it automatically
removes individual genotypes with missing phenotypic data, that it
converts switchgrass chromosome names to the format bigsnpr requires,
and that it can run GWAS on multiple phenotypes sequentially.
}
