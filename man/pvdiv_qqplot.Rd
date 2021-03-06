% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvdiv_gwas.R
\name{pvdiv_qqplot}
\alias{pvdiv_qqplot}
\title{Create a quantile-quantile plot with ggplot2.}
\usage{
pvdiv_qqplot(ps, ci = 0.95, lambdaGC = FALSE, tol = 1e-08)
}
\arguments{
\item{ps}{Numeric vector of p-values.}

\item{ci}{Numeric. Size of the confidence interval, 0.95 by default.}

\item{lambdaGC}{Logical. Add the Genomic Control coefficient as subtitle to
the plot?}

\item{tol}{Numeric. Tolerance for optional Genomic Control coefficient.}
}
\value{
A ggplot2 plot.
}
\description{
Assumptions for this quantile quantile plot:
Expected P values are uniformly distributed.
Confidence intervals assume independence between tests.
We expect deviations past the confidence intervals if the tests are
not independent.
For example, in a genome-wide association study, the genotype at any
position is correlated to nearby positions. Tests of nearby genotypes
will result in similar test statistics.
}
