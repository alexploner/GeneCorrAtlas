% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/GeneCorrAtlas.R
\name{gcatlas}
\alias{gcatlas}
\title{Atlas of genetic correlations}
\format{A data frame with 1176 rows and 6 variables}
\source{
\url{http://www.nature.com/ng/journal/vaop/ncurrent/extref/ng.3406-S2.csv}
}
\description{
A data frame containing the genetic correlations between 50 different traits
and diseases presented in Bulik-Sullivan et al. (2015).
}
\details{
\itemize{
  \item Trait1, name of first phenotype
  \item Trait2, name of second phenotype
  \item rg, the estimated correlation coefficient
  \item se, the standard error of the estimate
  \item z, the corresponding Wald statistic
  \item p, the corresponding p-value
}
}
\references{
Bulik-Sullivan B, Finucane HK, Anttila V, Gusev A, Day FR, Loh P-R,
   et al. An atlas of genetic correlations across human diseases and traits.
   Nat Genet [Internet]. 2015 Sep 28 [cited 2015 Oct 5];advance online publication
}
