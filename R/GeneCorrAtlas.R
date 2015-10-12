#' Display the Atlas of Genetic Correlations 
#'
#' Bulik-Sullivan et al. recently (2105) published an atlas of genetic
#' correlations for 50 different traits. This package wraps the more than 1000
#' correlation coefficients involved and offers some convenient functions to
#' extract and display the data
#'
#'
#' @source \url{http://www.nature.com/ng/journal/vaop/ncurrent/extref/ng.3406-S2.csv}
#' @references Bulik-Sullivan B, Finucane HK, Anttila V, Gusev A, Day FR, Loh P-R,
#'    et al. An atlas of genetic correlations across human diseases and traits.
#'    Nat Genet [Internet]. 2015 Sep 28 [cited 2015 Oct 5];advance online publication
#'
#' @docType package
#' @name GeneCorrAtlas
NULL


#' Atlas of genetic correlations
#'
#' A data frame containing the genetic correlations between 50 different traits
#' and diseases presented in Bulik-Sullivan et al. (2015).
#'
#' \itemize{
#'   \item Trait1, name of first phenotype
#'   \item Trait2, name of second phenotype
#'   \item rg, the estimated correlation coefficient
#'   \item se, the standard error of the estimate
#'   \item z, the corresponding Wald statistic
#'   \item p, the corresponding p-value
#' }
#'
#' @format A data frame with 1176 rows and 6 variables. The correlations are
#' given in the order as the appear in the lower triangle of the full correlation
#' matrix, by column. 
#' @source \url{http://www.nature.com/ng/journal/vaop/ncurrent/extref/ng.3406-S2.csv}
#' @references Bulik-Sullivan B, Finucane HK, Anttila V, Gusev A, Day FR, Loh P-R,
#'    et al. An atlas of genetic correlations across human diseases and traits.
#'    Nat Genet [Internet]. 2015 Sep 28 [cited 2015 Oct 5];advance online publication
#' @name gcatlas
NULL


#' Read Atlas of Genetic Correlations
#'
#' Read the supplementary file from disk or from the web. This is the function
#' originally used to generate the pre-packaged data object \code{\link{gcatlas}}.
#' Only use this function if you want to roll your own.
#'
#' @param x name of the supplementary file (by default, a link to the original
#'    publication)
#'
#' @return The data frame described at \code{\link{gcatlas}}
#'
#' @source \url{http://www.nature.com/ng/journal/vaop/ncurrent/extref/ng.3406-S2.csv}
#' @references Bulik-Sullivan B, Finucane HK, Anttila V, Gusev A, Day FR, Loh P-R,
#'    et al. An atlas of genetic correlations across human diseases and traits.
#'    Nat Genet [Internet]. 2015 Sep 28 [cited 2015 Oct 5];advance online publication
read_gcatlas = function(x)
{
	if (missing(x)) x = "http://www.nature.com/ng/journal/vaop/ncurrent/extref/ng.3406-S2.csv"
	read.csv(x)
}	


#' Select traits from the full atlas
#'
#' Subsets the a data frame with the same format as \code{\link{gcatlas}} to
#' include only traits from a specified list of trait names. These names can be
#' abbreviated.
#'
#' @param x a vector of trait names (possibly abbreviated)
#' @param data the name of the object to subset; by default the full data set
#' @param drop logical expression indicating whether to return the specified traits
#' (\code{FALSE}, default) or to exlclude them from the data (\code{TRUE}).
#'
#' @seealso \code{\link{gcatlas}}
sel_trt = function(x, data=gcatlas, drop=FALSE)
{


}

