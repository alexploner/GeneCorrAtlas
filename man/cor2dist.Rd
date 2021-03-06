% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneCorrAtlas.R
\name{cor2dist}
\alias{cor2dist}
\alias{cor2dist_abs}
\title{Convert a correlation matrix to a distance matrix}
\usage{
cor2dist(x)

cor2dist_abs(x)
}
\arguments{
\item{x}{A square matrix of correlations}
}
\value{
An object of class \code{dist}: \code{cor2dist} returns \eqn{1-r},
        \code{cor2dist2} returns \eqn{1-|r|}.
}
\description{
Two simple hepler functions that convert a matrix of correlations to
a distance/dissimilarity matrix
}
\seealso{
\code{\link{stats::dist}}, \code{\link{gcheatmap}}
}

