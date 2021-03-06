% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{makeScanoneThresholdObject}
\alias{makeScanoneThresholdObject}
\title{Make \code{scanone} threshold object.}
\usage{
makeScanoneThresholdObject(threshold, lodcolnames = "lod", alpha = NULL,
  fdr = NULL)
}
\arguments{
\item{threshold}{LOD threshold value.}

\item{lodcolnames}{Character vector of LOD column names
to be used in the \code{scanone} threshold object.}

\item{alpha}{Significance level (alpha) associated with the specified LOD
threshold. (Incompatible with \code{fdr}.)}

\item{fdr}{False-discovery rate associated with the specified LOD
threshold. (Incompatible with \code{alpha}.)}
}
\value{
If \code{alpha} is specified, a \code{summary.scanoneperm} object is
returned, with the given threshold at the specified \code{alpha} for each LOD
column. If \code{fdr} is specified, a \code{summary.scanonebins} object is
returned, with the given threshold at the specified \code{fdr} in the single
LOD column for all phenotypes. In any case, because no permutations are used
in generating the threshold, the returned object has attribute \code{'n.perm'}
containing the value \code{0}.
}
\description{
Make \code{scanone} threshold object.
}
\keyword{internal}
