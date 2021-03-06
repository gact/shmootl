% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getQTLPeaks.R
\name{testQTLPeakSignificance}
\alias{testQTLPeakSignificance}
\title{Test if QTL peaks are significant.}
\usage{
testQTLPeakSignificance(x, threshold, qtl.peaks, ...)
}
\arguments{
\item{x}{A \code{scanone} object.}

\item{threshold}{A single \code{numeric} LOD significance threshold, or an
object (e.g. \code{summary.scanoneperm}) containing one such threshold
and its associated significance level.}

\item{qtl.peaks}{Locus \code{mapframe} containing QTL peaks.}

\item{...}{Further arguments. These are passed to \code{\link{getLODProfile}},
and may include a \code{lodcolumn} for a \code{scanone} object with multiple 
LOD columns, or a \code{qtl.index} for a \code{qtl} object with multiple QTLs.}
}
\value{
Logical vector indicating which QTL peaks are significant.
}
\description{
Given an object containing a LOD profile and a \code{mapframe} of QTL peaks, 
this function tests which peaks are significant with respect to the specified 
LOD threshold.
}
\keyword{internal}
