% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map.R
\name{findLocusRowIndices}
\alias{findLocusRowIndices}
\title{Find row index of closest \code{mapframe} locus.}
\usage{
findLocusRowIndices(x, loc)
}
\arguments{
\item{x}{A \code{mapframe} object.}

\item{loc}{Locus \code{mapframe} specifying map positions.}
}
\value{
Vector of row indices, each element of which indicates the row of the
target \code{mapframe} object containing the closest locus to the map position
specified in the corresponding row of the locus \code{mapframe}.
}
\description{
Find row index of closest \code{mapframe} locus.
}
\keyword{internal}
