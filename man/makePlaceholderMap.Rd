% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map.R
\name{makePlaceholderMap}
\alias{makePlaceholderMap}
\title{Make a placeholder map.}
\usage{
makePlaceholderMap(locus.seqs, locus.ids = NULL, map.unit = "cM",
  step = 5)
}
\arguments{
\item{locus.seqs}{Sequences corresponding to individual loci.}

\item{locus.ids}{Individual locus IDs.}

\item{map.unit}{Map unit.}

\item{step}{Map step size.}
}
\description{
Make a placeholder map from the given combination of locus sequences, map
unit, and map step size. If locus IDs are specified, the locus IDs of the
placeholder map will be set from these. Otherwise, locus IDs are set to
pseudomarker IDs or default marker IDs, according to whether map units
are for a genetic or physical map, respectively.
}
\details{
Locus sequences and locus IDs are assumed to be ordered. As in \pkg{R/qtl},
any sequence with too few loci is extended, to ensure that all sequences
have enough loci.
}
\section{Map}{


The \pkg{R/qtl} \code{map} object - described in detail in Section 2.6 of
Broman and Sen (2009) - comprises a nonempty list, each element of which
contains two or more marker positions for a given chromosome. Each list
element is named for its corresponding chromosome, while the names of each
chromosome's marker position vector comprise the marker locus IDs.

Within \pkg{shmootl}, a \code{map} must represent either a genetic or a
physical map. These are distinguished by the \code{'map.unit'} attribute
of an object, which can contain either genetic or physical map units.
}

\keyword{internal}
