% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{loadSeqTables}
\alias{loadSeqTables}
\title{Load genome sequence info tables.}
\usage{
loadSeqTables()
}
\value{
A \code{list} of \code{data.frame} objects, each element named for a
given genome and containing sequence info for that genome. Each such
\code{data.frame} object can be used to create, but is distinct from,
a \pkg{GenomeInfoDb} \code{Seqinfo} object.
}
\description{
This function loads the sequence info tables of genomes for which package
data is available.
}
\keyword{internal}
