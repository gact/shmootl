% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gff.R
\name{readFeaturesGFF}
\alias{readFeaturesGFF}
\title{Read annotation from GFF file.}
\usage{
readFeaturesGFF(annofile)
}
\arguments{
\item{annofile}{Input annotation GFF file.}
}
\value{
A \code{data.frame} of primary genome features, such that each
row contains at least the reference sequence (\code{'chr'}), start and
end positions (\code{'start'} and \code{'end'} respectively), and \code{'ID'}
of a feature. Additional feature annotation is included, if available.
}
\description{
Read annotation from GFF file.
}
