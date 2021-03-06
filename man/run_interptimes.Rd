% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_interptimes.R
\name{run_interptimes}
\alias{run_interptimes}
\title{Interpolate time-series phenotypes in \pkg{R/qtl} CSV file.}
\usage{
run_interptimes(datafile = NA_character_, tol = 1e-05)
}
\arguments{
\item{datafile}{cross CSV file [required]}

\item{tol}{time-step equality tolerance}
}
\description{
Read data in an \pkg{R/qtl} cross CSV file. If the input \code{cross} has
time-series phenotypes, interpolate phenotype values in any gaps in the
time series, then write the resulting \code{cross} to the output file.
}
\seealso{
Other pipeline functions: \code{\link{run_annoqtl}},
  \code{\link{run_digest}}, \code{\link{run_estimap}},
  \code{\link{run_makecross}}, \code{\link{run_makegeno}},
  \code{\link{run_prep}}, \code{\link{run_pullmap}},
  \code{\link{run_pushmap}}, \code{\link{run_recode}},
  \code{\link{run_report}}, \code{\link{run_scanone}},
  \code{\link{run_scantwo}}
}
\author{
Yue Hu
}
\concept{
shmootl:utilities
}
