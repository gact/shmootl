% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdf5.R
\name{readResultsOverviewHDF5}
\alias{readResultsOverviewHDF5}
\title{Read QTL analysis results overview from an input HDF5 file.}
\usage{
readResultsOverviewHDF5(infile)
}
\arguments{
\item{infile}{An input HDF5 file.}
}
\value{
A \code{data.frame} containing a QTL analysis results
overview. Returns \code{NULL} if no results overview is found.
}
\description{
Read QTL analysis results overview from an input HDF5 file.
}
\seealso{
Other HDF5 functions: \code{\link{readCrossHDF5}},
  \code{\link{readMapHDF5}}, \code{\link{readResultHDF5}},
  \code{\link{writeCrossHDF5}}, \code{\link{writeMapHDF5}},
  \code{\link{writeResultHDF5}},
  \code{\link{writeResultsOverviewHDF5}}
}
