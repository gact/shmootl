% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/H5Stack-class.R
\docType{methods}
\name{openStack}
\alias{openStack}
\alias{openStack,H5Stack-method}
\title{Open stack of HDF5 objects.}
\usage{
openStack(h5stack, file, h5name = NULL)

\S4method{openStack}{H5Stack}(h5stack, file, h5name = NULL)
}
\arguments{
\item{h5stack}{A \code{H5Stack} object.}

\item{file}{A HDF5 file path.}

\item{h5name}{Name of HDF5 object to access.}
}
\value{
Opened \code{H5Stack} object.
}
\description{
Open stack of HDF5 objects.
}
\keyword{internal}
