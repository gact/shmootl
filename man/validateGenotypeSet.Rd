% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{validateGenotypeSet}
\alias{validateGenotypeSet}
\title{Validate a set of genotypes.}
\usage{
validateGenotypeSet(x, strict = FALSE)
}
\arguments{
\item{x}{Vector of genotypes.}

\item{strict}{Restrict validity to complete genotypes only.}
}
\value{
\code{TRUE} if genotype set is valid; otherwise raises first error.
}
\description{
Validate a set of genotypes.
}
\keyword{internal}
