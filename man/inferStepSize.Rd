% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{inferStepSize}
\alias{inferStepSize}
\title{Infer step size from step values.}
\usage{
inferStepSize(steps, tol = .Machine$double.eps^0.5)
}
\arguments{
\item{steps}{Numeric vector of steps (i.e. differences between consecutive 
values).}

\item{tol}{Tolerance for step equality.}
}
\value{
Inferred step size. Returns \code{NULL} if step size cannot be 
inferred.
}
\description{
Infer step size from step values.
}
\keyword{internal}
