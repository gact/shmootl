% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{parseFilenames}
\alias{parseFilenames}
\title{Parse filenames by pattern matching.}
\usage{
parseFilenames(filenames, pattern)
}
\arguments{
\item{filenames}{Character vector of filenames to parse.}

\item{filename.pattern}{Filename pattern, which must be a valid Perl regex
with named capture groups. Neither the capture groups nor the regex itself
are required to match any given filename, but all capture groups must have
a unique name.}
}
\value{
Character matrix in which each row contains the values of capture
groups for a given filename, and each column contains values of a given
capture group across files. Unmatched capture groups are represented by
\code{NA} values.
}
\description{
Parse filenames by pattern matching.
}
\keyword{internal}
