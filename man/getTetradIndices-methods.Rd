% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CrossInfo-class.R
\docType{methods}
\name{getTetradIndices}
\alias{getTetradIndices}
\alias{getTetradIndices,CrossInfo-method}
\title{Get tetrad indices for the given samples.}
\usage{
getTetradIndices(cross.info, samples = NULL)

\S4method{getTetradIndices}{CrossInfo}(cross.info, samples = NULL)
}
\arguments{
\item{cross.info}{A \code{\link{CrossInfo}} object.}

\item{samples}{Vector indicating samples to consider. This may be a 
logical vector of samples, an integer vector of sample indices, or 
a character vector that contains sample IDs or their equivalent 
syntactically valid name.}
}
\value{
Integer vector of tetrad indices for the given samples. If no samples 
are specified, all tetrad indices are returned. Returns \code{NULL} if tetrad 
indices are not present.
}
\description{
Get tetrad indices for the given samples.
}
\seealso{
Other CrossInfo methods: \code{\link{getAlleles}},
  \code{\link{getCrosstype}}, \code{\link{getGenotypes}},
  \code{\link{getMarkerIndices}},
  \code{\link{getMarkerNames}},
  \code{\link{getMarkerSeqs}}, \code{\link{getMarkers}},
  \code{\link{getNumMarkers}},
  \code{\link{getNumPhenotypes}},
  \code{\link{getNumSamples}}, \code{\link{getNumSeqs}},
  \code{\link{getPhenotypeIndices}},
  \code{\link{getPhenotypeNames}},
  \code{\link{getPhenotypes}},
  \code{\link{getSampleIndices}},
  \code{\link{getSampleNames}}, \code{\link{getSamples}},
  \code{\link{getSeqIndices}}, \code{\link{getSeqMarkers}},
  \code{\link{getSeqNames}}, \code{\link{getSequences}},
  \code{\link{getStrainIndices}},
  \code{\link{hasMarkerSeqs}}, \code{\link{hasSampleIDs}},
  \code{\link{hasStrainIndices}},
  \code{\link{hasTetradIndices}}, \code{\link{setAlleles}},
  \code{\link{setCrosstype}}, \code{\link{setGenotypes}},
  \code{\link{setMarkerSeqs}}, \code{\link{setMarkers}},
  \code{\link{setPhenotypes}}, \code{\link{setSamples}},
  \code{\link{setSequences}},
  \code{\link{setStrainIndices}},
  \code{\link{setTetradIndices}}
}
\keyword{internal}
