% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CrossInfo-class.R
\docType{methods}
\name{getMarkerNames}
\alias{getMarkerNames}
\alias{getMarkerNames,CrossInfo-method}
\title{Get marker names.}
\usage{
getMarkerNames(cross.info, markers = NULL)

\S4method{getMarkerNames}{CrossInfo}(cross.info, markers = NULL)
}
\arguments{
\item{cross.info}{A \code{\link{CrossInfo}} object.}

\item{markers}{Vector indicating markers to consider. This may be a logical 
vector of markers, an integer vector of marker indices, or a character vector
that contains marker IDs or their equivalent syntactically valid name. If 
none are specified, results are returned for all markers.}
}
\value{
Character vector of marker names.
}
\description{
Get marker names.
}
\seealso{
Other CrossInfo methods: \code{\link{getAlleles}},
  \code{\link{getCrosstype}}, \code{\link{getGenotypes}},
  \code{\link{getMarkerIndices}},
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
  \code{\link{getTetradIndices}},
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
