% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map.R
\name{as.map}
\alias{as.map}
\alias{as.map.data.frame}
\alias{as.map.list}
\alias{as.map.map}
\alias{as.map.mapframe}
\alias{as.map.scanone}
\title{Convert object to \code{map}.}
\usage{
as.map(from, map.unit = NULL)

\method{as.map}{data.frame}(from, map.unit = NULL)

\method{as.map}{list}(from, map.unit = NULL)

\method{as.map}{map}(from, map.unit = NULL)

\method{as.map}{mapframe}(from, map.unit = NULL)

\method{as.map}{scanone}(from, map.unit = NULL)
}
\arguments{
\item{from}{Object containing map data.}

\item{map.unit}{Explicitly sets \code{'map.unit'} attribute. This must be
specified if an object does not have an existing \code{'map.unit'} attribute.}
}
\value{
An \pkg{R/qtl} \code{map} corresponding to the input object.
}
\description{
Convert object to \code{map}.
}
\section{Map}{


The \pkg{R/qtl} \code{map} object - described in detail in Section 2.6 of
Broman and Sen (2009) - comprises a nonempty list, each element of which
contains two or more marker positions for a given chromosome. Each list
element is named for its corresponding chromosome, while the names of each
chromosome's marker position vector comprise the marker locus IDs.

Within \pkg{shmootl}, a \code{map} must represent either a genetic or a
physical map. These are distinguished by the \code{'map.unit'} attribute
of an object, which can contain either genetic or physical map units.
}

\seealso{
Other map utility functions: \code{\link{as.data.frame.mapframe}},
  \code{\link{as.data.frame.map}},
  \code{\link{as.mapframe}}, \code{\link{extractMarkers}},
  \code{\link{extractPseudomarkers}},
  \code{\link{findFlanking}}, \code{\link{findLoci}},
  \code{\link{findMarkers}}, \code{\link{getMapSteps}},
  \code{\link{getMapUnit}}, \code{\link{gmapframe}},
  \code{\link{hasMapCSV}}, \code{\link{inMapOrder}},
  \code{\link{inferMapStep}}, \code{\link{mapframe}},
  \code{\link{mapsEqual}}, \code{\link{orderMap}},
  \code{\link{readMapCSV}}, \code{\link{readMapHDF5}},
  \code{\link{readMapframeCSV}}, \code{\link{setMapUnit}},
  \code{\link{writeMapCSV}}, \code{\link{writeMapHDF5}},
  \code{\link{writeMapframeCSV}}
}
