% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wasa-data.R
\docType{data}
\name{wasa}
\alias{wasa}
\title{WAxSA Cross Data}
\format{An \pkg{R/qtl} \code{cross} object, with cross information contained
in attribute \code{'info'}. For more information on the \code{cross} object,
see the \pkg{R/qtl} function \code{read.cross}. For more information on the
\code{'info'} attribute, see the \code{\linkS4class{CrossInfo}} class
documentation.}
\source{
Supporting information from \href{http://dx.doi.org/10.1111/j.1365-294X.2011.05005.x}{Cubillos \emph{et al.} (2011)}
}
\usage{
data(wasa)
}
\description{
Data from a QTL analysis of haploid tetrad spores from an F1
cross between two strains of \emph{S. cerevisiae}: the West African (WA)
strain \code{'DBVPG6044'} and the Sake (SA) strain \code{'Y12'} (Cubillos
\emph{et al.} 2011).
}
\details{
This \code{cross} object contains 96 segregants from the WAxSA F1
cross, with 39 phenotypes, and 183 markers across the 16 nuclear chromosomes
of the S288C reference genome.

Segregant strains were growth phenotyped in 13 different environmental
conditions, and 3 growth variables were extracted from the growth curve
of each strain: the population doubling time during the exponential phase
(rate), the total growth on reaching the stationary phase (efficiency),
and the time lag required for the strain culture to adapt to the given
environmental condition and reach the exponential growth phase (adaptation)
 (Cubillos \emph{et al.} 2011; Warringer \emph{et al.} 2008).
}
\references{
Cubillos FA, Billi E, Zorgo E, Parts L, Fargier P, Omholt S,
Blomberg A, Warringer J, Louis EJ, Liti G (2011) Assessing the complex
architecture of polygenic traits in diverged yeast populations.
\emph{Molecular Ecology} \bold{20}(7):1401-13.
(\href{http://www.ncbi.nlm.nih.gov/pubmed/21261765}{PubMed})

Warringer J, Anevski D, Liu B, Blomberg A (2008) Chemogenetic
fingerprinting by analysis of cellular growth dynamics.
\emph{BMC Chemical Biology} \bold{8}:3.
(\href{http://www.ncbi.nlm.nih.gov/pubmed/18721464}{PubMed})
}
\keyword{data}
