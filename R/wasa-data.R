# Start of wasa-data.R #########################################################

# WAxSA Cross Data -------------------------------------------------------------
#' WAxSA Cross Data
#' 
#' @description  Data from a QTL analysis of haploid tetrad spores from an F1
#' cross between two strains of \emph{S. cerevisiae}: the West African (WA)
#' strain \code{'DBVPG6044'} and the Sake (SA) strain \code{'YPS128'}
#' (Cubillos \emph{et al.} 2011).
#' 
#' @details This \code{cross} object contains 96 segregants from the WAxSA F1
#' cross, with 39 phenotypes, and 183 markers across 16 nuclear chromosomes.
#' 
#' Segregant strains were growth phenotyped in 13 different environmental
#' conditions, and 3 growth variables were extracted from the growth curve
#' of each strain: the population doubling time during the exponential phase
#' (rate), the total growth on reaching the stationary phase (efficiency),
#' and the time lag required for the strain culture to adapt to the given
#' environmental condition and reach the exponential growth phase (adaptation)
#'  (Cubillos \emph{et al.} 2011; Warringer \emph{et al.} 2008).
#' 
#' @usage data(wasa)
#' 
#' @format An \pkg{R/qtl} \code{cross} object, with cross information contained
#' in attribute \code{'info'}. For more information on the \code{cross} object,
#' see the \pkg{R/qtl} function \code{read.cross}. For more information on the
#' \code{'info'} attribute, see the \code{\linkS4class{CrossInfo}} class
#' documentation.
#' 
#' @source Supporting information from \href{http://dx.doi.org/10.1111/j.1365-294X.2011.05005.x}{Cubillos \emph{et al.} (2011)}
#' 
#' @template ref-cubillos-2011
#' @template ref-warringer-2008
#' 
#' @docType data
#' @family datasets
#' @keywords data
'wasa'

# End of wasa-data.R ###########################################################