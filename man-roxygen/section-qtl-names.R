#' @section QTL Names:
#' 
#' A QTL name can be any valid item ID (see
#' \link[=shmootl-package]{package overview}). Default QTL names are of the
#' form \code{'04@33.0'}. As with pseudomarker IDs, this indicates the
#' chromosome and genetic map position of the QTL, where \code{'04'} is a
#' chromosome ID and \code{'33.0'} is the position of the given QTL in
#' centiMorgans.
#'  
#' In a default QTL name, the minimum precision of the map position is one
#' tenth of a centiMorgan. As in \pkg{R/qtl}, if the map step size is smaller
#' than this, the default QTL name is adjusted to match the appropriate level
#' of precision.
