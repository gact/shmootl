#' @section Mapframe:
#' 
#' A \code{mapframe} object is a class of \code{data.frame} that contains mapped
#' information, or a \code{data.frame} with equivalent properties. The leftmost
#' two columns of every \code{mapframe} object must be as follows:
#' 
#' \itemize{
#' \item{\code{chr}:}   {Chromosome/sequence labels in
#' the order defined for the current reference genome.}
#' \item{\code{pos}:}   {Map positions in increasing order within each
#' chromosome/sequence. The map unit used in a \code{mapframe} position
#' column must be indicated by setting its \code{'map.unit'} attribute.}
#' }
#' 
#' Any number of additional data columns can be included. If present, the row
#' names indicate the locus IDs, which should be unique within the given
#' \code{mapframe}.
