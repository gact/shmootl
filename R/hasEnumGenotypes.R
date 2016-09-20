# Start of hasEnumGenotypes.R ##################################################

# hasEnumGenotypes -------------------------------------------------------------
#' Test if object contains enumerated genotypes.
#' 
#' @param x Object containing genotype data.
#' 
#' @return \code{TRUE} if object contains enumerated genotypes;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname hasEnumGenotypes
hasEnumGenotypes <- function(x) {
    UseMethod('hasEnumGenotypes', x)
}

# hasEnumGenotypes.cross -------------------------------------------------------
#' @export
#' @rdname hasEnumGenotypes
hasEnumGenotypes.cross <- function(x) {
    stopifnot( 'cross' %in% class(x) )
    stopifnot( attr(x, 'crosstype') == 'bc' ) # TODO: support other cross types
    return( all( isEnumAllele( attr(x, 'alleles') ) ) )
}

# End of hasEnumGenotypes.R ####################################################