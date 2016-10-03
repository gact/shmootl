# Start of hasFounderGenotypes.R ###############################################

# hasFounderGenotypes ----------------------------------------------------------
#' Test if object contains founder genotypes.
#' 
#' @param x Object containing genotype data.
#' 
#' @return \code{TRUE} if object contains founder genotypes;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname hasFounderGenotypes
hasFounderGenotypes <- function(x) {
    UseMethod('hasFounderGenotypes', x)
}

# hasFounderGenotypes.cross ----------------------------------------------------
#' @export
#' @rdname hasFounderGenotypes
hasFounderGenotypes.cross <- function(x) {
    stopifnot( 'cross' %in% class(x) )
    stopifnot( attr(x, 'crosstype') == 'bc' ) # TODO: support other cross types
    return( all( isFounderAllele( attr(x, 'alleles') ) ) )
}

# hasFounderGenotypes.geno -----------------------------------------------------
#' @export
#' @rdname hasFounderGenotypes
hasFounderGenotypes.geno <- function(x) {
    stopifnot( 'geno' %in% class(x) )
    return( all( isFounderAllele( attr(x, 'alleles') ) ) )
}

# End of hasFounderGenotypes.R #################################################