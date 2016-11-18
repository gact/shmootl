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
    return( all( isFounderAllele( pull.alleles(x) ) ) )
}

# End of hasFounderGenotypes.R #################################################