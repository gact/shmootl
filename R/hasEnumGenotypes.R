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
    return( all( isEnumAllele( pull.alleles(x) ) ) )
}

# End of hasEnumGenotypes.R ####################################################