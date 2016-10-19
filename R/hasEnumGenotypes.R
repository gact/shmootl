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
    crosstype <- pull.crosstype(x)
    alleles <- pull.alleles(x)
    genotypes <- makeGenotypeSet(alleles, crosstype)
    return( all( isEnumGenotype(genotypes) ) )
}

# hasEnumGenotypes.geno --------------------------------------------------------
#' @export
#' @rdname hasEnumGenotypes
hasEnumGenotypes.geno <- function(x) {
    alleles <- pull.alleles(x)
    return( all( isEnumAllele(alleles) ) )
}

# End of hasEnumGenotypes.R ####################################################