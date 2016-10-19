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
    crosstype <- pull.crosstype(x)
    alleles <- pull.alleles(x)
    genotypes <- makeGenotypeSet(alleles, crosstype)
    return( all( isFounderGenotype(genotypes) ) )
}

# hasFounderGenotypes.geno -----------------------------------------------------
#' @export
#' @rdname hasFounderGenotypes
hasFounderGenotypes.geno <- function(x) {
    alleles <- pull.alleles(x)
    return( all( isFounderAllele(alleles) ) )
}

# End of hasFounderGenotypes.R #################################################