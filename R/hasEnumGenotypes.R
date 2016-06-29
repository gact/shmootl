# Start of hasEnumGenotypes.R ##################################################

# hasEnumGenotypes -------------------------------------------------------------
#' Test if object contains enumerated genotypes.
#' 
#' @param x Object containing genotype data.
#' 
#' @return TRUE if object contains enumerated genotypes; FALSE otherwise.
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
    
    # Get cross alleles.
    alleles <- attr(x, 'alleles')
    
    # Check alleles are for enumerated genotypes.
    if ( ! all( isEnumAllele(alleles) ) ) {
        return(FALSE)
    }
    
    # Check alleles are consecutive integers starting from 1.
    allele.numbers <- sort( as.integer(alleles) )
    if ( allele.numbers[1] != 1 || any( diff(allele.numbers) != 1 ) ) {
        return(FALSE)
    }
    
    # Drop markers that don't have at least one genotype.
    x <- qtl::drop.nullmarkers(x)
    
    # Get genotype matrix.
    geno.matrix <- qtl::pull.geno(x)
    
    # Replace encoded genotypes with actual genotype values.
    for ( i in getIndices(alleles) ) {
        geno.matrix[ geno.matrix == i ] <- alleles[i]
    }
    
    # Get column indices of genotype matrix.
    cols <- getColIndices(geno.matrix)
    
    # Get row indices of first genotype values.
    rows <- sapply(cols, function(j) which( ! is.na(geno.matrix[, j]) )[1])
    
    # Get first genotype value for each marker.
    prime.geno <- sapply(getIndices(cols), function(i) geno.matrix[rows[i], cols[i]])
    
    prime.geno <- unique(prime.geno)
    
    # If there is not a unique first genotype with the
    # symbol '1', this is not enumerated genotype data.
    if ( length(prime.geno) != 1 || prime.geno != '1' ) {
        return(FALSE)
    }
 
    return(TRUE)   
}

# End of hasEnumGenotypes.R ####################################################