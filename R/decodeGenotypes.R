# Start of decodeGenotypes.R ###################################################

# decodeGenotypes --------------------------------------------------------------
#' Decode genotype integers to genotype strings.
#' 
#' @param x Object containing genotype integers.
#' @param genotypes Sorted character vector of genotype symbols.
#' Input genotype integers should be valid indices of this vector.
#' 
#' @return Object of same class and dimensions as input object,
#' but with genotype integers decoded to genotype strings.
#' 
#' @export
#' @keywords internal
#' @rdname decodeGenotypes
decodeGenotypes <- function(x, genotypes) {
    UseMethod('decodeGenotypes', x)
}

# decodeGenotypes.integer ------------------------------------------------------
#' @export
#' @method decodeGenotypes integer
#' @rdname decodeGenotypes
decodeGenotypes.integer <- function(x, genotypes) {
    
    validateGenotypeSet(genotypes)
    
    exrange <- x[ ! is.na(x) & ! x %in% seq_along(genotypes) ]
    if ( length(exrange) > 0 ) {
        stop("encoded genotypes are out of range - '", toString(exrange), "'")
    }
    
    x <- genotypes[x]
    
    return(x)
}

# decodeGenotypes.data.frame ---------------------------------------------------
#' @export
#' @method decodeGenotypes data.frame
#' @rdname decodeGenotypes
decodeGenotypes.data.frame <- function(x, genotypes) {
    stopifnot( all( sapply(x, class) == 'integer' ) )
    for ( i in getColIndices(x) ) {
        x[, i] <- decodeGenotypes(x[, i], genotypes)
    }
    return(x)
}

# decodeGenotypes.matrix -------------------------------------------------------
#' @export
#' @method decodeGenotypes matrix
#' @rdname decodeGenotypes
decodeGenotypes.matrix <- function(x, genotypes) {
    stopifnot( typeof(x) == 'integer' )
    x <- apply(x, 2, decodeGenotypes, genotypes)
    return(x)
}

# End of decodeGenotypes.R #####################################################