# Start of decodeGenotypes.R ###################################################

# decodeGenotypes --------------------------------------------------------------
#' Decode genotype integers to genotype strings.
#' 
#' @param x Object containing genotype integers.
#' @param genotypes Sorted character vector of genotype symbols.
#' Input genotype integers should be valid indices of this vector.
#' @param missing.value If set to a missing value symbol (\code{'-'}),
#' \code{NA} values are converted to the given symbol.
#' 
#' @return Object of same class and dimensions as input object,
#' but with genotype integers decoded to genotype strings.
#' 
#' @export
#' @keywords internal
#' @rdname decodeGenotypes
decodeGenotypes <- function(x, genotypes, missing.value=NA_character_) {
    UseMethod('decodeGenotypes', x)
}

# decodeGenotypes.integer ------------------------------------------------------
#' @export
#' @method decodeGenotypes integer
#' @rdname decodeGenotypes
decodeGenotypes.integer <- function(x, genotypes, missing.value=NA_character_) {
    
    validateGenotypeSet(genotypes)
    
    exrange <- x[ ! is.na(x) & ! x %in% seq_along(genotypes) ]
    if ( length(exrange) > 0 ) {
        stop("encoded genotypes are out of range - '", toString(exrange), "'")
    }
    
    x <- genotypes[x]
    
    # If specified, convert NA values to given missing value.
    if ( ! is.na(missing.value) ) {
        stopifnot( missing.value == const$missing.value )
        x[ is.na(x) ] <- missing.value
    }
    
    return(x)
}

# decodeGenotypes.data.frame ---------------------------------------------------
#' @export
#' @method decodeGenotypes data.frame
#' @rdname decodeGenotypes
decodeGenotypes.data.frame <- function(x, genotypes,
    missing.value=NA_character_) {
    stopifnot( all( sapply(x, class) == 'integer' ) )
    for ( i in getColIndices(x) ) {
        x[, i] <- decodeGenotypes(x[, i], genotypes,
            missing.value=missing.value)
    }
    return(x)
}

# decodeGenotypes.matrix -------------------------------------------------------
#' @export
#' @method decodeGenotypes matrix
#' @rdname decodeGenotypes
decodeGenotypes.matrix <- function(x, genotypes, missing.value=NA_character_) {
    stopifnot( typeof(x) == 'integer' )
    x <- apply(x, 2, decodeGenotypes, genotypes, missing.value=missing.value)
    return(x)
}

# End of decodeGenotypes.R #####################################################