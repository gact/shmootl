# Start of encodeGenotypes.R ###################################################

# encodeGenotypes --------------------------------------------------------------
#' Encode genotype strings as genotype integers.
#' 
#' @param x Object containing genotype strings.
#' @param genotypes Sorted character vector of genotype symbols.
#' Input genotype strings should be elements of this vector.
#' 
#' @return Object of same class and dimensions as input object,
#' but with genotype strings encoded as genotype integers.
#' 
#' @export
#' @keywords internal
#' @rdname encodeGenotypes
encodeGenotypes <- function(x, genotypes) {
    UseMethod('encodeGenotypes', x)
}

# encodeGenotypes.character ----------------------------------------------------
#' @export
#' @method encodeGenotypes character
#' @rdname encodeGenotypes
encodeGenotypes.character <- function(x, genotypes) {
    
    validateGenotypeSet(genotypes)
    
    unknown <- x[ ! x %in% c(genotypes, NA_character_) ]
    if ( length(unknown) > 0 ) {
        stop("cannot encode unknown symbols as genotypes - '", unknown, "'")
    }
    
    x <- match(x, genotypes)
    
    return(x)
}

# encodeGenotypes.data.frame ---------------------------------------------------
#' @export
#' @method encodeGenotypes data.frame
#' @rdname encodeGenotypes
encodeGenotypes.data.frame <- function(x, genotypes) {
    for ( i in getColIndices(x) ) {
        x[, i] <- encodeGenotypes(as.character(x[, i]), genotypes)
    }
    return(x)
}

# encodeGenotypes.matrix -------------------------------------------------------
#' @export
#' @method encodeGenotypes matrix
#' @rdname encodeGenotypes
encodeGenotypes.matrix <- function(x, genotypes) {
    x <- apply(x, 2, as.character)
    x <- apply(x, 2, encodeGenotypes, genotypes)
    return(x)
}

# End of encodeGenotypes.R #####################################################