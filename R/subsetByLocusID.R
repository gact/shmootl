# Start of subsetByLocusID.R ###################################################

# subsetByLocusID (S3) ---------------------------------------------------------
#' Subset object by locus ID.
#' 
#' @param x Object with locus IDs.
#' @param predicate Function applied to locus IDs that determines which loci are
#' retained in the output object.
#'     
#' @return Input object containing only loci whose locus IDs are evaluated as 
#' true by the predicate function.
#' 
#' @keywords internal
#' @rdname subsetByLocusID
subsetByLocusID <- function(x, predicate) {
    stopifnot( is.function(predicate) )
    UseMethod('subsetByLocusID', x)
}

# subsetByLocusID.data.frame ---------------------------------------------------
#' @rdname subsetByLocusID
subsetByLocusID.data.frame <- function(x, predicate) { 
    
    if ( nrow(x) > 0 ) {
        
        if ( ! hasRownames(x) ) {
            stop("cannot subset by locus ID - no locus IDs found")
        }
        
        # Extract any loci satisfying the predicate function.
        x <- x[ sapply(rownames(x), predicate), ]
    }
    
    return(x)
}

# subsetByLocusID.map ----------------------------------------------------------
#' @rdname subsetByLocusID
subsetByLocusID.map <- function(x, predicate) {
    
    # Extract filtered map data for each sequence.
    for ( i in 1:length(x) ) {
        x[[i]] <- x[[i]][ sapply(names(x[[i]]), predicate) ]
    }
    
    # Keep only sequences with enough markers.
    x <- subsetMap(x, lengths(x) >= const$min.lps)
    
    if ( length(x) < const$min.spm ) {
        stop("subsetted map has too few sequences (min=", const$min.spm, ")")
    }
    
    return(x)
}

# subsetByLocusID.DNAStringSet -------------------------------------------------
#' @rdname subsetByLocusID
subsetByLocusID.DNAStringSet <- function(x, predicate) {
    
    stopifnot( 'loci' %in% names(x@metadata) )
    stopifnot( 'mapframe' %in% class(x@metadata[['loci']]) )
    
    if ( ! hasRownames(x@metadata[['loci']]) ) {
        stop("cannot subset by locus ID - no locus IDs found")
    }
    
    indices <- which( sapply(rownames(x@metadata[['loci']]), predicate) )
    
    x@metadata[['loci']] <- x@metadata[['loci']][indices, ]
        
    for ( i in 1:length(x) ) {
        x[[i]] <- x[[i]][indices]
    }
    
    return(x)
}

# subsetByLocusID.QualityScaledDNAStringSet ------------------------------------
#' @rdname subsetByLocusID
subsetByLocusID.QualityScaledDNAStringSet <- function(x, predicate) {

    stopifnot( 'loci' %in% names(x@metadata) )
    stopifnot( 'mapframe' %in% class(x@metadata[['loci']]) )
    
    if ( ! hasRownames(x@metadata[['loci']]) ) {
        stop("cannot subset by locus ID - no locus IDs found")
    }
    
    indices <- which( sapply(rownames(x@metadata[['loci']]), predicate) )
    
    x@metadata[['loci']] <- x@metadata[['loci']][indices, ]
    
    for ( i in 1:length(x) ) {
        x[[i]] <- x[[i]][indices]
        x@quality[[i]] <- x@quality[[i]][indices]
    }
    
    return(x)
}

# subsetByLocusID (S4) ---------------------------------------------------------
#' @rdname subsetByLocusID
setGeneric('subsetByLocusID', subsetByLocusID)

# DNAStringSet::subsetByLocusID ------------------------------------------------
#' @rdname subsetByLocusID
setMethod('subsetByLocusID', signature='DNAStringSet', 
    definition = subsetByLocusID.DNAStringSet)

# QualityScaledDNAStringSet::subsetByLocusID -----------------------------------
#' @rdname subsetByLocusID
setMethod('subsetByLocusID', signature='QualityScaledDNAStringSet', 
    definition = subsetByLocusID.QualityScaledDNAStringSet)

# End of subsetByLocusID.R #####################################################