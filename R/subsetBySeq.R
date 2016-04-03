# Start of subsetBySeq.R #######################################################

# subsetBySeq (S3) -------------------------------------------------------------
#' Subset object by the specified sequences.
#' 
#' @param x Object subsettable by sequence.
#' @param sequences Sequences with which to subset object.
#'     
#' @return Input object subsetted by the specified sequences.
#' 
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings QualityScaledDNAStringSet
#' @keywords internal
#' @rdname subsetBySeq
subsetBySeq <- function(x, sequences=NULL) {
    if ( ! is.null(sequences) && length(sequences) == 0 ) {
        stop("cannot subset object by sequence - no sequences specified")
    }
    UseMethod('subsetBySeq', x)
}

# subsetBySeq.character --------------------------------------------------------
#' @rdname subsetBySeq
subsetBySeq.character <- function(x, sequences=NULL) {
    
    if ( ! is.null(sequences) ) {
        sequences <- normSeq(sequences)
        x.seqs <- normSeq(x)
        x <- x[ x.seqs %in% sequences ]
    }
    
    return(x)
}

# subsetBySeq.data.frame -------------------------------------------------------
#' @rdname subsetBySeq
subsetBySeq.data.frame <- function(x, sequences=NULL) {
    
    if ( ! is.null(sequences) ) {
        indices <- matchSeqRowIndices(x, sequences, simplify=TRUE)
        x <- x[indices, ]
    }
    
    return(x)
}

# subsetBySeq.map --------------------------------------------------------------
#' @rdname subsetBySeq
subsetBySeq.map <- function(x, sequences=NULL) {
    
    if ( ! is.null(sequences) ) {
        sequences <- normSeq(sequences)
        norm.map.seqs <- normSeq( names(x) )
        x <- subsetMap(x, norm.map.seqs %in% sequences)
    }
    
    if ( length(x) < const$min.spm ) {
        stop("subsetted map has too few sequences (min=", const$min.spm, ")")
    }
    
    return(x)
}

# subsetBySeq.DNAStringSet -----------------------------------------------------
#' @rdname subsetBySeq
subsetBySeq.DNAStringSet <- function(x, sequences=NULL) {
    
    stopifnot( 'loci' %in% names(x@metadata) )
    stopifnot( 'mapframe' %in% class(x@metadata[['loci']]) )

    if ( ! is.null(sequences) ) {
        
        others <- otherattributes(x)
        meta <- x@metadata
        elem.meta <-x@elementMetadata
        
        indices <- matchSeqRowIndices(meta[['loci']], sequences, 
            simplify=TRUE)
        
        meta[['loci']] <- meta[['loci']][indices, ]
        
        bases <- lapply(Biostrings::as.list(x), function(b) b[indices])
        
        x <- Biostrings::DNAStringSet(bases)
        
        x@elementMetadata <- elem.meta
        x@metadata <- meta
        otherattributes(x) <- others         
    }
    
    return(x)
}

# subsetBySeq.QualityScaledDNAStringSet ----------------------------------------
#' @rdname subsetBySeq
subsetBySeq.QualityScaledDNAStringSet <- function(x, sequences=NULL) {
    
    stopifnot( 'loci' %in% names(x@metadata) )
    stopifnot( 'mapframe' %in% class(x@metadata[['loci']]) )
    
    if ( ! is.null(sequences) ) {

        others <- otherattributes(x)
        meta <- x@metadata
        elem.meta <-x@elementMetadata
        
        indices <- matchSeqRowIndices(meta[['loci']], sequences, 
            simplify=TRUE)
        
        meta[['loci']] <- meta[['loci']][indices, ]

        bases <- lapply(Biostrings::as.list(x), function(b) b[indices])
        quals <- lapply(Biostrings::as.list(x@quality), function(q) q[indices])
        
        x <- Biostrings::QualityScaledDNAStringSet( 
            Biostrings::DNAStringSet(bases), 
            Biostrings::PhredQuality( Biostrings::BStringSet(quals) )
        )
        
        x@elementMetadata <- elem.meta
        x@metadata <- meta
        otherattributes(x) <- others        
    }
    
    return(x)
}

# subsetBySeq (S4) -------------------------------------------------------------
#' @rdname subsetBySeq
setGeneric('subsetBySeq', subsetBySeq)

# DNAStringSet::subsetBySeq ----------------------------------------------------
#' @rdname subsetBySeq
setMethod('subsetBySeq', signature='DNAStringSet', 
    definition = subsetBySeq.DNAStringSet)

# QualityScaledDNAStringSet::subsetBySeq ---------------------------------------
#' @rdname subsetBySeq
setMethod('subsetBySeq', signature='QualityScaledDNAStringSet', 
    definition = subsetBySeq.QualityScaledDNAStringSet)

# End of subsetBySeq.R #########################################################