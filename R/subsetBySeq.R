# Start of subsetBySeq.R #######################################################

# subsetBySeq ------------------------------------------------------------------
#' Subset object by the specified sequences.
#' 
#' @param x Object subsettable by sequence.
#' @param sequences Sequences with which to subset object.
#'     
#' @return Input object subsetted by the specified sequences.
#' 
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

# subsetBySeq.qtlintervals -----------------------------------------------------
#' @rdname subsetBySeq
subsetBySeq.qtlintervals <- function(x, sequences=NULL) {
    
    if ( ! is.null(sequences) ) {
        
        sequences <- normSeq(sequences)
        
        interval.seqs <- sapply(x, function(obj) unique(obj[, 'chr']))
        
        for ( i in rev( getIndices(x) ) ) { # NB: must delete in reverse order
            if ( ! interval.seqs[i] %in% sequences ) {
                x[i] <- NULL
            }
        }
    }
    
    return(x)
}

# End of subsetBySeq.R #########################################################