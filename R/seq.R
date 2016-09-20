# Start of seq.R ###############################################################

# formatSeq --------------------------------------------------------------------
#' Format sequence labels.
#'   
#' @param x Vector of sequence labels.
#' @param prefix Chromosome prefix.
#' @param use.roman Use the Roman numeral form of the chromosome.
#'          
#' @return Vector of formatted sequence labels.
#' 
#' @export
#' @family chromosome/sequence functions
#' @rdname formatSeq
formatSeq <- function(x, prefix=c('', 'c', 'chr'), use.roman=TRUE) {
    
    if ( length(x) > 0 ) {
        lu <- splitSeqLabels(x)
        lu$chr <- formatChr(lu$chr, prefix=prefix, use.roman=use.roman)
        x <- joinSeqLabels(lu)
    }
    
    return(x)
}

# isNormSeq --------------------------------------------------------------------
#' Test for normalised sequence labels.
#' 
#' @param x Vector of sequence labels.
#'      
#' @return Logical vector indicating which sequence labels are normalised.
#' 
#' @export
#' @family chromosome/sequence functions
#' @rdname isNormSeq
isNormSeq <-function(x) {
    stopifnot( is.vector(x) || is.factor(x) )
    genome <- genomeOpt()
    return( ! is.na(x) & x %in% const$seqtab[[genome]]$seqids )
}

# joinSeqLabels ----------------------------------------------------------------
#' Join sequence labels.
#'   
#' @param x A \code{data.frame} representing split sequence labels.
#'      
#' @return Character vector of joined sequence labels.
#'  
#' @keywords internal
#' @rdname joinSeqLabels
joinSeqLabels <- function(x) {
    
    stopifnot( is.data.frame(x) )
    stopifnot( colnames(x) == c('chr', 'seq') )
    stopifnot( sapply(x, class) == 'character' )
    
    res <- rep(NA_character_, nrow(x))
    
    if ( nrow(x) > 0 ) {
        
        has.chr <- ! is.na(x$chr)
        is.seq <- has.chr & ! is.na(x$seq)
        is.chr <- has.chr & ! is.seq
        
        res[is.seq] <- paste(x$chr[is.seq], x$seq[is.seq], sep='_')
        res[is.chr] <- x$chr[is.chr]
    } 
    
    return(res)
}

# normSeq ----------------------------------------------------------------------
#' Normalise sequence labels.
#' 
#' @param x Object containing sequence labels.
#' 
#' @return Input object with normalised sequence labels.
#' 
#' @template section-chr-seq
#' 
#' @export
#' @family chromosome/sequence functions
#' @rdname normSeq
normSeq <- function(x) {
    UseMethod('normSeq', x)
}

# normSeq.character ------------------------------------------------------------
#' @export
#' @rdname normSeq
normSeq.character <- function(x) {
    
    if ( length(x) > 0 ) {
    
        lu <- splitSeqLabels(x)
        lu$chr <- normChr(lu$chr)
        res <- joinSeqLabels(lu)
        
        genome <- genomeOpt()
        
        res[ ! res %in% const$seqtab[[genome]]$seqids ] <- NA
        
        # Check all elements were normalised to a known sequence.
        unresolved <- unique(x[ is.na(res) ])
        if ( length(unresolved) > 0 ) {
            stop("cannot normalise sequence labels - '", toString(unresolved), "'")
        }
        
        x <- res
    } 
    
    return(x)
}

# normSeq.data.frame -----------------------------------------------------------
#' @export
#' @rdname normSeq
normSeq.data.frame <- function(x) {
    
    if ( nrow(x) > 0 ) {
        seqcol.index <- getSeqColIndex(x)
        x[, seqcol.index] <- normSeq(x[, seqcol.index])
    }
    
    return(x)
}

# normSeq.factor ---------------------------------------------------------------
#' @export
#' @rdname normSeq
normSeq.factor <- function(x) {
    return( normSeq( as.character(x) ) )
}

# normSeq.integer --------------------------------------------------------------
#' @export
#' @rdname normSeq
normSeq.integer <- function(x) {
    return( normSeq( as.character(x) ) )
}

# normSeq.map ------------------------------------------------------------------
#' @export
#' @rdname normSeq
normSeq.map <- function(x) {
    
    names(x) <- normSeq( names(x) )
    
    if( anyDuplicated( names(x) ) ) {
        stop("map has inconsistent sequence labels")
    }
    
    return(x)
}

# normSeq.numeric --------------------------------------------------------------
#' @export
#' @rdname normSeq
normSeq.numeric <- function(x) {
    return( normSeq( as.character(x) ) )
}

# orderSeq ---------------------------------------------------------------------
#' Order sequence labels. 
#'  
#' @param x Vector of sequence labels.
#'          
#' @return Vector of indices for the input sequence labels,
#' ordered with respect to their normalised form.
#' 
#' @export
#' @family chromosome/sequence functions
#' @rdname orderSeq
orderSeq <- function(x) {
    stopifnot( is.vector(x) || is.factor(x) )
    genome <- genomeOpt()
    return( order( match(normChr(x), const$seqtab[[genome]]$seqids) ) )
}

# rankSeq ----------------------------------------------------------------------
#' Rank sequence labels. 
#'  
#' @param x Vector of sequence labels.
#'          
#' @return Vector of ranks for the normalised
#' form of the input sequence labels.
#' 
#' @export
#' @family chromosome/sequence functions
#' @rdname rankSeq
rankSeq <- function(x) {
    stopifnot( is.vector(x) || is.factor(x) )
    genome <- genomeOpt()
    return( match(normChr(x), const$seqtab[[genome]]$seqids))
}

# sortSeq ----------------------------------------------------------------------
#' Sort sequence labels. 
#'  
#' @param x Vector of sequence labels.
#'          
#' @return Input vector of sequence labels, sorted
#' with respect to their normalised form.
#' 
#' @export
#' @family chromosome/sequence functions
#' @rdname sortSeq
sortSeq <- function(x) {
    stopifnot( is.vector(x) || is.factor(x) )
    genome <- genomeOpt()
    return( x[ order( match(normChr(x), const$seqtab[[genome]]$seqids) ) ] )
}

# splitSeqLabels ---------------------------------------------------------------
#' Split sequence labels.
#'   
#' @param x Vector of sequence labels.
#'      
#' @return A \code{data.frame} representing split sequence labels.
#' 
#' @keywords internal
#' @rdname splitSeqLabels
splitSeqLabels <- function(x) {
    
    stopifnot( is.vector(x) || is.factor(x) )
    
    res <- data.frame(chr=character( length(x) ), 
        seq=character( length(x) ), stringsAsFactors=FALSE)
    
    if ( length(x) > 0 ) {
    
        # Get sequences as whitespace-stripped strings.
        strip.strings <- stripWhite( as.character(x) )
        
        # Split by underscore.
        split.strings <- strsplit(strip.strings, '_', fixed=TRUE)
        
        unsplit <- unique(x[ ! lengths(split.strings) %in% 1:2 ])
        if ( length(unsplit) > 0 ) {
            stop("cannot split sequence labels - '", toString(unsplit), "'")
        }
        
        res$chr <- sapply(split.strings, function(x) 
            if ( length(x) >= 1 ) { x[1] } else { NA_character_ })
        
        res$seq <- sapply(split.strings, function(x) 
            if ( length(x) == 2 ) { x[2] } else { NA_character_ })
    }
    
    return(res)
}

# End of seq.R #################################################################