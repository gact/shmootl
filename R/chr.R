# Start of chr.R ###############################################################

# formatChr --------------------------------------------------------------------
#' Format chromosome labels.
#'   
#' @param x Vector of chromosome labels.
#' @param prefix Chromosome prefix.
#' @param use.roman Use the Roman numeral form of the chromosome.
#'          
#' @return Vector of formatted chromosome labels.
#' 
#' @keywords internal
#' @rdname formatChr
formatChr <- function(x, prefix=c('', 'c', 'chr'), use.roman=TRUE) {
    
    stopifnot( is.vector(x) || is.factor(x) )
    prefix <- match.arg(prefix)
    stopifnot( isBOOL(use.roman) )
    
    if ( length(x) > 0 ) {
        
        # Ensure all chromosomes are normalised.
        unresolved <- is.na(x) | ! x %in% const$chrtab$seqids
        x[unresolved] <- normChr(x[unresolved])
        
        if (use.roman) {
            indices <- match(x, const$chrtab$seqids)
            x <- const$chrtab$seqnames[indices]
        }
        
        x <- paste0(prefix, x)
    }
    
    return(x)
}

# isNormChr --------------------------------------------------------------------
#' Test for normalised chromosome labels.
#' 
#' @param x Vector of chromosome labels.
#'      
#' @return Logical vector indicating which chromosome labels are normalised.
#' 
#' @keywords internal
#' @rdname isNormChr
isNormChr <-function(x) {
    stopifnot( is.vector(x) || is.factor(x) )
    return( ! is.na(x) & x %in% const$chrtab$seqids )
}

# normChr ----------------------------------------------------------------------
#' Normalise chromosome labels.
#'   
#' @param x Vector of chromosome labels.
#'          
#' @return Vector of normalised chromosome labels.
#' 
#' @keywords internal
#' @rdname normChr
normChr <- function(x) {
    
    stopifnot( is.vector(x) || is.factor(x) )
    
    if ( length(x) > 0 ) {
        
        # Get chromosomes as whitespace-stripped strings.
        chr.strings <- stripWhite( as.character(x) )
        
        # Strip any extraneous parts of chromosome name.
        m <- regexec(const$pattern$chromosome, chr.strings, ignore.case=TRUE) 
        matches <- regmatches(chr.strings, m)
        chr.strings <- sapply(matches, function(x) 
            if ( length(x) > 0 ) { x[2] } else { NA } )
        
        # Resolve any chromosome aliases.
        chr.upper <- toupper(chr.strings)
        alias.mask <- ! is.na(chr.upper) & chr.upper %in% names(const$alias2chrom)
        chr.strings[alias.mask] <- const$alias2chrom[ chr.upper[alias.mask] ]
        
        # Get chromosome strings as integers.
        chr.numbers <- suppressWarnings( as.integer(chr.strings) )
        
        # Get indices of known chromosome seqnames matching chromosome strings.
        chr.str.indices <- match(chr.strings, const$chrtab$seqnames)
        
        # Get indices of known chromosome seqids matching chromosome strings.
        chr.num.indices <- match(chr.numbers, as.integer(const$chrtab$seqids))
        
        # Get combined indices of known chromosomes matching those specified.
        indices <- rep(NA_integer_, length(x))
        for ( chr.indices in list(chr.str.indices, chr.num.indices) ) {
            indices[ ! is.na(chr.indices) ] <- chr.indices[ ! is.na(chr.indices) ]
        }
        
        # Check all elements were normalised.
        unresolved <- unique(x[ is.na(indices) ])
        if ( length(unresolved) > 0 ) {
            stop("cannot normalise chromosomes - '", toString(unresolved), "'")
        }
        
        res <- const$chrtab$seqids[indices]
    
    } else {
        
        res <- character()
    }
    
    return(res)
}

# orderChr ---------------------------------------------------------------------
#' Order chromosome labels. 
#'  
#' @param x Vector of chromosome labels.
#'          
#' @return Vector of ordered indices for the input chromosome labels.
#' 
#' @keywords internal
#' @rdname orderChr
orderChr <- function(x) {
    stopifnot( is.vector(x) || is.factor(x) )
    return( order( match(normChr(x), const$chrtab$seqids) ) )
}

# rankChr ----------------------------------------------------------------------
#' Rank chromosome labels. 
#'  
#' @param x Vector of chromosome labels.
#'          
#' @return Vector of ranks for the input chromosome labels.
#' 
#' @keywords internal
#' @rdname rankChr
rankChr <- function(x) {
    stopifnot( is.vector(x) || is.factor(x) )
    return( match(normChr(x), const$chrtab$seqids))
}

# sortChr ----------------------------------------------------------------------
#' Sort chromosome labels. 
#'  
#' @param x Vector of chromosome labels.
#'          
#' @return Sorted vector of chromosome labels.
#' 
#' @keywords internal
#' @rdname sortChr
sortChr <- function(x) {
    stopifnot( is.vector(x) || is.factor(x) )
    return( x[ order( match(normChr(x), const$chrtab$seqids) ) ] )
}

# End of chr.R #################################################################