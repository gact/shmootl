# Start of getQTLPeaks.R #######################################################

# getQTLPeaks ------------------------------------------------------------------
#' Get QTL peaks.
#'
#' @usage 
#' ## Generic method.
#' getQTLPeaks(x, chr=NULL, ...)
#'
#' @param x An \pkg{R/qtl} \code{scanone} or \code{qtl} object.
#' @param chr Vector indicating which sequences to consider. If no sequences
#' are specified, QTL peaks are returned for all available sequences.
#' @param ... Further arguments (see below).
#' @param threshold In a \code{scanone} or equivalent object, this indicates 
#' the LOD significance threshold for QTL peaks.
#' @template param-lodcolumn
#' @param qtl.indices In a \code{qtl} object, this option indicates the QTLs 
#' for which a QTL peak should be returned. 
#'   
#' @return A \code{mapframe} object in which each row indicates the locus of a 
#' QTL peak, with row names containing respective QTL names. Returns \code{NULL}
#' if there are no significant QTLs.
#'   
#' @export
#' @rdname getQTLPeaks
getQTLPeaks <- function(x, chr=NULL, ...) {
    UseMethod('getQTLPeaks', x)
}

# getQTLPeaks.mapframe ---------------------------------------------------------
#' @rdname getQTLPeaks
getQTLPeaks.mapframe <- function(x, chr=NULL, threshold=NULL, lodcolumn=NULL) {
    
    stopifnot( getMapUnit(x) == 'cM' )
    stopifnot( nrow(x) > 0 )
    threshold <- as.numeric(threshold)
    stopifnot( isSingleNonNegativeNumber(threshold) )

    seqcol.index <- getSeqColIndex(x)
    poscol.index <- getPosColIndex(x)
    lodcol.index <- getDatColIndices(x, datcolumns=lodcolumn)
    
    if ( length(lodcol.index) > 1 ) {
        stop("cannot get QTL peaks for multiple LOD columns - please choose one")
    } else if ( length(lodcol.index) == 0 ) {
        stop("no LOD column found")
    }

    # Get sequences corresponding to individual map loci.
    x.seqs <- pullLocusSeq(x)
    
    # Get map sequences.
    map.seqs <- unique(x.seqs)
    
    # Get specified sequences.
    chr <- subsetBySeq(map.seqs, chr)
    stopifnot( length(chr) > 0 )
    
    # Return NULL if no LOD data.
    if( allNA(x[, lodcol.index]) ) {
        return(NULL)
    }
    
    # Create LOD character mask. Regions with significant LOD values are marked
    # with the sequence name, while all other loci are marked with NA. 
    # NB: this is used instead of a simple logical mask so as to prevent
    # significant regions 'spilling over' across sequence boundaries.
    lod.mask <- vector('character', nrow(x))
    for ( chr.seq in chr ) {
        indices <- which( x.seqs == chr.seq )
        seq.lod <- x[indices, lodcol.index]
        bool.mask <- ! is.na(seq.lod) & seq.lod >= threshold
        char.mask <- sapply(bool.mask, function(sig.lod) 
            if (sig.lod) {chr.seq} else { NA })
        lod.mask[indices] <- char.mask
    }

    # Get list of vectors, each containing the row 
    # indices of the loci in a significant region.
    lod.runs <- rle(lod.mask)
    J <- cumsum(lod.runs$lengths)
    I <- c( 1, sapply(J[1:(length(J)-1)], function(j) j + 1) )
    I <- I[ ! is.na(lod.runs$values) ]
    J <- J[ ! is.na(lod.runs$values) ]
    significant.rows <- mapply(function(i, j) i:j, I, J, SIMPLIFY=FALSE)
   
    # Return NULL if no significant regions found.
    if ( length(significant.rows) == 0 ) {
        return(NULL)
    }
    
    # Get row index of peak for each significant 
    # region, taking first peak in case of ties.
    peak.indices <- unlist( lapply(significant.rows, function(rows) 
        rows[ which.max(x[rows, lodcol.index]) ]  ) )

    # Create genetic mapframe of QTL peaks.
    qtl.peaks <- gmapframe(chr=x[peak.indices, seqcol.index], 
        pos=x[peak.indices, poscol.index])
    
    # If mapframe has row names, get locus IDs at QTL 
    # peaks, identify unnamed loci (i.e. pseudomarkers)..
    if ( hasRownames(x) ) {
        peak.ids <- rownames(x)[peak.indices]
        unnamed <- which( isPseudomarkerID(rownames(x)[peak.indices]) )
    } else { # ..otherwise all loci are unnamed.
        peak.ids <- character( length(peak.indices) )
        unnamed <- 1:length(peak.indices)
    }
    
    # Ensure every QTL peak has a name.
    if ( length(unnamed) > 0 ) {
        peak.ids[unnamed] <- makeDefaultQTLNames( qtl.peaks[unnamed, ], 
            step=inferMapStep(x) )
    }
    
    # Set QTL peak names.
    rownames(qtl.peaks) <- peak.ids
    
    return(qtl.peaks)
}

# getQTLPeaks.qtl --------------------------------------------------------------
#' @rdname getQTLPeaks
getQTLPeaks.qtl <- function(x, chr=NULL, qtl.indices=NULL) {
    
    # NB: if you're thinking of assuming QTL-object 
    # sequences are in order, DON'T!
    
    # Get sequences corresponding to individual QTLs.
    norm.qtl.seqs <- normSeq(x$chr)
    
    # Get specified sequences.
    chr <- subsetBySeq(sortSeq( unique(norm.qtl.seqs) ), chr)
    stopifnot( length(chr) > 0 )
    
    qtl.indices <- resolveQtlIndices(x, qtl.indices)
    
    # Filter QTL indices by sequence. 
    qtl.indices <- which( norm.qtl.seqs %in% chr & 1:x$n.qtl %in% qtl.indices )
    
    if ( length(qtl.indices) == 0 ) {
        stop("QTL object has no QTLs on sequences '", toString(chr), "'")
    }
    
    # Get QTL peaks for the given sequence and QTL indices.
    qtl.peaks <- gmapframe(chr=norm.qtl.seqs[qtl.indices], 
        pos=x$pos[qtl.indices], row.names=x$name[qtl.indices])
        
    return(qtl.peaks)
}

# getQTLPeaks.scanone ----------------------------------------------------------
#' @rdname getQTLPeaks
getQTLPeaks.scanone <- function(x, chr=NULL, threshold=NULL, lodcolumn=NULL) {
    return( getQTLPeaks(as.mapframe(x), chr=chr, threshold=threshold, 
        lodcolumn=lodcolumn) )
}

# testQTLPeakSignificance ------------------------------------------------------
#' Test if QTL peaks are significant.
#' 
#' Given an object containing a LOD profile and a \code{mapframe} of QTL peaks, 
#' this function tests which peaks are significant with respect to the specified 
#' LOD threshold. 
#' 
#' @param x A \code{scanone} object.
#' @param threshold LOD significance threshold.
#' @param qtl.peaks Locus \code{mapframe} containing QTL peaks.
#' @param ... Further arguments. These are passed to \code{\link{getLODProfile}},
#' and may include a \code{lodcolumn} for a \code{scanone} object with multiple 
#' LOD columns, or a \code{qtl.index} for a \code{qtl} object with multiple QTLs.
#' 
#' @return Logical vector indicating which QTL peaks are significant.
#'   
#' @keywords internal
#' @rdname testQTLPeakSignificance
testQTLPeakSignificance <- function(x, threshold, qtl.peaks, ...) {
    stopifnot( isSingleNonNegativeNumber(threshold) )
    lod.profile <- getLODProfile(x, ...)
    peak.lods <- getLODValues(lod.profile, qtl.peaks)
    return( ! is.na(peak.lods) & peak.lods >= threshold )
}
    
# End of getQTLPeaks.R #########################################################