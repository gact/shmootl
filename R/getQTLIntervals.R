# Start of getQTLIntervals.R ###################################################

# TODO: implement optional Bayesian credible interval.

# getQTLIntervals --------------------------------------------------------------
#' Get QTL intervals.
#'
#' @usage 
#' ## Generic method.
#' getQTLIntervals(x, chr=NULL, drop=1.5, expandtomarkers=FALSE, ...)
#'
#' @param x A \code{scanone} or \code{qtl} object.
#' @param chr Vector indicating which sequences to consider. If no sequences
#' are specified, QTL intervals are returned for all available sequences.
#' @param drop LOD units that the LOD profile must drop to form the interval.
#' @param expandtomarkers Expand the LOD interval to the nearest flanking 
#' markers, or to the respective terminal loci.
#' @param ... Further arguments (see below).
#' @param threshold In a \code{scanone} or equivalent object, this indicates 
#' the LOD significance threshold for QTL intervals.
#' @template param-lodcolumn
#' @param qtl.peaks Locus \code{mapframe} indicating the location of the QTL 
#' peaks in a \code{scanone} result. If not specified, these are found using 
#' the LOD profile of the \code{scanone} object.
#' @param qtl.indices In a \code{qtl} object, this option indicates the QTLs 
#' for which a QTL interval should be returned. 
#'   
#' @return An object of class \code{qtlintervals}, which is essentially a list
#' of \code{mapframe} objects, each containing three rows of information about
#' the lower QTL interval limit, QTL peak, and upper QTL interval limit, 
#' respectively. Returns \code{NULL} if there are no significant QTLs.
#'  
#' @export
#' @rdname getQTLIntervals
getQTLIntervals <- function(x, chr=NULL, drop=1.5, expandtomarkers=FALSE, ...) {
    UseMethod('getQTLIntervals', x)
}

# getQTLIntervals.mapframe -----------------------------------------------------
#' @export
#' @rdname getQTLIntervals
getQTLIntervals.mapframe <- function(x, chr=NULL, drop=1.5, expandtomarkers=FALSE, 
    threshold=NULL, lodcolumn=NULL, qtl.peaks=NULL, ...) {

    stopifnot( getMapUnit(x) == 'cM' )
    stopifnot( nrow(x) > 0 )
    stopifnot( isSingleNonNegativeNumber(threshold) )
    getThresholdAlpha(threshold) # validate threshold alpha
    stopifnot( isSingleNonNegativeNumber(drop) )
    stopifnot( isBOOL(expandtomarkers) )
    
    # Set standard QTL interval indices.
    ci <- c(low=1, peak=2, high=3)
    
    # Get LOD profile for the given LOD column.
    x <- getLODProfile(x, lodcolumn=lodcolumn)
    
    # Get sequences corresponding to individual map loci.
    x.seqs <- pullLocusSeq(x)
    
    # Get map sequences.
    map.seqs <- unique(x.seqs)
    
    # Get specified sequences.
    chr <- subsetBySeq(map.seqs, chr)
    stopifnot( length(chr) > 0 )
    
    # Get list of sequence ranges.
    index.ranges <- lapply(chr, function (chr.seq) 
        range( which( x.seqs == chr.seq ) ) )
    names(index.ranges) <- chr
    
    # If QTL peaks specified, filter by sequence and test significance..
    if ( ! is.null(qtl.peaks) ) {
        
        stopifnot( 'mapframe' %in% qtl.peaks )
        qtl.peaks <- subsetBySeq(qtl.peaks, chr)
        
        if ( nrow(qtl.peaks) > 0 ) {
            
            sig.peaks <- testQTLPeakSignificance(x, threshold, qtl.peaks)
            subthreshold <- rownames(qtl.peaks)[ ! sig.peaks ]
            if ( length(subthreshold) > 0 ) { 
                warning("QTL peaks below threshold - '", toString(subthreshold), "'")
            }
            
        } else {
            
            qtl.peaks <- NULL
        }
        
    } else { # ..otherwise find QTL peaks from LOD profile.
        
        qtl.peaks <- getQTLPeaks(x, chr=chr, threshold=threshold)
    }
    
    # Return NULL if there are no significant peaks.
    if ( is.null(qtl.peaks) ) {
        return(NULL)
    }

    # Get row indices of QTL peaks.
    peak.indices <- findLocusRowIndices(x, qtl.peaks)
    
    # Init matrix of raw intervals, each row of which will contain the row 
    # indices of the lower limit, peak, and upper limit of a significant region.
    raw.intervals <- matrix( nrow=length(peak.indices), ncol=3, 
        dimnames=list(NULL, names(ci)) )

    for ( i in getIndices(peak.indices) ) {
        
        # Init interval indices to peak index.
        l <- u <- peak.index <- peak.indices[i]
        
        # Get sequence endpoint indices.
        s <- as.character(x.seqs[peak.index])
        seq.start <- index.ranges[[s]][1]
        seq.end <- index.ranges[[s]][2]
        
        # Set LOD cutoff.
        cutoff <- max( x$lod[peak.index] - drop, 0 )
        
        # Decrement lower limit within this sequence while LOD score meets cutoff.
        while ( x$lod[l] >= cutoff && l > seq.start && ! is.na(x$lod[l-1]) ) { 
            l <- l - 1 
        }
        
        # Increment upper limit within this sequence while LOD score meets cutoff.
        while ( x$lod[u] >= cutoff && u < seq.end && ! is.na(x$lod[u+1]) ) { 
            u <- u + 1 
        }
        
        raw.intervals[i, ] <- c(l, peak.index, u)
    }

    # Create matrix of merged intervals, formed by combining overlapping raw intervals.
    merged.intervals <- matrix( nrow=0, ncol=3, dimnames=list(NULL, names(ci)) )
    i <- 1
    while ( i <= nrow(raw.intervals) ) {
        
        # Merged range will include this interval.
        j <- i
        
        # Add any intervals whose limits overlap the merged range.
        while ( j < nrow(raw.intervals) && raw.intervals[j+1, 'low'] <= raw.intervals[j, 'high'] ) {
            j <- j + 1
        }
        
        # Set the lower limit to that of the first interval.
        l <- raw.intervals[i, 'low']
        
        # Set the merged interval peak to the highest peak across the  
        # merged intervals, taking the first peak in case of ties.
        peak.indices <- raw.intervals[i:j, 'peak']
        peak.index <- peak.indices[ which.max( x$lod[peak.indices] ) ]
        
        # Set the upper limit to that of the last interval.
        u <- raw.intervals[j, 'high']
        
        # Add to merged intervals.
        merged.intervals <- rbind(merged.intervals, c(l, peak.index, u))
        
        i <- j + 1 
    }

    # If expanding to markers, expand each interval limit 
    # to next flanking marker or sequence endpoint.
    if (expandtomarkers) {
        
        # Create mask indicating which loci are markers.
        if ( hasRownames(x) ) {
            marker.mask <- ! isPseudomarkerID( rownames(x) )
        } else {
            marker.mask <- rep_len(FALSE, nrow(x))
        }
        
        # Expand each interval to next flanking marker or sequence endpoint.
        for ( r in getRowIndices(merged.intervals) ) {
            
            # Get interval limit indices.
            l <- merged.intervals[r, 'low']
            u <- merged.intervals[r, 'high']
            
            # Get sequence endpoint indices.
            s <- as.character(x.seqs[l]) # NB: sequence identical across interval.
            seq.start <- index.ranges[[s]][1]
            seq.end <- index.ranges[[s]][2]
            
            # Decrement lower limit within this sequence while locus is not a marker.
            while ( ! marker.mask[l] && l > seq.start && ! is.na(x$lod[l-1]) ) { 
                l <- l - 1 
            }
            
            # Increment upper limit within this sequence while locus is not a marker.
            while ( ! marker.mask[u] && u < seq.end && ! is.na(x$lod[u+1]) ) { 
                u <- u + 1 
            }
            
            # Set expanded interval limit indices.
            merged.intervals[r, 'low'] <- l
            merged.intervals[r, 'high'] <- u
        }
    }    

    # Create list of interval mapframes.
    intervals <- list()
    
    for ( r in getRowIndices(merged.intervals) ) {
        
        merged.interval <- merged.intervals[r, ]
        
        interval <- as.mapframe(x[merged.interval, ], map.unit='cM')
        
        # If peak coincides with start of interval, 
        # rename start of interval to 'ci.low'.
        if ( merged.interval['low'] == merged.interval['peak'] ) {
            temp.id <- rownames(interval)[ ci['low'] ]
            rownames(interval)[ ci['low'] ] <- 'ci.low'
            rownames(interval)[ ci['peak'] ] <- temp.id
        } 
        
        # If peak coincides with end of interval, 
        # rename end of interval to 'ci.high'.
        if ( merged.interval['peak'] == merged.interval['high'] ) {
            rownames(interval)[ ci['high'] ] <- 'ci.high'
        }    
        
        attr(interval, 'threshold') <- threshold
        attr(interval, 'drop') <- drop
        
        intervals[[r]] <- interval
    }

    # Get interval peaks.
    interval.peaks <- as.mapframe( do.call( rbind, lapply(intervals, 
        function(interval) interval[ ci['peak'], ]) ), map.unit='cM' )
     
    # Get locus IDs at interval peaks.
    peak.ids <- rownames(interval.peaks)
    
    # Identify unnamed interval peaks.
    unnamed <- which( isPseudomarkerID(peak.ids) )
   
    # Ensure every QTL interval has a name.
    if ( length(unnamed) > 0 ) {
        peak.ids[unnamed] <- makeDefaultQTLNames( interval.peaks[unnamed, ], 
            step=inferMapStep(x) )
    }

    # Set QTL interval names.
    names(intervals) <- peak.ids
    
    class(intervals) <- c('qtlintervals', 'list')
    
    attr(intervals, 'threshold') <- threshold
    attr(intervals, 'drop') <- drop
    
    return(intervals) 
}

# getQTLIntervals.qtl ----------------------------------------------------------
#' @export
#' @rdname getQTLIntervals
getQTLIntervals.qtl <- function(x, chr=NULL, drop=1.5, expandtomarkers=FALSE, 
    qtl.indices=NULL, ...) {
    
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
    
    # Get list of QTL peak mapframes for given QTL indices.
    qtl.peaks <- lapply(qtl.indices, function(i) gmapframe(chr=norm.qtl.seqs[i], 
        pos=x$pos[i], row.names=x$name[i]))
    
    # Get QTL intervals corresponding to QTL peaks.
    intervals <- list( length(qtl.peaks) )
    for ( i in getIndices(qtl.peaks) ) {
        lod.profile <- getLODProfile(x, qtl.index=qtl.indices[i])
        interval.result <- getQTLIntervals(lod.profile, drop=drop, 
            expandtomarkers=expandtomarkers, qtl.peaks=qtl.peaks[[i]])
        intervals[[i]] <- interval.result[[1]]
        names(intervals)[i] <- names(interval.result)[1]
    }
    
    class(intervals) <- c('qtlintervals', 'list')
    
    attr(intervals, 'drop') <- drop
    
    return(intervals)
}

# getQTLIntervals.scanone ------------------------------------------------------
#' @export
#' @rdname getQTLIntervals
getQTLIntervals.scanone <- function(x, chr=NULL, drop=1.5, expandtomarkers=FALSE,
    threshold=NULL, lodcolumn=NULL, qtl.peaks=NULL, ...) {
    return( getQTLIntervals(as.mapframe(x), chr=chr, threshold=threshold, drop=drop, 
        expandtomarkers=expandtomarkers, lodcolumn=lodcolumn, qtl.peaks=qtl.peaks) )
}
    
# End of getQTLIntervals.R #####################################################