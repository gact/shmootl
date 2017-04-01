# Start of getQTLIntervals.R ###################################################

# getQTLIntervals --------------------------------------------------------------
#' Get QTL intervals.
#' 
#' @description Given a specified LOD \code{threshold} and associated
#' significance level \code{alpha} (or alternatively, false-discovery rate
#' \code{fdr}), this function gets approximate QTL intervals for a given
#' \code{scanone} object, or a \code{qtl} object with a \code{'lodprofile'}
#' attribute (as output by \pkg{R/qtl} function \code{refineqtl}).
#' 
#' @usage
#' ## Generic method.
#' getQTLIntervals(x, chr=NULL, ci.function=c('lodint', 'bayesint'), drop=1.5,
#'     prob=0.95, expandtomarkers=FALSE, ...)
#' 
#' @param x A \code{scanone} or \code{qtl} object.
#' @param chr Vector indicating which sequences to consider. If no sequences
#' are specified, QTL intervals are returned for all available sequences.
#' @param ci.function Option to indicate which function should be used for
#' estimating approximate confidence intervals for QTL location. Set to
#' \code{'lodint'} for LOD support intervals (adjusting stringency with the
#' \code{drop} parameter), or to \code{'bayesint'} for Bayesian credible
#' intervals (adjusting stringency with the \code{prob} parameter). For more
#' information on the QTL interval methods used, see functions \code{lodint} and
#' \code{bayesint} in the \pkg{R/qtl} manual, as well as Section 4.5 of Broman
#' and Sen (2009).
#' @param drop LOD units that the LOD profile must drop to form the interval.
#' This is used only if \code{ci.function} is set to \code{'lodint'}.
#' @param prob Desired probability coverage for the Bayesian credible interval.
#' This is used only if \code{ci.function} is set to \code{'bayesint'}.
#' @param expandtomarkers Expand the LOD interval to the nearest flanking
#' markers, or to the respective terminal loci.
#' @param ... Further arguments (see below).
#' @param threshold For a \code{scanone} or equivalent object, this
#' contains a single \code{numeric} LOD significance threshold, or
#' an object (e.g. \code{summary.scanoneperm}) containing one such
#' threshold and its associated significance level.
#' @template param-lodcolumn
#' @param qtl.peaks Locus \code{mapframe} indicating the location of the QTL
#' peaks in a \code{scanone} result. If not specified, these are inferred
#' from the LOD profile of the \code{scanone} object.
#' @param qtl.indices In a \code{qtl} object, this option indicates the QTLs
#' for which a QTL interval should be returned.
#'   
#' @return An object of class \code{qtlintervals}, which is essentially a list
#' of \code{data.frame} objects, each containing three rows of information about
#' the lower QTL interval limit, QTL peak, and upper QTL interval limit,
#' respectively. Returns an empty \code{qtlintervals} object if there are no
#' significant QTLs.
#' 
#' @template author-thomas-walsh
#' @template author-yue-hu
#' @template ref-broman-2003
#' @template ref-broman-2009
#' @template seealso-rqtl-manual
#' 
#' @export
#' @family QTL functions
#' @rdname getQTLIntervals
getQTLIntervals <- function(x, chr=NULL, ci.function=c('lodint', 'bayesint'),
    drop=1.5, prob=0.95, expandtomarkers=FALSE, ...) {
    UseMethod('getQTLIntervals', x)
}

# getQTLIntervals.mapframe -----------------------------------------------------
#' @export
#' @rdname getQTLIntervals
getQTLIntervals.mapframe <- function(x, chr=NULL, ci.function=c('lodint',
    'bayesint'), drop=1.5, prob=0.95, expandtomarkers=FALSE, threshold=NULL,
    lodcolumn=NULL, qtl.peaks=NULL, ...) {

    stopifnot( getMapUnit(x) == 'cM' )
    stopifnot( nrow(x) > 0 )
    stopifnot( isBOOL(expandtomarkers) )
    stopif( is.null(threshold) )
    
    ci.function <- match.arg(ci.function)
    
    # Ensure only relevant parameter is stored in `qtlintervals` object.
    if ( ci.function == 'lodint' ) {
        prob <- NULL
    } else { # ci.function == 'bayesint'
        drop <- NULL
    }
    
    # Init `qtlintervals` object.
    intervals <- qtlintervals(threshold=threshold, drop=drop, prob=prob)
    
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
        }
        
    } else { # ..otherwise find QTL peaks from LOD profile.
        
        qtl.peaks <- getQTLPeaks(x, chr=chr, threshold=threshold)
    }
    
    # Return empty intervals list if there are no significant peaks.
    if ( nrow(qtl.peaks) == 0 ) {
        return(intervals)
    }

    # Get row indices of QTL peaks.
    peak.indices <- findLocusRowIndices(x, qtl.peaks)
    
    # Init matrix of raw intervals, each row of which will contain the row 
    # indices of the lower limit, peak, and upper limit of a significant region.
    raw.intervals <- matrix( nrow=length(peak.indices), ncol=3, 
        dimnames=list(NULL, names(ci)) )
    
    # Estimate LOD support intervals, if specified..
    if ( ci.function == 'lodint' ) {
        
        for ( i in seq_along(peak.indices) ) {
            
            # Init interval indices to peak index.
            l <- u <- peak.index <- peak.indices[i]
            
            # Get sequence endpoint indices.
            peak.seq <- x.seqs[peak.index]
            seq.start <- index.ranges[[peak.seq]][1]
            seq.end <- index.ranges[[peak.seq]][2]
            
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
    
    # ..otherwise estimate Bayesian credible intervals.
    } else { # ci.function == 'bayesint'
        
        for ( i in seq_along(peak.indices) ) {
            
            # Init interval indices to peak index.
            l <- u <- peak.index <- peak.indices[i]
            
            # Get sequence endpoint indices.
            peak.seq <- x.seqs[peak.index]
            seq.start <- index.ranges[[peak.seq]][1]
            seq.end <- index.ranges[[peak.seq]][2]
            
            # Call R/qtl function to do Bayesian credible interval.
            bayes.interval <- qtl::bayesint(results=x, chr=peak.seq,
                qtl.index=1, prob=prob, expandtomarkers=FALSE)
            
            # Decrement lower limit within this sequence until position matches
            # Bayesian interval lower endpoint (within numeric tolerance).
            while ( ! isTRUE( all.equal(x$lod[l], bayes.interval$pos[1]) ) &&
                l > seq.start && ! is.na(x$lod[l-1]) ) {
                l <- l - 1
            }
            
            # Increment upper limit within this sequence until position matches
            # Bayesian interval upper endpoint (within numeric tolerance).
            while ( ! isTRUE( all.equal(x$lod[u], bayes.interval$pos[3]) ) &&
                u < seq.end && ! is.na(x$lod[u+1]) ) {
                u <- u + 1
            }
            
            raw.intervals[i, ] <- c(l, peak.index, u)
        }
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
    
    for ( r in getRowIndices(merged.intervals) ) {
        
        merged.interval <- merged.intervals[r, ]
        
        interval <- as.data.frame(x[merged.interval, ])
        colnames(interval)[2] <- 'pos (cM)'
        
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
        
        intervals[[r]] <- interval
    }

    # Get interval peaks.
    interval.peaks <- do.call( rbind, lapply(intervals,
        function(interval) interval[ ci['peak'], ]) )
     
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
    
    return(intervals) 
}

# getQTLIntervals.qtl ----------------------------------------------------------
#' @export
#' @rdname getQTLIntervals
getQTLIntervals.qtl <- function(x, chr=NULL, ci.function=c('lodint', 'bayesint'),
    drop=1.5, prob=0.95, expandtomarkers=FALSE, qtl.indices=NULL, ...) {
    
    ci.function <- match.arg(ci.function)
    
    # NB: if you're thinking of assuming QTL-object 
    # sequences are in order, DON'T!
    
    # Get sequences corresponding to individual QTLs.
    norm.qtl.seqs <- normSeq(x$chr)
    
    # Get specified sequences.
    chr <- subsetBySeq(sortSeq( unique(norm.qtl.seqs) ), chr)
    stopifnot( length(chr) > 0 )
    
    # Filter QTL indices by sequence.
    seq.mask <- norm.qtl.seqs %in% chr
    qtl.mask <- getMask(x, requested=qtl.indices)
    qtl.indices <- which( seq.mask & qtl.mask )
    
    if ( length(qtl.indices) == 0 ) {
        stop("QTL object has no QTLs on sequences '", toString(chr), "'")
    }
    
    # Get list of QTL peak mapframes for given QTL indices.
    qtl.peaks <- lapply(qtl.indices, function(i) gmapframe(chr=norm.qtl.seqs[i], 
        pos=x$pos[i], row.names=x$name[i]))
    
    # Get QTL intervals corresponding to QTL peaks.
    intervals <- list( length(qtl.peaks) )
    for ( i in seq_along(qtl.peaks) ) {
        lod.profile <- getLODProfile(x, qtl.index=qtl.indices[i])
        interval.result <- getQTLIntervals(lod.profile, ci.function=ci.function,
            drop=drop, prob=prob, expandtomarkers=expandtomarkers,
            qtl.peaks=qtl.peaks[[i]])
        intervals[[i]] <- interval.result[[1]]
        names(intervals)[i] <- names(interval.result)[1]
    }
    
    class(intervals) <- c('qtlintervals', 'list')
    
    if ( ci.function == 'lodint' ){
        attr(intervals, 'drop') <- drop
    } else { # ci.function == 'bayesint'
        attr(intervals, 'prob') <- prob
    }
    
    return(intervals)
}

# getQTLIntervals.scanone ------------------------------------------------------
#' @export
#' @rdname getQTLIntervals
getQTLIntervals.scanone <- function(x, chr=NULL, ci.function=c('lodint',
    'bayesint'), drop=1.5, prob=0.95, expandtomarkers=FALSE, threshold=NULL,
    lodcolumn=NULL, qtl.peaks=NULL, ...) {
    return( getQTLIntervals(as.mapframe(x), chr=chr, ci.function=ci.function,
        drop=drop, prob=prob, expandtomarkers=expandtomarkers,
        threshold=threshold, lodcolumn=lodcolumn, qtl.peaks=qtl.peaks) )
}
    
# End of getQTLIntervals.R #####################################################