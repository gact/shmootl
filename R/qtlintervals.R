# Start of qtlIntervals.R ######################################################

# estPhysicalPositions ---------------------------------------------------------
#' Add estimated physical positions to QTL intervals.
#' 
#' @param qtl.intervals  A \code{qtlintervals} object.
#' @param map.key An optional \code{mapkey} function for converting genetic map
#' positions to physical map positions. If none is specified, the default
#' \code{mapkey} for the current genome is used to estimate physical positions.
#' 
#' @return A copy of the input \code{qtlintervals} object, with an additional
#' column - \code{'pos (bp)'} - containing physical map positions of the QTL
#' intervals, as estimated from their genetic map positions.
#' 
#' @export
#' @family QTL functions
#' @rdname estPhysicalPositions
estPhysicalPositions <- function(qtl.intervals, map.key=NULL) {
    
    # Remove any existing physical locus positions.
    if ( hasPhysicalPositions(qtl.intervals) ) {
        
        for ( i in seq_along(qtl.intervals) ) {
            qtl.intervals[[i]] <- deleteColumn(qtl.intervals[[i]],
                col.name='pos (bp)')
        }
    }
    
    # Set mapkey option within this function, if specified.
    if ( ! is.null(map.key) ) {
        prev.mapkey <- mapkeyOpt(map.key)
        on.exit( mapkeyOpt(prev.mapkey) )
    }
    
    # Get copy of QTL intervals rescaled in physical (base-pair) units.
    physical.intervals <- lapply(qtl.intervals, setMapUnit, 'bp')
    
    # Insert physical locus positions in each QTL interval.
    for ( i in seq_along(qtl.intervals) ) {
        
        physical.positions <- pullLocusPos(physical.intervals[[i]])
        
        qtl.intervals[[i]] <- insertColumn(qtl.intervals[[i]], col.index=3,
            col.name='pos (bp)', data=physical.positions)
    }
    
    return(qtl.intervals)
}

# hasPhysicalPositions ---------------------------------------------------------
#' Test if QTL intervals contain physical positions.
#' 
#' @param qtl.intervals A \code{qtlintervals} object.
#' 
#' @return \code{TRUE} if \code{qtlintervals} object contains physical
#' positions in terms of base-pair units; \code{FALSE} otherwise.
#' 
#' @export
#' @rdname hasPhysicalPositions
hasPhysicalPositions <- function(qtl.intervals) {
    
    stopifnot( 'qtlintervals' %in% class(qtl.intervals) )
    
    poscol.index.list <- lapply(qtl.intervals, function(qtl.interval)
        which( colnames(qtl.interval) == 'pos (bp)' ) )
    
    poscol.index.lengths <- lengths(poscol.index.list)
    
    if ( any( poscol.index.lengths > 1 ) ) {
        stop("QTL intervals have multiple physical position columns")
    }
    
    if ( all( poscol.index.lengths == 1 ) ) {
        
        if ( all( unlist(poscol.index.list) == 3 ) ) {
            
            status <- TRUE
            
        } else {
            
            stop("QTL intervals have physical positions in unexpected columns")
        }
        
    } else if ( any( poscol.index.lengths == 1 )) {
        
        stop("QTL intervals have incomplete physical position data")
        
    } else {
        
        status <- FALSE
    }
    
    return(status)
}

# qtlintervals -----------------------------------------------------------------
#' Create a \code{qtlintervals} object.
#' 
#' This function creates a \code{qtlintervals} object. Currently it
#' can only create a \code{qtlintervals} object of length zero.
#' 
#' @param drop LOD units that the LOD profile must drop to form a QTL interval.
#' This should only be specified if the QTL intervals are LOD support intervals,
#' and is incompatible with the \code{prob} parameter.
#' @param prob The probability coverage of the Bayesian credible interval.
#' This should only be specified if the QTL intervals are Bayesian credible
#' intervals, and is incompatible with the \code{drop} parameter.
#' @param threshold LOD significance threshold for QTL intervals.
#' @param alpha Significance level of QTL intervals threshold. (Mutually
#' exclusive with \code{fdr}.)
#' @param fdr False discovery rate (FDR) of QTL intervals threshold. (Mutually
#' exclusive with \code{alpha}.)
#' @param ... Unused arguments.
#' 
#' @return An empty \code{qtlintervals} list with the given attributes.
#' 
#' @template author-thomas-walsh
#' @template author-yue-hu
#' 
#' @keywords internal
#' @rdname qtlintervals
qtlintervals <- function(drop=NULL, prob=NULL, threshold=NULL, alpha=NULL,
    fdr=NULL, ...) {
    
    # TODO: create full qtlintervals object
    
    stopifnot( emptyArgs(...) )
    
    intervals <- list()
    
    if ( ! is.null(drop) && ! is.null(prob) ) {
        stop("cannot set both LOD support interval drop",
            "and Bayesian credible interval probability")
    } else if ( ! is.null(drop) ) { # LOD support interval
        stopifnot( isSingleNonNegativeNumber(drop) )
        attr(intervals, 'drop') <- drop
    } else if ( ! is.null(prob) ) { # Bayes credible interval
        stopifnot( isSingleProbability(prob) )
        attr(intervals, 'prob') <- prob
    }
    
    if ( ! is.null(threshold) ) {
        stopifnot( isSingleNonNegativeNumber(threshold) )
        attr(intervals, 'threshold') <- threshold
    }    
    
    if ( ! is.null(alpha) && ! is.null(fdr) ) {
        stop("cannot set both significance level (alpha) and FDR")
    } else if ( ! is.null(alpha) ) {
        stopifnot( isSingleProbability(alpha) )
        attr(intervals, 'alpha') <- unname(alpha)
    } else if ( ! is.null(fdr) ) {
        stopifnot( isSingleFiniteNumber(fdr) )
        stopifnot( fdr > 0 & fdr < 1 )
        attr(intervals, 'fdr') <- unname(fdr)
    }
    
    class(intervals) <- c('qtlintervals', 'list')
    
    return(intervals)
}

# End of qtlIntervals.R ########################################################