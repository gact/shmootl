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
#' @include map.R
#' @rdname estPhysicalPositions
estPhysicalPositions <- function(qtl.intervals, map.key=NULL) {
    
    stopifnot( 'qtlintervals' %in% class(qtl.intervals) )
    
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

# qtlintervals -----------------------------------------------------------------
#' Create an empty \code{qtlintervals} object.
#' 
#' @param drop LOD units that the LOD profile must drop to form a QTL interval.
#' @param threshold LOD significance threshold for QTL intervals.
#' @param alpha Significance level of QTL intervals threshold. (Mutually
#' exclusive with \code{fdr}.)
#' @param fdr False discovery rate (FDR) of QTL intervals threshold. (Mutually
#' exclusive with \code{alpha}.)
#' @param ... Unused arguments.
#' 
#' @return An empty \code{qtlintervals} list with the given attributes.
#' 
#' @keywords internal
#' @rdname qtlintervals
qtlintervals <- function(drop=NULL, threshold=NULL, alpha=NULL, fdr=NULL, ...) {
    
    # TODO: create full qtlintervals object
    
    stopifnot( emptyArgs(...) )
    
    intervals <- list()
    
    if ( ! is.null(drop) ) {
        stopifnot( isSingleNonNegativeNumber(drop) )
        attr(intervals, 'drop') <- drop
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