# Start of qtlIntervals.R ######################################################

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