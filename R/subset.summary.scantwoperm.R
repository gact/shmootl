# Start of subset.summary.scantwoperm.R ########################################

# subset.summary.scantwoperm ---------------------------------------------------
#' Subset a \code{summary.scantwoperm} object.
#'
#' @param x A \code{summary.scantwoperm} object.
#' @param alphas Significance levels by which to subset the input object.
#' @template param-lodcolumns
#' @param ... Unused arguments.
#' 
#' @return The input \code{summary.scantwoperm} object,
#' subsetted with respect to the specified constraints.
#'  
#' @export subset.summary.scantwoperm
#' @keywords internal
#' @method subset summary.scantwoperm
#' @rdname subset.summary.scantwoperm
subset.summary.scantwoperm <- function(x, alphas=NULL, lodcolumns=NULL, ...) {
    
    stopifnot( identical(names(x), const$scantwo.lodtypes$scantwoperm) )
    
    lodcol.mask <- NULL
    alpha.mask <- NULL
    
    for ( i in seq_along(x) ) {
        
        mask <- getColMask(x[[i]], requested=lodcolumns)
        if ( ! is.null(lodcol.mask) ) {
            stopifnot( identical(mask, lodcol.mask) )
        } else {
            lodcol.mask <- mask
        }
        
        mask <- getRowMask(x[[i]], requested=alphas)
        if ( ! is.null(alpha.mask) ) {
            stopifnot( identical(mask, alpha.mask) )
        } else {
            alpha.mask <- mask
        }
    }
    
    for ( i in seq_along(x) ) {
        x[[i]] <- x[[i]][alpha.mask, lodcol.mask, drop=FALSE]
    }
    
    return(x)
}

# End of subset.summary.scantwoperm.R ##########################################