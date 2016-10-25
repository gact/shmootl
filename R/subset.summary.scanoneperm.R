# Start of subset.summary.scanoneperm.R ########################################

# subset.summary.scanoneperm ---------------------------------------------------
#' Subset a \code{summary.scanoneperm} object.
#'
#' @param x A \code{summary.scanoneperm} object.
#' @param alphas Significance levels by which to subset the input object.
#' @template param-lodcolumns
#' @param ... Unused arguments.
#' 
#' @return The input \code{summary.scanoneperm} object,
#' subsetted with respect to the specified constraints.
#'  
#' @export subset.summary.scanoneperm
#' @keywords internal
#' @method subset summary.scanoneperm
#' @rdname subset.summary.scanoneperm
subset.summary.scanoneperm <- function(x, alphas=NULL, lodcolumns=NULL, ...) {
    
    alpha.mask <- getRowMask(x, requested=alphas)
    lodcol.mask <- getColMask(x, requested=lodcolumns)
    
    others <- otherattributes(x)
    
    x <- x[alpha.mask, lodcol.mask, drop=FALSE]
    
    if ( ncol(x) == 1 ) {
        colnames(x) <- 'lod'
    }
    
    class(x) <- 'summary.scanoneperm'
    
    otherattributes(x) <- others
    
    return(x)
}

# End of subset.summary.scanoneperm.R ##########################################