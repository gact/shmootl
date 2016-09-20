# Start of getLODValues.R ######################################################

# getLODValues -----------------------------------------------------------------
#' Get LOD values.
#' 
#' Get the LOD value at each of the specified
#' map positions in the given LOD profile.
#'
#' @param x A \code{scanone} object or equivalent \code{data.frame}
#' in which the data columns contain LOD score values.
#' @param loc Locus \code{mapframe} specifying map positions.
#' @template param-lodcolumn
#' 
#' @return Vector of LOD values. Returns \code{NA} if
#' LOD score is not available at a given position.
#' 
#' @export
#' @importFrom stats approx
#' @rdname getLODValues
getLODValues <- function(x, loc, lodcolumn=NULL) {
    
    stopifnot( 'mapframe' %in% class(x) || 'scanone' %in% class(x) )
    stopifnot( getMapUnit(x) == 'cM' )
    stopifnot( nrow(x) > 0 )
    stopifnot( 'mapframe' %in% class(loc) )
    stopifnot( getMapUnit(loc) == 'cM' )
    
    # Get pos column index.
    poscol.index <- getPosColIndex(x)

    # Get LOD column index.
    lodcol.index <- getDatColIndices(x, datcolumns=lodcolumn)
    
    if ( length(lodcol.index) > 1 ) {
        stop("cannot get LOD values for multiple LOD columns - please choose one")
    } else if ( length(lodcol.index) == 0 ) {
        stop("no LOD column found")
    }
    
    if ( nrow(loc) == 0 ) {
        stop("cannot get LOD values - no loci specified")
    }
    
    lod.values <- numeric( nrow(loc) )
    
    # Find LOD value at each specified map position.
    for ( i in getRowIndices(loc) ) {
        
        # Get row indices of flanking loci.
        row.indices <- findFlankingRowIndices(x, loc[i, ])
        
        # Get row index of the closest flanking locus, 
        # taking the first locus in case of ties.
        row.index <- row.indices[ which.min( abs(loc$pos[i] - x[row.indices, poscol.index]) ) ]
        
        # Get distance to closest locus.
        proximity <- abs(loc$pos[i] - x[row.index, poscol.index])
        
        # If closest locus is within numeric tolerance, take its LOD value..
        if ( proximity < .Machine$double.eps^0.5 ) {
            
            lod.values[i] <- x[row.index, lodcol.index]
            
        } else { # ..otherwise interpolate LOD value from closest flanking loci.
            
            x <- x[row.indices, poscol.index]
            y <- x[row.indices, lodcol.index]
            
            if ( anyNA(y) ) {
                lod.values[i] <- NA
            } else {
                interpolation <- stats::approx(x, y, loc$pos, 'linear')
                lod.values[i] <- interpolation$y
            }
        }
    }
    
    return(lod.values)
}

# End of getLODValues.R ########################################################