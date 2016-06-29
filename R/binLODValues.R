# Start of scanonebins.R #######################################################

# binLODValues -----------------------------------------------------------------
#' Bin LOD values of \code{scanone} object.
#'
#' @param x A \code{scanone} object or equivalent \code{data.frame} in which the
#' data columns contain LOD score values.
#' @template param-lodcolumns
#' 
#' @return A \code{scanonebins} array of LOD bin counts, with a single row, each
#' column corresponding to a bin spanning an interval of LOD values, and each
#' slice corresponding to a LOD column. Each column name indicates the half-open
#' interval that defines the given bin (e.g. \code{'LOD[1,2)'} for a bin that
#' contains the number of loci with LOD scores greater than or equal to one and
#' less than two). Returned bins cover the range of LOD values from 1.0 until the
#' lowest integer that is greater than the maximum observed LOD value.
#' 
#' @export
#' @include scanonebins.R
#' @rdname binLODValues
binLODValues <- function(x, lodcolumns=NULL) {
    
    stopifnot( 'scanone' %in% class(x) )
    
    lodcol.indices <- getDatColIndices(x, datcolumns=lodcolumns, strict=TRUE)
    
    num.lodcols <- length(lodcol.indices)
    
    stopifnot( num.lodcols > 0 )
    
    lodcol.names <- colnames(x)[lodcol.indices]
    
    lod.bin.list <- vector('list', num.lodcols)
    
    for ( i in 1:num.lodcols ) {
        
        lods <- x[, lodcol.indices[i]]
        
        lods <- lods[ ! is.na(lods) & lods >= 1.0 ]
        
        if ( length(lods) > 0 ) {
            
            floored.lods <- floor(lods)
            
            num.bins <- max(floored.lods)
            
            bin.starts <- seq_len(num.bins)
            
            bin.labels <- makeLODBinLabels(bin.starts)
            
            binned.lods <- rep_len(0, num.bins)
            names(binned.lods) <- bin.labels
            
            bin.freqs <- table(floored.lods)
            
            binned.lods[ bin.starts %in% names(bin.freqs) ] <- bin.freqs
            
        } else {
            
            binned.lods <- integer()
        }
        
        lod.bin.list[[i]] <- binned.lods
    }
    
    max.idx <- which.max( lengths(lod.bin.list) )
    
    bin.labels <- names(lod.bin.list[[max.idx]])
    
    num.bins <- length(bin.labels)
    
    for ( i in getIndices(lod.bin.list) ) {
        length(lod.bin.list[[i]]) <- num.bins
    }
    
    lod.bins <- array( unlist(lod.bin.list), dim=c(1, num.bins, num.lodcols),
        dimnames=list(NULL, bin.labels, lodcol.names) )
    
    lod.bins[ is.na(lod.bins) ] <- 0
    
    attr(lod.bins, 'n.loci') <- nrow(x)
    
    class(lod.bins) <- c('scanonebins', 'array')
    
    return(lod.bins)
}

# End of scanonebins.R #########################################################