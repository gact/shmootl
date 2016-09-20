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
#' interval that defines the given bin (e.g. \code{'LOD[1.1,1.2)'} for a bin that
#' contains the number of loci with LOD scores greater than or equal to 1.1 and
#' less than 1.2).
#' 
#' @export
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
        
        lods <- lods[ ! is.na(lods) & lods >= const$lod.bin$min.start ]
        
        if ( length(lods) > 0 ) {
            
            # Floor LOD values to nearest bin start.
            floored.lods <- lods - (lods %% const$lod.bin$size)
            
            lod.bin.max.start <- max(floored.lods)
            
            bin.starts <- seq(from=const$lod.bin$min.start,
                to=lod.bin.max.start, by=const$lod.bin$size)
            
            bin.labels <- makeLODBinLabels(bin.starts)
            
            binned.lods <- rep_len(0, length(bin.labels))
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
    
    class(lod.bins) <- c('scanonebins', 'array')
    
    return(lod.bins)
}

# End of scanonebins.R #########################################################