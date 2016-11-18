# Start of getLODMatrix.R ######################################################

# getLODMatrix -----------------------------------------------------------------
#' Get LOD matrix.
#' 
#' Given a \code{scantwo} object, this function extracts a LOD matrix, in which
#' the upper and lower triangle contain LOD values of a type specified by the
#' parameters \code{upper} and \code{lower}, respectively. The method used to
#' obtain derived LOD values is that of \pkg{R/qtl} function \code{plot.scantwo}
#' (Broman \emph{et al.} 2003).
#' 
#' @param x An \pkg{R/qtl} \code{scantwo} object.
#' @param chr Vector indicating which sequences to include in the LOD matrix.
#' @param incl.markers When \code{scantwo} has been run with a fixed non-zero
#' step width, this indicates if markers (as distinct from pseudomarkers) should
#' be included in the returned LOD matrix.
#' @template param-lodcolumn
#' @param lower Type of LOD score to put in the lower triangle of the LOD matrix.
#' @param upper Type of LOD score to put in the upper triangle of the LOD matrix.
#' 
#' @return A LOD matrix taken from the input \code{scantwo} object, containing
#' LOD scores for a single phenotype specified by the \code{lodcolumn} parameter.
#' 
#' @template ref-broman-2003
#' 
#' @seealso Function \code{plot.scantwo} in the
#' \href{http://www.rqtl.org/manual/qtl-manual.pdf}{R/qtl manual}.
#' 
#' @export
#' @rdname getLODMatrix
getLODMatrix <- function(x, chr=NULL, incl.markers=FALSE, lodcolumn=NULL,
    lower=c('full', 'cond-int', 'int', 'cond-add', 'add'),
    upper=c('add', 'cond-add', 'int', 'cond-int', 'full') ) {
    
    stopifnot( 'scantwo' %in% class(x) )
    
    if ( isSingleString(lower) ) {
        lower <- resolveScantwoLodtypes(lower, to='plot.scantwo')
    }
    lower <- match.arg(lower)
    
    if ( isSingleString(upper) ) {
        upper <- resolveScantwoLodtypes(upper, to='plot.scantwo')
    }
    upper <- match.arg(upper)
    
    # Get specified sequences.
    x.seqs <- pull.chr(x$map)
    chr <- x.seqs <- subsetBySeq(x.seqs, chr)
    stopifnot( length(chr) > 0 )
    
    # Subset scantwo object by specified sequences.
    if ( ! is.null(chr) ) {
        x <- subsetBySeq(x, chr)
    }
    
    # Get LOD column index.
    lodcol.index <- getLodColIndex(x, lodcolumn=lodcolumn)
    
    # Get scantwo LOD matrix for given LOD column.
    if ( ! is.matrix(x$lod) ) {
        scantwo.matrix <- x$lod[,, lodcol.index]
    } else {
        scantwo.matrix <- x$lod
    }
    
    # Remove markers if requested.
    if ( ! incl.markers && any( x$map[, 'eq.spacing'] == 0 ) ) {
        indices  <- which( x$map[, 'eq.spacing'] == 1 )
        scantwo.matrix <- scantwo.matrix[indices, indices]
    }
    
    # Get scantwo map per-locus sequences and LOD values.
    locus.seqs <- pullLocusSeq(x$map)
    locus.lods <- diag(scantwo.matrix)
    
    if ( all( locus.lods < .Machine$double.eps ) ) {
        lodtypes <- c(lower, upper)
        incalculable <- lodtypes[ lodtypes %in% c('cond-int', 'cond-add') ]
        if ( length(incalculable) > 0 ) {
            stop("without scanone results, cannot get LOD matrix for LOD types - '",
                toString(incalculable), "'")
        }
    }
    
    L <- lower.tri(scantwo.matrix)
    U <- upper.tri(scantwo.matrix)
    
    # Init LOD matrix to input scantwo matrix.
    lod.matrix <- scantwo.matrix
    
    # If requested LOD type for lower triangle is not the same
    # as the input scantwo matrix, calculate derived LOD values.
    if ( lower != 'full' ) {
        
        if ( lower == 'add' ) {
            
            lod.matrix[L] <- t(scantwo.matrix)[L]
            
        } else if ( lower == 'int' ) {
            
            lod.matrix[L] <- scantwo.matrix[L] - t(scantwo.matrix)[L]
            
        } else if ( lower %in% c('cond-int', 'cond-add') ) {
            
            if ( lower == 'cond-add' ) {
                lod.matrix[L] <- t(scantwo.matrix)[L]
            }
            
            max.lod.one <- tapply(locus.lods, locus.seqs, max, na.rm=TRUE)
            
            for ( i in seq_along(x.seqs) ) {
                
                cols <- matchSeqRowIndices(x$map, x.seqs[i], simplify=TRUE)
                
                for ( j in seq( i, length(x.seqs) ) ) {
                    
                    rows <- matchSeqRowIndices(x$map, x.seqs[j], simplify=TRUE)
                    
                    seq.matrix <- lod.matrix[rows, cols] - max(max.lod.one[c(i,j)])
                    seq.matrix[ isNegativeNumber(seq.matrix) ] <- 0
                    
                    if ( i == j ) {
                        SL <- lower.tri(seq.matrix)
                        lod.matrix[rows, cols][SL] <- seq.matrix[SL]
                    } else {
                        lod.matrix[rows, cols] <- seq.matrix
                    }
                }
            }
        }
    }
    
    # If requested LOD type for upper triangle is not the same
    # as the input scantwo matrix, calculate derived LOD values.
    if ( upper != 'add' ) {
        
        if ( upper == 'full' ) {
            
            lod.matrix[U] <- t(scantwo.matrix)[U]
            
        } else if ( upper == 'int' ) {
            
            lod.matrix[U] <- t(scantwo.matrix)[U] - scantwo.matrix[U]
            
        } else if ( upper %in% c('cond-int', 'cond-add') ) {
            
            if ( upper == 'cond-int' ) {
                lod.matrix[U] <- t(scantwo.matrix)[U]
            }
            
            max.lod.one <- tapply(locus.lods, locus.seqs, max, na.rm=TRUE)
            
            for ( i in seq_along(x.seqs) ) {
                
                rows <- matchSeqRowIndices(x$map, x.seqs[i], simplify=TRUE)
                
                for ( j in seq( i, length(x.seqs) ) ) {
                    
                    cols <- matchSeqRowIndices(x$map, x.seqs[j], simplify=TRUE)
                    
                    seq.matrix <- lod.matrix[rows, cols] - max(max.lod.one[c(i,j)])
                    seq.matrix[ isNegativeNumber(seq.matrix) ] <- 0
                    
                    if ( i == j ) {
                        SU <- upper.tri(seq.matrix)
                        lod.matrix[rows, cols][SU] <- seq.matrix[SU]
                    } else {
                        lod.matrix[rows, cols] <- seq.matrix
                    }
                }
            }
        }
    }
    
    return(lod.matrix)
}

# End of getLODMatrix.R ########################################################