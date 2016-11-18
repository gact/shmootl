# Start of getLODProfile.R #####################################################

# getLODProfile ----------------------------------------------------------------
#' Get LOD profile.
#' 
#' @usage 
#' ## Generic method.
#' getLODProfile(x, ...)
#'
#' @param x Object containing single-QTL LOD profile data.
#' @param ... Further arguments (see below).
#' @template param-lodcolumn
#' @param qtl.index In a \code{qtl} object, this indicates
#' the QTL for which a LOD profile should be returned.
#' 
#' @return A \code{scanone} object containing
#' a LOD profile for a single phenotype.
#' 
#' @export
#' @rdname getLODProfile
getLODProfile <- function(x, ...) {
    UseMethod('getLODProfile', x)
}

# getLODProfile.data.frame -----------------------------------------------------
#' @export
#' @rdname getLODProfile
getLODProfile.data.frame <- function(x, lodcolumn=NULL, ...) {
    
    stopifnot( getMapUnit(x) == 'cM' )
    stopifnot( nrow(x) > 0 )
    
    seqcol.index <- getSeqColIndex(x)
    poscol.index <- getPosColIndex(x)
    lodcol.index <- getLodColIndex(x, lodcolumn=lodcolumn)
    
    others <- otherattributes(x)
    
    # Create LOD profile with specified LOD column.
    lod.profile <- x[, c(seqcol.index, poscol.index, lodcol.index)]
    colnames(lod.profile) <- c('chr', 'pos', 'lod')
    class(lod.profile) <- c('scanone', 'data.frame')
    
    otherattributes(lod.profile) <- others
    
    return(lod.profile)
}

# getLODProfile.qtl ------------------------------------------------------------
#' @export
#' @rdname getLODProfile
getLODProfile.qtl <- function(x, qtl.index=NULL, ...) {

    stopifnot( 'lodprofile' %in% names( attributes(x) ) )
    
    qtl.index <- getIndices(x, requested=qtl.index)
    
    if ( length(qtl.index) > 1 ) {
        stop("cannot get multiple LOD profiles from QTL object - please choose one")
    } else if ( length(qtl.index) == 0 ) {
        stop("no QTL index specified")
    }
        
    # Get LOD profiles of QTL object: a 
    # per-QTL list of scanone objects.
    lod.profiles <- attr(x, 'lodprofile')    
        
    # Get LOD profile for given QTL.
    lod.profile <- lod.profiles[[qtl.index]]
    
    return(lod.profile)
}

# getLODProfile.scantwo --------------------------------------------------------
#' @export
#' @rdname getLODProfile
getLODProfile.scantwo <- function(x, lodcolumn=NULL, ...) {
    
    lodcol.index <- getLodColIndex(x, lodcolumn=lodcolumn)
    
    if ( ! is.matrix(x$lod) ) {
        scantwo.matrix <- x$lod[,, lodcol.index]
    } else {
        scantwo.matrix <- x$lod
    }
    
    locus.ids <- pullLocusIDs(x$map)
    locus.seqs <- pullLocusSeq(x$map)
    locus.pos <- pullLocusPos(x$map)
    locus.lods <- diag(scantwo.matrix)
    
    lod.profile <- data.frame(
        row.names=locus.ids,
        chr=locus.seqs,
        pos=locus.pos,
        lod=locus.lods,
    stringsAsFactors = default.stringsAsFactors() )
    class(lod.profile) <- c('scanone', 'data.frame')
    
    return(lod.profile)
}

# End of getLODProfile.R #######################################################