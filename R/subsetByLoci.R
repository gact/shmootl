# Start of subsetByLoci.R ######################################################

# subsetByLoci -----------------------------------------------------------------
#' Subset object by the specified loci.
#' 
#' @param x Object subsettable by loci.
#' @param loc Locus \code{mapframe} specifying map loci.
#'     
#' @return Input object subsetted by the specified loci.
#' 
#' @keywords internal
#' @rdname subsetByLoci
subsetByLoci <- function(x, loc=NULL) {
    UseMethod('subsetByLoci', x)
}

# subsetByLoci.map -------------------------------------------------------------
#' @rdname subsetByLoci
subsetByLoci.map <- function(x, loc=NULL) {
    
    if ( ! is.null(loc) ) {
    
        stopifnot( 'mapframe' %in% class(loc) )
        
        # Get map unit.
        map.unit <- getMapUnit(x)
        
        if ( is.na(map.unit) ) {
            stop("map unit not found")
        }
        
        if ( getMapUnit(loc) != map.unit ) {
            stop("map unit mismatch")
        }
        
        if ( nrow(loc) == 0 ) {
            stop("locus mapframe must have at least one row")
        }
        
        # Get sequences and positions of specified loci.
        loc.seqs <- pullLocusSeq(loc)
        loc.pos <- pullLocusPos(loc)
        
        # Get map sequences.
        map.seqs <- names(x)
        
        # Get normalised map sequences.
        norm.map.seqs <- normSeq(map.seqs)
        
        stopifnot( all( loc.seqs %in% norm.map.seqs ) )
        
        # Get vector mapping original sequence labels to their normalised form.
        map2res <- structure(norm.map.seqs, names=map.seqs)
        
        # Get subset sequences.
        subset.seqs <- unique(loc.seqs)
        
        # Get pseudo-map for the specified loci.
        pseudo.map <- lapply(subset.seqs, function(subset.seq) 
            loc.pos[ loc.seqs == subset.seq ])
        names(pseudo.map) <- subset.seqs
        
        # Subset the map by matching sequences.
        x <- subsetMap(x, norm.map.seqs %in% subset.seqs)
        
        # Subset the map by matching locus positions.
        for ( map.seq in map.seqs ) {
            subset.pos <- pseudo.map[[ map2res[map.seq] ]]
            subset.mask <- sapply(x[[map.seq]], function (pos) 
                min( abs(pos - subset.pos) ) <= .Machine$double.eps^0.5 )
            x[[map.seq]] <- x[[map.seq]][subset.mask]
        }
        
        # Keep only sequences with enough markers.
        x <- subsetMap(x, lengths(x) >= const$min.lps)
        
        if ( length(x) < const$min.spm ) {
            stop("subsetted map has too few sequences (min=", const$min.spm, ")")
        }
    }
    
    return(x)
}

# subsetByLoci.data.frame ------------------------------------------------------
#' @rdname subsetByLoci
subsetByLoci.mapframe <- function(x, loc=NULL) {
    
    if ( ! is.null(loc) ) {
        indices <- matchLocusRowIndices(x, loc)
        x <- x[indices, ]
    }
    
    return(x)
}

# End of subsetByLoci.R ########################################################