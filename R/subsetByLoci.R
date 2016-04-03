# Start of subsetByLoci.R ######################################################

# subsetByLoci (S3) ------------------------------------------------------------
#' Subset object by the specified loci.
#' 
#' @param x Object subsettable by loci.
#' @param loc Locus \code{mapframe} specifying map loci.
#'     
#' @return Input object subsetted by the specified loci.
#' 
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings QualityScaledDNAStringSet
#' @keywords internal
#' @rdname subsetByLoci
subsetByLoci <- function(x, loc=NULL) {
    UseMethod('subsetByLoci', x)
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

# subsetByLoci.DNAStringSet ----------------------------------------------------
#' @rdname subsetByLoci
subsetByLoci.DNAStringSet <- function(x, loc=NULL) {
    
    stopifnot( 'loci' %in% names(x@metadata) )
    stopifnot( 'mapframe' %in% class(x@metadata[['loci']]) )
    
    if ( ! is.null(loc) ) {
    
        others <- otherattributes(x)
        meta <- x@metadata
        elem.meta <-x@elementMetadata
        
        indices <- matchLocusRowIndices(meta[['loci']], loc)

        meta[['loci']] <- meta[['loci']][indices, ]
        
        bases <- lapply(Biostrings::as.list(x), function(b) b[indices])
        
        x <- Biostrings::DNAStringSet(bases)
    
        x@elementMetadata <- elem.meta
        x@metadata <- meta
        otherattributes(x) <- others    
    }
    
    return(x)
}

# subsetByLoci.QualityScaledDNAStringSet ---------------------------------------
#' @rdname subsetByLoci
subsetByLoci.QualityScaledDNAStringSet <- function(x, loc=NULL) {

    stopifnot( 'loci' %in% names(x@metadata) )
    stopifnot( 'mapframe' %in% class(x@metadata[['loci']]) )
    
    if ( ! is.null(loc) ) {

        others <- otherattributes(x)
        meta <- x@metadata
        elem.meta <-x@elementMetadata

        indices <- matchLocusRowIndices(meta[['loci']], loc)

        meta[['loci']] <- meta[['loci']][indices, ]

        bases <- lapply(Biostrings::as.list(x), function(b) b[indices])
        quals <- lapply(Biostrings::as.list(x@quality), function(q) q[indices])

        x <- Biostrings::QualityScaledDNAStringSet( 
            Biostrings::DNAStringSet(bases), 
            Biostrings::PhredQuality( Biostrings::BStringSet(quals) )
        )

        x@elementMetadata <- elem.meta
        x@metadata <- meta
        otherattributes(x) <- others
    }
    
    return(x)
}

# subsetByLoci (S4) ------------------------------------------------------------
#' @rdname subsetByLoci
setGeneric('subsetByLoci', subsetByLoci)

# DNAStringSet::subsetByLoci ---------------------------------------------------
#' @rdname subsetByLoci
setMethod('subsetByLoci', signature='DNAStringSet', 
    definition = subsetByLoci.DNAStringSet)

# QualityScaledDNAStringSet::subsetByLoci --------------------------------------
#' @rdname subsetByLoci
setMethod('subsetByLoci', signature='QualityScaledDNAStringSet', 
    definition = subsetByLoci.QualityScaledDNAStringSet)

# End of subsetByLoci.R ########################################################