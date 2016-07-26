# Start of map.R ###############################################################

# Map Utilities ----------------------------------------------------------------
#' Genetic and physical map utilities.
#' 
#' These utilities include functions for \pkg{R/qtl} \code{map} objects,
#' \pkg{shmootl} \code{mapframe} objects. A \code{mapframe} can be created
#' in a similar manner to a \code{data.frame} with \code{\link{mapframe}}.
#' All arguments are passed to the \code{data.frame} constructor, except for
#' an addition \code{map.unit} parameter, which can be used to set the map unit
#' of the new \code{mapframe}. The function \code{\link{gmapframe}} creates a
#' genetic \code{mapframe} with map unit set automatically to \code{'cM'}.
#' Coercion functions exist to convert a \code{map} to a \code{mapframe}
#' (\code{\link{as.mapframe}}), and vice versa (\code{\link{as.map}}).
#' 
#' @template map
#' @template mapframe
#' @template mapunits
#' 
#' @docType package
#' @name Map Utilities
NULL

# as.data.frame.map ------------------------------------------------------------
#' Coerce \code{map} to \code{data.frame}.
#' 
#' @param x An \pkg{R/qtl} \code{map} object.
#' @param ... Unused arguments.
#' @param map.unit Explicitly sets \code{'map.unit'}.
#' 
#' @return A \code{data.frame} corresponding to the input \code{map} object.
#' 
#' @export
#' @family map utilities
#' @method as.data.frame map
#' @rdname as.data.frame.map
as.data.frame.map <- function(x, ..., map.unit=NULL) {

    # If map unit specified, update that of object..
    if ( ! is.null(map.unit) ) {
        x <- setMapUnit(x, map.unit)
    } else { # ..otherwise get map unit from object.
        map.unit <- getMapUnit(x)
    }
    
    validateMap(x)
    
    # Get sequences corresponding to individual map loci.
    x.seqs <- pullLocusSeq(x)
    
    # Get individual locus positions.
    x.pos <- pullLocusPos(x)
    
    # Get individual locus IDs.
    x.ids <- pullLocusIDs(x)
    
    to <- data.frame(chr=x.seqs, pos=x.pos, row.names=x.ids)
    
    to <- orderMap(to)
    
    if ( ! is.na(map.unit) ) {
        to <- setMapUnit(to, map.unit)
    }
    
    return(to)
}

# as.data.frame.mapframe -------------------------------------------------------
#' Coerce \code{mapframe} to \code{data.frame}.
#' 
#' @param x A \code{mapframe} object.
#' @param ... Unused arguments.
#'     
#' @return A \code{data.frame} corresponding to the input \code{mapframe} object.
#' 
#' @export
#' @family map utilities
#' @method as.data.frame mapframe
#' @rdname as.data.frame.mapframe
as.data.frame.mapframe <- function(x, ...) {
    class(x) <- 'data.frame'
    return(x)
}

# as.map -----------------------------------------------------------------------
#' Coerce object to \code{map}.
#' 
#' @param from Object containing map data.
#' @param map.unit Explicitly sets \code{'map.unit'}. This must be specified if
#' an object does not have an existing \code{'map.unit'} attribute.
#' 
#' @return An \pkg{R/qtl} \code{map} corresponding to the input object.
#' 
#' @export
#' @family map utilities
#' @rdname as.map
as.map <- function(from, map.unit=NULL) {
    UseMethod('as.map', from)
}

# as.map.data.frame ------------------------------------------------------------
#' @export
#' @method as.map data.frame
#' @rdname as.map
as.map.data.frame <- function(from, map.unit=NULL) {

    if ( nrow(from) == 0 ) {
        stop("cannot create map from empty data frame")
    }
    
    # If map unit specified, update that of object..
    if ( ! is.null(map.unit) ) {
        from <- setMapUnit(from, map.unit)
    } else { # ..otherwise get map unit from object.
        map.unit <- getMapUnit(from)
    }
    
    validateMapframe(from)
    
    # Get normalised sequences corresponding to individual map loci.
    norm.from.seqs <- normSeq( pullLocusSeq(from) )
    
    # Get individual locus positions.
    from.pos <- pullLocusPos(from) # NB: strips any map-unit suffixes.
    
    # Get individual locus IDs.
    from.ids <- pullLocusIDs(from)
    
    if ( is.null(from.ids) ) {
        stop("cannot create map from data frame - no locus IDs found")
    }
    
    # Get normalised map sequences.
    norm.map.seqs <- sortSeq( unique(norm.from.seqs) )
    
    if ( length(norm.map.seqs) < const$min.spm  ) {
        stop("cannot create map from data frame - too few sequences (min=",
            const$min.spm, ")")
    }
    
    to <- list()
    
    # Set map data for each sequence.
    for ( norm.map.seq in norm.map.seqs ) {
        
        # Get row indices of data frame for this sequence.
        indices <- which( norm.from.seqs == norm.map.seq )
        
        # Get map positions for this sequence.
        seq.pos <- from.pos[indices]
        
        if ( length(seq.pos) < const$min.lps ) {
            stop("cannot create map from data frame - sequence has too few loci - '", seq, "'")
        }
        
        # Set map position names from locus IDs.
        names(seq.pos) <- from.ids[indices]
        
        # Sort map positions in numeric order.
        seq.pos <- sort.int(seq.pos)
        
        # Set class of map position vector.
        class(seq.pos) <- 'A' # NB: assumes no 'X' chromosomes.
        
        # Set map positions for sequence.
        to[[norm.map.seq]] <- seq.pos
    }
    
    to <- setMapUnit(to, map.unit)
    
    class(to) <- 'map'
    
    return(to)
}

# as.map.list ------------------------------------------------------------------
#' @export
#' @method as.map list
#' @rdname as.map
as.map.list <- function(from, map.unit=NULL) {

    # If map unit specified, update that of object..
    if ( ! is.null(map.unit) ) {
        from <- setMapUnit(from, map.unit)
    } else { # ..otherwise get map unit from object.
        map.unit <- getMapUnit(from)
    }
    
    validateMap(from)
    
    # Get map sequences.
    map.seqs <- names(from)
    
    # Get normalised map sequences.
    norm.map.seqs <- normSeq(map.seqs)
    
    # Get vector mapping normalised sequence labels to their original form.
    res2map <- structure(map.seqs, names=norm.map.seqs)

    to <- list()
    
    # Set map data for each sequence.
    for ( norm.map.seq in sortSeq(norm.map.seqs) ) {
        
        # Get map positions for this sequence.
        seq.pos <- from[[ res2map[norm.map.seq] ]]
        
        # Sort map positions in numeric order.
        seq.pos <- sort.int(seq.pos)
        
        # Set class of map position vector.
        class(seq.pos) <- 'A' # NB: assumes no 'X' chromosomes.
        
        # Set map positions for sequence.
        to[[norm.map.seq]] <- seq.pos
    }
    
    to <- setMapUnit(to, map.unit)
    
    class(to) <- 'map'
    
    return(to)
}


# as.map.map -------------------------------------------------------------------
#' @export
#' @method as.map map
#' @rdname as.map
as.map.map <- function(from, map.unit=NULL) {
    
    if ( ! is.null(map.unit) ) {
        from <- setMapUnit(from, map.unit)
    }
    
    validateMap(from)
    
    return(from)
}

# as.map.mapframe --------------------------------------------------------------
#' @export
#' @method as.map mapframe
#' @rdname as.map
as.map.mapframe <- function(from, map.unit=NULL) {
    return( as.map.data.frame(from, map.unit=map.unit) )
}

# as.map.scanone ---------------------------------------------------------------
#' @export
#' @method as.map scanone
#' @rdname as.map
as.map.scanone <- function(from, map.unit=NULL) {
    stopifnot( map.unit == 'cM' )
    return( as.map.data.frame(from, map.unit='cM') )
}

# as.mapframe ------------------------------------------------------------------
#' Coerce object to \code{mapframe}.
#' 
#' @param from Object containing map data.
#' @param map.unit Explicitly sets \code{'map.unit'}. This must be specified if
#' an object does not have an existing \code{'map.unit'} attribute.
#' 
#' @return A \code{mapframe} corresponding to the input object.
#' 
#' @export
#' @family map utilities
#' @rdname as.mapframe
as.mapframe <- function(from, map.unit=NULL) {
    UseMethod('as.mapframe', from)
}

# as.mapframe.data.frame -------------------------------------------------------
#' @export
#' @method as.mapframe data.frame
#' @rdname as.mapframe
as.mapframe.data.frame <- function(from, map.unit=NULL) {
    
    # If map unit specified, update that of object..
    if ( ! is.null(map.unit) ) {
        from <- setMapUnit(from, map.unit)
    } else { # ..otherwise get map unit from object.
        map.unit <- getMapUnit(from)
    }
    
    validateMapframe(from)
    
    seqcol.index <- getSeqColIndex(from)
    poscol.index <- getPosColIndex(from)
    
    # Ensure data frame rows are valid, if present.
    if ( nrow(from) > 0 ) {
        
        # Get normalised sequences corresponding to individual map loci.
        norm.from.seqs <- normSeq( pullLocusSeq(from) )
        
        # Get individual locus positions.
        from.pos <- pullLocusPos(from) # NB: strips any map-unit suffixes.

        from[[seqcol.index]] <- norm.from.seqs
        from[[poscol.index]] <- from.pos
        
        # Order data frame by sequence, then map position.
        from <- from[ order(rankSeq(norm.from.seqs), from.pos), ]
    }
    
    # Ensure column headings are as expected for a mapframe.
    colnames(from)[ c(seqcol.index, poscol.index) ] <- c('chr', 'pos')
    
    # Ensure sequence and positions columns are leftmost.
    if ( seqcol.index != 1 || poscol.index != 2 ) {
        column.indices <- getColIndices(from)
        other.indices <- column.indices[ -c(seqcol.index, poscol.index) ]
        from[, column.indices] <- from[, c( seqcol.index, poscol.index, other.indices)]
    }
    
    from <- setMapUnit(from, map.unit)
    
    class(from) <- c('mapframe', 'data.frame')
    
    return(from)
}

# as.mapframe.list -------------------------------------------------------------
#' @export
#' @method as.mapframe list
#' @rdname as.mapframe
as.mapframe.list <- function(from, map.unit=NULL) {
    return( as.mapframe( as.map(from, map.unit=map.unit) ) )
}

# as.mapframe.map --------------------------------------------------------------
#' @export
#' @method as.mapframe map
#' @rdname as.mapframe
as.mapframe.map <- function(from, map.unit=NULL) {
    return( as.mapframe( as.data.frame(from), map.unit=map.unit ) )
}

# as.mapframe.mapframe ---------------------------------------------------------
#' @export
#' @method as.mapframe mapframe
#' @rdname as.mapframe
as.mapframe.mapframe <- function(from, map.unit=NULL) {
    
    if ( ! is.null(map.unit) ) {
        from <- setMapUnit(from, map.unit)
    }
    
    validateMapframe(from)
    
    return(from)
}

# as.mapframe.scanone ----------------------------------------------------------
#' @export
#' @method as.mapframe scanone
#' @rdname as.mapframe
as.mapframe.scanone <- function(from, map.unit=NULL) {
    stopifnot( map.unit == 'cM' )
    return( as.mapframe.data.frame(from, map.unit='cM') )
}

# convertMapUnit ---------------------------------------------------------------
#' Convert map unit and rescale map positions.
#' 
#' @template mapunits
#' 
#' @param x Object containing map data.
#' @param map.unit New map unit.
#' 
#' @return Input object with converted map unit and rescaled map positions.
#' 
#' @keywords internal
#' @rdname convertMapUnit
convertMapUnit <- function(x, map.unit) {
    
    # If object is a data-frame, check that map units have
    # been removed from position column names and data.
    if ( is.data.frame(x) ) {
        
        if ( ! is.na( getPosColNameMapUnit(x) ) ) {
            stop("function 'convertMapUnit' cannot take object with map unit ",
                "in column name - consider calling 'setMapUnit' instead")
        }
        
        if ( ! is.na( getPosColDataMapUnit(x) ) ) {
            stop("function 'convertMapUnit' cannot take object with map unit ",
                "in map positions - consider calling 'setMapUnit' instead")
        }
    }
    
    # Get current mapkey.
    map.key <- mapkeyOpt()
    
    # Pull locus mapframe from object.
    loc <- pullLoci(x)
    
    # Apply mapkey to loci.
    loc <- map.key(loc, map.unit)
    
    # Set map unit directly.
    # NB: map unit of object must match the locus mapframe, but we
    # cannot set this through setMapUnit as it calls this function.
    attr(x, 'map.unit') <- map.unit
    
    # Replace loci of object.
    x <- pushLoci(x, loc)
    
    return(x)
}

# extractMarkers ---------------------------------------------------------------
#' Extract marker loci.
#' 
#' @param x Object containing map data.
#' 
#' @return Input object with pseudomarkers removed.
#' 
#' @export
#' @family map utilities
#' @rdname extractMarkers
extractMarkers <- function(x) {
    return( subsetByLocusID(x, isMarkerID) )
}

# extractPseudomarkers ---------------------------------------------------------
#' Extract pseudomarker loci.
#' 
#' @param x Object containing map data.
#' 
#' @return Input object with markers removed.
#' 
#' @export
#' @family map utilities
#' @rdname extractPseudomarkers
extractPseudomarkers <- function(x) {
    return( subsetByLocusID(x, isPseudomarkerID) )
}

# findFlanking -----------------------------------------------------------------
#' Find map loci flanking map positions.
#' 
#' @param x Object containing map data.
#' @param loc Locus \code{mapframe} specifying map positions.
#' @param expandtomarkers Expand the flanking interval to the nearest flanking
#' markers, or to the respective terminal loci.
#' 
#' @return A \code{mapframe} containing two loci marking respectively the
#' lower and upper limit of the flanking interval.
#' 
#' @export
#' @family map utilities
#' @rdname findFlanking
findFlanking <- function(x, loc, expandtomarkers=FALSE) {
    UseMethod('findFlanking', x)
}

# findFlanking.map -------------------------------------------------------------
#' @export
#' @rdname findFlanking
findFlanking.map <- function(x, loc, expandtomarkers=FALSE) {
    
    stopifnot( 'mapframe' %in% class(loc) )
    stopifnot( isBOOL(expandtomarkers) )
    
    # Get map unit.
    map.unit <- getMapUnit(x)
    
    if ( is.na(map.unit) ) {
        stop("map unit not found")
    }
    
    if ( getMapUnit(loc) != map.unit ) {
        stop("map unit mismatch")
    }
    
    # Get locus mapframe sequence.
    loc.seq <- unique( pullLocusSeq(loc) )
    
    # Check that locus mapframe refers to one sequence.
    # NB: implicitly checks that loc has one or more rows.
    stopifnot( length(loc.seq) == 1 )
    
    # Get map sequences.
    map.seqs <- names(x)
    
    # Get normalised map sequences.
    norm.map.seqs <- normSeq(map.seqs)
    
    stopifnot( loc.seq %in% norm.map.seqs )
    
    # Get vector mapping normalised sequence labels to their original form.
    res2map <- structure(map.seqs, names=norm.map.seqs)

    # Get all map positions for this sequence.
    seq.pos <- x[[ res2map[loc.seq] ]]
    
    exrange <- loc$pos[ ! sapply( loc$pos, inRange, range(seq.pos) ) ]
    if ( length(exrange) > 0 ) {
        stop("locus positions out of range - ", toString(exrange), "'")
    }

    if (expandtomarkers) {
        # If expanding to markers, create mask in which
        # markers and sequence endpoints are TRUE..
        loc.mask <- isMarkerID( names(seq.pos) ) |
            getIndices(seq.pos) %in% c(1, length(seq.pos))
    } else {
        # ..otherwise create mask with all TRUE.
        loc.mask <- rep_len( TRUE, length(seq.pos) )
    }
    
    # Init flanking locus info.
    flank.seq <- rep_len(loc.seq, 2)
    flank.pos <- loc$pos[ c(1, nrow(loc)) ]
    flank.ids <- character(2)
    
    # Get lower flanking locus info.
    lower.flanking <- seq.pos[ loc.mask & flank.pos[1] >= seq.pos ]
    lower.diff <- flank.pos[1] - lower.flanking
    lower.indices <- which( lower.diff == min(lower.diff) )
    L <- lower.indices[ length(lower.indices) ]
    flank.pos[1] <- lower.flanking[L]
    flank.ids[1] <- names(lower.flanking)[L]
    
    # Get upper flanking locus info.
    upper.flanking <- seq.pos[ loc.mask & seq.pos >= flank.pos[2] ]
    upper.diff <- upper.flanking - flank.pos[2]
    upper.indices <- which( upper.diff == min(upper.diff) )
    U <- upper.indices[1]
    flank.pos[2] <- upper.flanking[U]
    flank.ids[2] <- names(upper.flanking)[U]
    
    return( mapframe(chr=flank.seq, pos=flank.pos,
        row.names=flank.ids, map.unit=map.unit) )
}

# findFlanking.mapframe --------------------------------------------------------
#' @export
#' @rdname findFlanking
findFlanking.mapframe <- function(x, loc, expandtomarkers=FALSE) {
    row.indices <- findFlankingRowIndices(x, loc, expandtomarkers=expandtomarkers)
    return( mapframe(chr=as.character(x$chr[row.indices]), pos=x$pos[row.indices],
        row.names=rownames(x)[row.indices], map.unit=getMapUnit(x)) )
}

# findFlanking.scanone ---------------------------------------------------------
#' @export
#' @rdname findFlanking
findFlanking.scanone <- function(x, loc, expandtomarkers=FALSE) {
    return( findFlanking.mapframe(x, loc, expandtomarkers=expandtomarkers) )
}

# findFlankingRowIndices -------------------------------------------------------
#' Find row index of flanking \code{mapframe} loci.
#' 
#' @param x A \code{mapframe} or \code{scanone} object.
#' @param loc Locus \code{mapframe} specifying map positions.
#' @param expandtomarkers Expand the flanking interval to the nearest flanking
#' markers, or to the respective terminal loci.
#' 
#' @return Vector containing the row indices of flanking loci.
#' 
#' @keywords internal
#' @rdname findFlankingRowIndices
findFlankingRowIndices <- function(x, loc, expandtomarkers=FALSE) {
    
    stopifnot( 'mapframe' %in% class(x) || 'scanone' %in% class(x) )
    stopifnot( 'mapframe' %in% class(loc) )
    stopifnot( isBOOL(expandtomarkers) )
    
    # Get map unit.
    map.unit <- getMapUnit(x)
    
    if ( is.na(map.unit) ) {
        stop("map unit not found")
    }
    
    if ( getMapUnit(loc) != map.unit ) {
        stop("map unit mismatch")
    }
    
    if ( nrow(x) == 0 ) {
        stop("cannot find flanking row indices in mapframe with zero rows")
    }
    
    # Get sequences corresponding to individual map loci.
    x.seqs <- pullLocusSeq(x)
    
    # Get locus mapframe sequence.
    loc.seq <- unique( pullLocusSeq(loc) )
    
    # Check that loc mapframe refers to one sequence.
    # NB: implicitly checks that loc has one or more rows.
    stopifnot( length(loc.seq) == 1 )
    
    # Get row indices of sequence loci.
    indices <- which( x.seqs == loc.seq )
    stopifnot( length(indices) > 0 )
    
    # Get map positions for this sequence.
    seq.pos <- x$pos[indices]

    exrange <- loc$pos[ ! sapply( loc$pos, inRange, range(seq.pos) ) ]
    if ( length(exrange) > 0 ) {
        stop("locus positions out of range - ", toString(exrange), "'")
    }
    
    if ( expandtomarkers && hasRownames(x) ) {
        # If expanding to markers and mapframe has locus IDs, create
        # a mask in which markers and sequence endpoints are TRUE..
        loc.mask <- isMarkerID(rownames(x)[indices]) |
            getIndices(seq.pos) %in% c( 1, length(seq.pos) )
    } else {
        # ..otherwise create mask with all TRUE.
        loc.mask <- rep_len( TRUE, length(seq.pos) )
    }
    
    # Get lower flanking locus index.
    lower.pos <- loc$pos[1]
    lower.flanking <- seq.pos[ loc.mask & lower.pos >= seq.pos ]
    lower.diff <- lower.pos - lower.flanking
    lower.indices <- which( lower.diff == min(lower.diff) )
    lower.index <- lower.indices[ length(lower.indices) ] + min(indices) - 1
    
    # Get upper flanking locus index.
    upper.pos <- loc$pos[ nrow(loc) ]
    upper.flanking <- seq.pos[ loc.mask & seq.pos >= upper.pos ]
    upper.diff <- upper.flanking - upper.pos
    upper.indices <- which( upper.diff == min(upper.diff) )
    upper.index <- upper.indices[1] + min(indices) - 1 + length(indices) - length(upper.flanking)
    
    return( c(lower.index, upper.index) )
}

# findLoci ---------------------------------------------------------------------
#' Find closest map loci.
#' 
#' @param x Object containing map data.
#' @param loc Locus \code{mapframe} specifying map positions.
#' 
#' @return A \code{mapframe} whose rows each contain the closest locus to the
#' map position specified in the corresponding row of the locus \code{mapframe}.
#' 
#' @export
#' @family map utilities
#' @rdname findLoci
findLoci <- function(x, loc) {
    
    stopifnot( 'mapframe' %in% class(loc) )
    
    prox.loci <- data.frame( matrix( ncol=2, nrow=nrow(loc),
        dimnames=list(NULL, c('chr', 'pos')) ) )
    
    # Find closest locus to each map position.
    for ( i in getRowIndices(loc) ) {
        
        # Get loci flanking this map position.
        flanking <- findFlanking(x, loc[i, ])
       
        # Get closest locus, taking the first in case of ties.
        prox.loci[i, ] <- flanking[ which.min( abs(loc$pos[i] - flanking$pos) ), ]
    }
    
    return( as.mapframe(prox.loci, getMapUnit(loc)) )
}

# findLocusRowIndices ----------------------------------------------------------
#' Find row index of closest \code{mapframe} locus.
#' 
#' @param x A \code{mapframe} object.
#' @param loc Locus \code{mapframe} specifying map positions.
#' 
#' @return Vector of row indices, each element of which indicates the row of the
#' target \code{mapframe} object containing the closest locus to the map position
#' specified in the corresponding row of the locus \code{mapframe}.
#' 
#' @keywords internal
#' @rdname findLocusRowIndices
findLocusRowIndices <- function(x, loc) {
    
    stopifnot( 'mapframe' %in% class(loc) )
    
    prox.indices <- integer( nrow(loc) )
    
    # Find row indices of closest locus to each specified map position.
    for ( i in getRowIndices(loc) ) {
        
        # Get row indices of flanking loci.
        row.indices <- findFlankingRowIndices(x, loc[i, ])
        
        # Get row index of closest locus, taking the first in case of ties.
        prox.indices[i] <- row.indices[ which.min( abs(loc$pos[i] - x$pos[row.indices]) ) ]
    }
    
    return(prox.indices)
}

# findMarkers ------------------------------------------------------------------
#' Find closest markers.
#' 
#' @param x Object containing map data.
#' @param loc Locus \code{mapframe} specifying map positions.
#' 
#' @return A \code{mapframe} whose rows each contain the closest marker to the
#' map position specified in the corresponding row of the locus \code{mapframe}.
#' 
#' @export
#' @family map utilities
#' @rdname findMarkers
findMarkers <- function(x, loc) {
    findLoci(extractMarkers(x), loc)
}

# getSeqColIndex ---------------------------------------------------------------
#' Get sequence column index.
#' 
#' Get sequence column index for a \code{mapframe} or equivalent
#' \code{data.frame}. The sequence column is taken to be the column
#' whose heading is 'chr'. An error is raised if there is not exactly
#' one sequence column.
#' 
#' @param x A \code{mapframe} or equivalent \code{data.frame}.
#' 
#' @return Sequence column index.
#' 
#' @keywords internal
#' @rdname getSeqColIndex
getSeqColIndex <- function(x) {
    
    stopifnot( is.data.frame(x) )
    stopif( anyDuplicated( colnames(x) ) )
    
    seqcol.indices <- which( colnames(x) == 'chr' )
    
    if ( length(seqcol.indices) == 1 ) {
        seqcol.index <- seqcol.indices[1]
    } else if ( length(seqcol.indices) == 0 ) {
        stop("no sequence columns found")
    } else {
        stop("multiple sequence columns found")
    }
    
    return(seqcol.index)
}

# getDatColIndices -------------------------------------------------------------
#' Get data column indices.
#' 
#' Get data column indices for a \code{mapframe} or equivalent \code{data.frame}.
#' Data columns are taken to be any columns that aren't a sequence or position
#' column.
#' 
#' @param x A \code{mapframe} or equivalent \code{data.frame}.
#' @param datcolumns Optional parameter indicating which data columns to
#' consider. This must be either a character vector of data column names
#' or a numeric vector of indices \emph{with respect to the set of data
#' columns}. If no data columns are specified, all are considered.
#' @param strict Option indicating that the specified data columns must
#' be in the same order as they are in the \code{mapframe} object.
#' 
#' @return Vector of data column indices.
#' 
#' @keywords internal
#' @rdname getDatColIndices
getDatColIndices <- function(x, datcolumns=NULL, strict=FALSE) {
    
    stopifnot( is.data.frame(x) )
    stopif( anyDuplicated( colnames(x) ) )
    stopifnot( isBOOL(strict) )
    
    column.indices <- getColIndices(x)
    seqcol.indices <- getSeqColIndex(x)
    poscol.indices <- getPosColIndex(x)
    datcol.indices <- column.indices[ column.indices != seqcol.indices &
        column.indices != poscol.indices ]
    
    if ( ! is.null(datcol.indices) ) {
        
        if ( ! is.null(datcolumns) ) {
            
            stopif( anyNA(datcolumns) )
            
            if ( isWholeNumber(datcolumns) ) {
                
                available.datcolumns <- getIndices(datcol.indices)
                
                exrange <- datcolumns[ ! inRange( datcolumns, range(available.datcolumns) ) ]
                if ( length(exrange) > 0 ) {
                    stop("data column indices out of range - '", toString(exrange), "'")
                }
                
                datcol.indices <- datcol.indices[ datcolumns ]
                
            } else if ( is.character(datcolumns) ) {
                
                available.headings <- colnames(x)[datcol.indices]
                
                unfound <- datcolumns[ ! datcolumns %in% available.headings ]
                if ( length(unfound) > 0 ) {
                    stop("data column names not found - '", toString(unfound), "'")
                }
                
                datcol.indices <- datcol.indices[ match(datcolumns, available.headings) ]
                
            } else {
                
                stop("data columns must be specified by index or name")
            }
            
            if (strict) { # NB: also ensures no duplicates
                if ( is.unsorted(datcol.indices, strictly=TRUE) ) {
                    stop("data columns not specified in strictly increasing order")
                }
            }
        }
        
    } else {
        
        datcol.indices <- integer()
    }
    
    return(datcol.indices)
}

# getMapSteps ------------------------------------------------------------------
#' Get steps between loci.
#' 
#' @param x Object containing map data.
#' 
#' @return List of vectors, each giving the steps between markers on a given
#' sequence.
#' 
#' @export
#' @family map utilities
#' @rdname getMapSteps
getMapSteps <- function(x) {
    UseMethod('getMapSteps', x)
}

# getMapSteps.mapframe ---------------------------------------------------------
#' @export
#' @rdname getMapSteps
getMapSteps.mapframe <- function(x) {
    
    # Get sequences corresponding to individual map loci.
    x.seqs <- pullLocusSeq(x)
    
    # Get individual locus positions.
    x.pos <- pullLocusPos(x)
    
    # Get map sequences.
    map.seqs <- unique(x.seqs)
    
    # Get map steps.
    map.steps <- lapply(map.seqs, function(map.seq)
        diff(x.pos[ x.seqs == map.seq ]) )
    names(map.steps) <- map.seqs
    
    return(map.steps)
}

# getMapSteps.map --------------------------------------------------------------
#' @export
#' @rdname getMapSteps
getMapSteps.map <- function(x) {
    return( lapply(x, function(seq.pos) diff( as.numeric(seq.pos) ) ) )
}

# getMapSteps.scanone ----------------------------------------------------------
#' @export
#' @rdname getMapSteps
getMapSteps.scanone <- function(x) {
    return( getMapSteps.mapframe(x) )
}

# getMapUnit -------------------------------------------------------------------
#' Get map unit.
#' 
#' Get the map unit of the object. In most cases, this is taken from the
#' \code{'map.unit'} attribute of the object. For a \code{data.frame}, the
#' position column is also checked for map unit information; an error is
#' raised if \code{data.frame} map unit information is inconsistent.
#' 
#' @template mapunits
#'  
#' @param x Object containing map data.
#' 
#' @return Map unit. Returns NA if map unit information not found.
#' 
#' @export
#' @family map utilities
#' @rdname getMapUnit
getMapUnit <- function(x) {
    UseMethod('getMapUnit', x)
}

# getMapUnit.data.frame --------------------------------------------------------
#' @export
#' @rdname getMapUnit
getMapUnit.data.frame <- function(x) {
    
    map.unit <- NA_character_
    
    existing.mapunit <- getMapUnit.default(x)
    heading.mapunit <- getPosColNameMapUnit(x)
    positions.mapunit <- getPosColDataMapUnit(x)
    
    for ( cue in list(existing.mapunit, heading.mapunit, positions.mapunit) ) {
        if ( ! is.null(cue) && ! is.na(cue) ) {
            if ( is.na(map.unit) ) {
                map.unit <- cue
            } else if ( cue != map.unit ) {
                stop("data frame has inconsistent map unit information")
            }
        }
    }
    
    return(map.unit)
}

# getMapUnit.default -----------------------------------------------------------
#' @export
#' @rdname getMapUnit
getMapUnit.default <- function(x) {
    
    map.unit <- attr(x, 'map.unit')
    
    if ( is.null(map.unit) ) {
        map.unit <- NA_character_
    }
    
    return(map.unit)
}

# getMapUnit.map ---------------------------------------------------------------
#' @export
#' @rdname getMapUnit
getMapUnit.map <- function(x) {
    
    map.unit <- getMapUnit.default(x)
    
    if ( is.na(map.unit) ) {
        map.unit <- 'cM'
    }
    
    return(map.unit)
}

# getMapUnit.mapframe ----------------------------------------------------------
#' @export
#' @rdname getMapUnit
getMapUnit.mapframe <- function(x) {
    return( getMapUnit.default(x) )
}

# getMapUnit.scanone -----------------------------------------------------------
#' @export
#' @rdname getMapUnit
getMapUnit.scanone <- function(x) {
    
    map.unit <- getMapUnit.default(x)
    
    if ( ! is.na(map.unit) && map.unit != 'cM' ) {
        stop("scanone result has unexpected map unit - '", map.unit, "'")
    }
    
    return('cM')
}

# getMapUnitSuffix -------------------------------------------------------------
#' Get map unit suffix from a vector of map positions.
#' 
#' @template mapunits
#' 
#' @param x A character vector of map positions.
#' 
#' @return Map unit indicated by the map positions. Returns NA if the map
#' position vector does not contain map unit suffixes.
#' 
#' @keywords internal
#' @rdname getMapUnitSuffix
getMapUnitSuffix <- function(x) {
    
    stopifnot( is.character(x) )
    
    map.unit <- NA_character_
    map.type <- NULL
    
    if ( length(x) > 0 ) {
        
        positions <- stripWhite( unname(x) )
        
        # Create genetic pos matrix, in which rows correspond to
        # loci and columns correspond to genetic map patterns.
        pos.matrix <- vapply( const$pattern$gmap, function(p)
            sapply(positions, grepl, pattern=p, ignore.case=TRUE),
            logical( length(positions) ) )
        
        # If genetic pos matrix was simplified to a vector (as
        # for a single map position), restore to a matrix.
        if ( ! is.matrix(pos.matrix) ) {
            pos.matrix <- matrix(pos.matrix, nrow=length(positions),
                dimnames=list(positions, const$pattern$gmap) )
        }
        
        # For each locus, get column indices of
        # first matching genetic map pattern.
        unit.indices <- unique( sapply(getIndices(positions), function(i)
            match(TRUE, pos.matrix[i, ]) ) )
        
        if ( length(unit.indices) == 1 ) {
            
            if ( ! is.na(unit.indices) ) {
                map.type <- 'gmap'
            }
        
        } else {
            
            stop("map positions have inconsistent/incomplete map unit information")
        }
        
        if ( is.null(map.type) ) {
            
            # Create physical pos matrix, in which rows correspond to
            # loci and columns correspond to physical map patterns.
            pos.matrix <- vapply( const$pattern$pmap, function(p)
                sapply(positions, grepl, pattern=p, ignore.case=TRUE),
                logical( length(positions) ) )
            
            # If physical pos matrix was simplified to a vector
            # (as for a single map position), restore to a matrix.
            if ( ! is.matrix(pos.matrix) ) {
                pos.matrix <- matrix(pos.matrix, nrow=length(positions),
                    dimnames=list(positions, const$pattern$pmap) )
            }
            
            # For each locus, get column indices of first matching genetic
            # map pattern. NB: must check in order (Mb > kb > bp).
            unit.indices <- unique( sapply(getIndices(positions), function(i)
                match(TRUE, pos.matrix[i, ])) )
            
            if ( length(unit.indices) == 1 ) {
                
                if ( ! is.na(unit.indices) ) {
                    map.type <- 'pmap'
                }
                
            } else {
                stop("map positions have inconsistent/incomplete map unit information")
            }
        }
        
        # Update map unit, if identified.
        if ( ! is.null(map.type) ) {
            map.unit <- rownames(const$map.info[[map.type]])[ unit.indices[1] ]
        }
    }
    
    return(map.unit)
}

# getPosColDataMapUnit ---------------------------------------------------------
#' Get map unit from position column values.
#' 
#' @template mapunits
#' 
#' @param x A \code{mapframe} or equivalent \code{data.frame},
#' or a vector of map positions.
#' 
#' @return Map unit indicated by the position column data. Returns NA if map
#' unit information not found in position column data.
#' 
#' @keywords internal
#' @rdname getPosColDataMapUnit
getPosColDataMapUnit <- function(x) {
    UseMethod('getPosColDataMapUnit', x)
}

# getPosColDataMapUnit.data.frame ----------------------------------------------
#' @rdname getPosColDataMapUnit
getPosColDataMapUnit.data.frame <- function(x) {
    return( getMapUnitSuffix( as.character(x[, getPosColIndex(x)]) ) )
}

# getPosColDataMapUnit.character -----------------------------------------------
#' @rdname getPosColDataMapUnit
getPosColDataMapUnit.character <- function(x) {
    return( getMapUnitSuffix(x) )
}

# getPosColDataMapUnit.numeric -------------------------------------------------
#' @rdname getPosColDataMapUnit
getPosColDataMapUnit.numeric <- function(x) {
    return(NA_character_)
}

# getPosColIndex ---------------------------------------------------------------
#' Get \code{mapframe} position column index.
#' 
#' Get position column index for a \code{mapframe} or equivalent \code{data.frame}.
#' The position column is taken to be the leftmost column whose heading contains
#' the word 'pos'. This can also be in uppercase (i.e. 'POS'), but can't be part
#' of a larger word, such as 'position'. An error is raised if a position column
#' cannot be found.
#' 
#' @param x A \code{mapframe} or equivalent \code{data.frame}.
#' 
#' @return Position column index.
#' 
#' @keywords internal
#' @rdname getPosColIndex
getPosColIndex <- function(x) {
    
    stopifnot( is.data.frame(x) )
    stopif( anyDuplicated( colnames(x) ) )
    
    poscol.indices <- which( grepl(const$pattern$poscol, colnames(x),
        ignore.case=TRUE) )
    
    if ( length(poscol.indices) == 0 ) {
        stop("pos column not found")
    }
   
    return(poscol.indices[1])
}

# getPosColNameMapUnit ---------------------------------------------------------
#' Get map unit from position column heading.
#' 
#' @template mapunits
#' 
#' @param x A \code{mapframe} or equivalent \code{data.frame},
#' or a position column heading.
#' 
#' @return Map unit indicated by the position column heading. Returns NA if map
#' unit information not found in position column heading.
#' 
#' @keywords internal
#' @rdname getPosColNameMapUnit
getPosColNameMapUnit <- function(x) {
    UseMethod('getPosColNameMapUnit', x)
}

# getPosColNameMapUnit.data.frame ----------------------------------------------
#' @rdname getPosColNameMapUnit
getPosColNameMapUnit.data.frame <- function(x) {
    return( getPosColNameMapUnit(colnames(x)[ getPosColIndex(x) ]) )
}

# getPosColNameMapUnit.character -----------------------------------------------
#' @rdname getPosColNameMapUnit
getPosColNameMapUnit.character <- function(x) {
    
    stopifnot( length(x) == 1 )
    stopifnot( grepl(const$pattern$poscol, x, ignore.case=TRUE) )
    
    map.unit <- NA_character_
    
    poscol.words <- unlist( strsplit(x, '[^[:alpha:]]+') )
    
    mapunit.mask <- const$known.map.units %in% poscol.words
    
    matching.map.units <- const$known.map.units[mapunit.mask]
    
    if ( length(matching.map.units) == 1 ) {
        map.unit <- matching.map.units[1]
    } else if ( length(matching.map.units) > 1 ) {
        stop("pos column heading contains multiple map units - '", x, "'")
    } 
    
    return(map.unit)
}

# gmapframe --------------------------------------------------------------------
#' Create a new genetic \code{mapframe}.
#' 
#' @template mapframe
#' 
#' @param ... Further arguments. These are passed to the \code{data.frame}
#' constructor. All arguments should be passed as keyword arguments.
#' 
#' @return A new \code{mapframe} with genetic map positions.
#' 
#' @export
#' @family map utilities
#' @rdname gmapframe
gmapframe <- function(...) {
    return( mapframe(map.unit='cM', ...) )
}

# inferMapStep -----------------------------------------------------------------
#' Infer map step size.
#' 
#' @param x Object containing map data.
#' @param tol Tolerance for step equality.
#' 
#' @return Inferred map step. Returns \code{NULL} if map step could not be
#' inferred.
#' 
#' @export
#' @family map utilities
#' @rdname inferMapStep
inferMapStep <- function(x, tol=1e-5) {
    map.steps <- unlist( getMapSteps(x) )
    map.step <- inferStepSize(map.steps)
    return(map.step)
}

# inMapOrder -------------------------------------------------------------------
#' Test if object is in map order.
#' 
#' @param x Object containing map data.
#' 
#' @return TRUE if map data of object is ordered by sequence, then map
#' position; FALSE otherwise.
#' 
#' @export
#' @family map utilities
#' @rdname inMapOrder
inMapOrder <- function(x) {
    UseMethod('inMapOrder', x)
}

# inMapOrder.data.frame --------------------------------------------------------
#' @export
#' @rdname inMapOrder
inMapOrder.data.frame <- function(x) {
    
    if ( nrow(x) > 0 ) {
    
        x.seqs <- pullLocusSeq(x)
        x.pos <- pullLocusPos(x)
        
        if ( is.unsorted( orderSeq(x.seqs) ) ) {
            return(FALSE)
        }
        
        map.seqs <- unique(x.seqs)
        
        for ( map.seq in map.seqs ) {
            if ( is.unsorted(x.pos[ x.seqs == map.seq ]) ) {
                return(FALSE)
            }
        }
    }
    
    return(TRUE)
}

# inMapOrder.map ---------------------------------------------------------------
#' @export
#' @rdname inMapOrder
inMapOrder.map <- function(x) {
    
    if ( is.unsorted( orderSeq( names(x) ) ) ) {
        return(FALSE)
    }
    
    if ( any( sapply(x, is.unsorted) ) ) {
        return(FALSE)
    }
    
    return(TRUE)
}

# intersectLoci ----------------------------------------------------------------
#' Find intersection set of loci.
#' 
#' @param ... One or more objects containing map data.
#' 
#' @return Locus \code{mapframe} with loci common to all input objects.
#' 
#' @keywords internal
#' @rdname intersectLoci
intersectLoci <- function(...) {
    
    loc.list <- list(...)
    
    stopifnot( length(loc.list) > 0 )
    
    intersection <- pullLoci(loc.list[[1]])
    
    if ( length(loc.list) > 1 ) {
        
        for ( i in 2:length(loc.list) ) {
            intersection <- matchLoci(intersection, loc.list[[i]])
        }
    }
    
    return(intersection)
}

# isGeneticMap -----------------------------------------------------------------
#' Test if object is a genetic \code{map}.
#' 
#' @param x Test object.
#' 
#' @return TRUE if object is a genetic \code{map} object; FALSE otherwise.
#' 
#' @keywords internal
#' @rdname isGeneticMap
isGeneticMap <- function(x) {
    
    if ( 'map' %in% class(x) ) {
        map.unit <- getMapUnit(x)
        if ( ! is.na(map.unit) && const$known.map.types[map.unit] == 'gmap' ) {
            return(TRUE)
        }
    }
    
    return(FALSE)
}

# isGeneticMapframe ------------------------------------------------------------
#' Test if object is a genetic \code{mapframe}.
#' 
#' @param x Test object.
#' 
#' @return TRUE if object is a \code{mapframe} whose pos column contains
#' genetic map positions; FALSE otherwise.
#' 
#' @keywords internal
#' @rdname isGeneticMapframe
isGeneticMapframe <- function(x) {
    
    if ( 'mapframe' %in% class(x) ) {
        map.unit <- getMapUnit(x)
        if ( ! is.na(map.unit) && const$known.map.types[map.unit] == 'gmap' ) {
            return(TRUE)
        }
    }
    
    return(FALSE)
}

# isGeneticMapUnit -------------------------------------------------------------
#' Test if object is a known genetic map unit.
#' 
#' @param x Test object.
#' 
#' @return TRUE if object is a known genetic map unit; FALSE otherwise.
#' 
#' @keywords internal
#' @rdname isGeneticMapUnit
isGeneticMapUnit <- function(x) {
    tryCatch({
        validateMapUnit(x, map.type='gmap')
    }, error=function(e) {
        return(FALSE)
    })
    return(TRUE)
}

# isPhysicalMap ----------------------------------------------------------------
#' Test if object is a physical \code{map}.
#' 
#' @param x Test object.
#' 
#' @return TRUE if object is a physical \code{map}; FALSE otherwise.
#' 
#' @keywords internal
#' @rdname isPhysicalMap
isPhysicalMap <- function(x) {
    
    if ( 'map' %in% class(x) ) {
        map.unit <- getMapUnit(x)
        if ( ! is.na(map.unit) && const$known.map.types[map.unit] == 'pmap' ) {
            return(TRUE)
        }
    }
    
    return(FALSE)
}

# isPhysicalMapframe -----------------------------------------------------------
#' Test if object is a physical \code{mapframe}.
#' 
#' @param x Test object.
#' 
#' @return TRUE if object is a \code{mapframe} whose pos column contains
#' physical map positions; FALSE otherwise.
#' 
#' @keywords internal
#' @rdname isPhysicalMapframe
isPhysicalMapframe <- function(x) {
    
    if ( 'mapframe' %in% class(x) ) {
        map.unit <- getMapUnit(x)
        if ( ! is.na(map.unit) && const$known.map.types[map.unit] == 'pmap' ) {
            return(TRUE)
        }
    }
    
    return(FALSE)
}

# isPhysicalMapUnit ------------------------------------------------------------
#' Test if object is a known physical map unit.
#' 
#' @param x Test object.
#' 
#' @return TRUE if object is a known physical map unit; FALSE otherwise.
#' 
#' @keywords internal
#' @rdname isPhysicalMapUnit
isPhysicalMapUnit <- function(x) {
    
    tryCatch({
        validateMapUnit(x, map.type='pmap')
    }, error=function(e) {
        return(FALSE)
    })
    
    return(TRUE)
}

# isValidMapUnit ---------------------------------------------------------------
#' Test if object is a known map unit.
#' 
#' @param x Test object.
#' 
#' @return TRUE if object is a known genetic or physical map unit;
#' FALSE otherwise.
#' 
#' @keywords internal
#' @rdname isValidMapUnit
isValidMapUnit <- function(x) {
    
    tryCatch({
        validateMapUnit(x)
    }, error=function(e) {
        return(FALSE)
    })
    
    return(TRUE)
}

# makePlaceholderMap -----------------------------------------------------------
#' Make a placeholder map.
#' 
#' Make a placeholder map from the given combination of locus sequences, map
#' unit, and map step size. If locus IDs are specified, the locus IDs of the
#' placeholder map will be set from these. Otherwise, locus IDs are set to
#' pseudomarker IDs or default marker IDs, according to whether map units
#' are for a genetic or physical map, respectively.
#' 
#' Locus sequences and locus IDs are assumed to be ordered. As in R/qtl, any
#' sequence with too few loci is extended, to ensure that all sequences have
#' enough loci.
#' 
#' @param locus.seqs Sequences corresponding to individual loci.
#' @param locus.ids Individual locus IDs.
#' @param map.unit Map unit.
#' @param step Map step size.
#' 
#' @keywords internal
#' @rdname makePlaceholderMap
makePlaceholderMap <- function(locus.seqs, locus.ids=NULL, map.unit='cM',
    step=5) {
    
    stopifnot( is.character(locus.seqs) )
    stopifnot( length(locus.seqs) > 0 )
    stopifnot( isValidMapUnit(map.unit) )
    stopifnot( isSinglePositiveNumber(step) )
    
    # Get normalised locus sequences.
    norm.locus.seqs <- normSeq(locus.seqs)
    
    # Test if map unit is for a genetic map.
    is.gmap <- isGeneticMapUnit(map.unit)
    
    # Get run-length encoding of normalised locus sequences.
    runs <- rle(norm.locus.seqs)
    
    # Check for any sequence label that is
    # not grouped together in a single run.
    if ( anyDuplicated(runs$values) ) {
        stop("map sequences are not ordered")
    }
    
    # Get placeholder-map sequences.
    map.seqs <- runs$values
    
    # Get index list for locus sequences.
    # NB: groups indices by sequence.
    index.list <- getRunIndexList(norm.locus.seqs)
    
    # Generate placeholder locus positions.
    locus.pos <- unlist( lapply( unname(index.list), function(indices)
        seq(0, length=length(indices), by=step)) )
    
    # Create placeholder mapframe.
    placeholder <- mapframe(chr=norm.locus.seqs, pos=locus.pos,
        map.unit=map.unit)
    
    # If locus IDs specified, set from those..
    if ( ! is.null(locus.ids) ) {
        
        stopifnot( is.character(locus.ids) )
        stopifnot( length(locus.ids) == length(locus.seqs) )
        
        rownames(placeholder) <- locus.ids
        
    } else { # ..otherwise generate default locus IDs.
        
        if (is.gmap) {
            rownames(placeholder) <- makePseudomarkerIDs(placeholder)
        } else {
            rownames(placeholder) <- makeDefaultMarkerIDs(placeholder)
        }
    }
    
    # Check placeholder mapframe for short sequences.
    short.seqs <- map.seqs[ sapply(getRunIndices(runs),
        function(i) runs$lengths[i] < const$min.lps) ]
    
    # If placeholder mapframe contains short sequences,
    # extend these until they have enough loci.
    if ( length(short.seqs) > 0 ) {
        
        short.flank.seqs <- character()
        short.flank.pos <- character()
        
        # Get flanking locus sequences and
        # positions for each short sequence.
        for ( short.seq in short.seqs ) {
            
            # Get loci in this sequence.
            loc <- placeholder[ placeholder$chr == short.seq, ]
            
            # Calculate shortfall in number of loci.
            shortfall <- const$min.lps - nrow(loc)
            
            # Get number of flanking steps needed to get enough loci.
            num.steps <- (shortfall + 1) %/% 2
            
            # Generate sequence of flanking positions in each direction.
            lower.flank.pos <- seq(from=loc$pos[1] - step, by=-step,
                length.out=num.steps)
            upper.flank.pos <- seq(from=loc$pos[ nrow(loc) ] + step,
                by=step, length.out=num.steps)
            
            # Set flanking locus info for this sequence.
            flank.pos <- c(lower.flank.pos, upper.flank.pos)
            flank.seqs <- rep(loc$chr[1], length(flank.pos))
            
            # Append to flanking locus info for all short sequences.
            short.flank.seqs <- c(short.flank.seqs, flank.seqs)
            short.flank.pos <- c(short.flank.pos, flank.pos)
        }
        
        # Create mapframe of flanking loci.
        flanking <- mapframe(chr=short.flank.seqs, pos=short.flank.pos,
            map.unit=map.unit)
        
        # Generate default locus IDs.
        if (is.gmap) {
            rownames(flanking) <- makePseudomarkerIDs(flanking)
        } else {
            rownames(flanking) <- makeDefaultMarkerIDs(flanking)
        }
        
        placeholder <- rbind(placeholder, flanking)
        
        placeholder <- orderMap(placeholder)
    }
    
    return( as.map(placeholder) )
}

# makeMapFromDefaultMarkerIDs --------------------------------------------------
#' Make map from default marker IDs.
#' 
#' @template map
#' 
#' @param marker.ids Vector of default marker IDs.
#' 
#' @return An \pkg{R/qtl} \code{map} object.
#' 
#' @export
#' @keywords internal
#' @rdname makeMapFromDefaultMarkerIDs
makeMapFromDefaultMarkerIDs <- function(marker.ids) {
    return( as.map( parseDefaultMarkerIDs(marker.ids) ) )
}

# makeMapFromDefaultQTLNames ---------------------------------------------------
#' Make map from default QTL names.
#' 
#' @template map
#' 
#' @param qtl.names Vector of default QTL names.
#' 
#' @return An \pkg{R/qtl} \code{map} object.
#' 
#' @export
#' @keywords internal
#' @rdname makeMapFromDefaultQTLNames
makeMapFromDefaultQTLNames <- function(qtl.names) {
    return( as.map( parseDefaultQTLNames(qtl.names) ) )
}

# makeMapFromPseudomarkerIDs ---------------------------------------------------
#' Make map from pseudomarker IDs.
#' 
#' @template map
#' 
#' @param loc.ids Vector of pseudomarker IDs.
#' 
#' @return An \pkg{R/qtl} \code{map} object.
#' 
#' @export
#' @keywords internal
#' @rdname makeMapFromPseudomarkerIDs
makeMapFromPseudomarkerIDs <- function(loc.ids) {
    return( as.map( parsePseudomarkerIDs(loc.ids) ) )
}

# mapframe ---------------------------------------------------------------------
#' Create a new \code{mapframe}.
#' 
#' @template mapframe
#' 
#' @param map.unit Map unit for the new \code{mapframe}.
#' @param ... Further arguments. These are passed to the \code{data.frame}
#' constructor. All arguments should be passed as keyword arguments.
#' 
#' @return A new \code{mapframe}.
#' 
#' @export
#' @family map utilities
#' @rdname mapframe
mapframe <- function(map.unit=NULL, ...) {
    
    stopifnot( allKwargs(...) )
    
    args <- list(...)
    
    mapframe.colnames <- names(args)[ ! names(args) %in% const$dataframe.keywords ]
    
    prev.state <- getOption('stringsAsFactors')
    options(stringsAsFactors=FALSE)
    
    if ( length(mapframe.colnames) == 0 ) {
        x <- data.frame(chr=character(), pos=numeric(), ...)
    } else {
        x <- data.frame(...)
    }
    
    options(stringsAsFactors=prev.state)
    
    return( as.mapframe(x, map.unit=map.unit) )
}

# mapkey -----------------------------------------------------------------------
#' Make new \code{mapkey} rescaling function.
#' 
#' This function takes as input one or more maps - one of each supported
#' map type - and returns a function that can be used to rescale a locus
#' \code{mapframe} of map positions. The input maps must be comparable,
#' having the same set of markers on the same chromosomes, although the
#' marker names may be different. As in \pkg{R/qtl}, map positions are
#' interpolated if they are within the bounds of the component maps, and
#' extrapolated otherwise.
#' 
#' @param ... One or more component maps, all with comparable marker loci.
#' 
#' @return A \code{mapkey} function similar to the \pkg{R/qtl} function
#' \code{interpPositions}, taking a locus \code{mapframe} and map unit
#' parameter, and returning a \code{mapframe} that has been rescaled in
#' terms of the given map unit.
#' 
#' @export
#' @rdname mapkey
mapkey <- function(...) {
    
    stopif( anyKwargs(...) )
    
    maps <- list(...)
    
    # If component maps specified, generate mapkey object..
    if ( length(maps) > 0 ) {
        
        # Ensure all maps are valid map objects,
        # converting to a map object if necessary.
        maps <- lapply(maps, as.map)
        
        # Ensure all maps have normalised sequence labels.
        maps <- lapply(maps, normSeq)
        
        # Get map units.
        map.units <- unname( sapply(maps, getMapUnit) )
        
        if ( anyNA(map.units) ) {
            stop("mapkey component maps must have valid map units")
        }
        
        # TODO: relax limitation of maps to basic map units.
        if ( ! all(map.units %in% const$basic.map.unit) ) {
            stop("mapkey component maps must be in basic map units")
        }
        
        # Get map types from map unit information.
        map.types <- unname( sapply(map.units, function(map.unit)
            const$known.map.types[[map.unit]]) )
        
        # Check that each map type is represented only once.
        stopif( anyDuplicated(map.types) )
        
        # Set mapkey names from map units.
        names(maps) <- map.units # NB: assumes all maps use different units
        
        # Check that all maps contain the same sequences in the same order.
        if ( ! all( sapply(maps[-1], function(m)
            identical(names(m), names(maps[[1]])) ) ) ) {
            stop("sequence mismatch between mapkey component maps")
        }
        
        # Check that all maps have a consistent number of markers in each sequence.
        if ( ! all( sapply(maps[-1], function(m)
            identical(lengths(m), lengths(maps[[1]])) ) ) ) {
            stop("marker count mismatch between mapkey component maps")
        }
        
        map.key <- function(loc, map.unit) {
            
            stopifnot( 'mapframe' %in% class(loc) )
            validateMapUnit(map.unit)
            
            seqcol.index <- getSeqColIndex(loc)
            poscol.index <- getPosColIndex(loc)
            
            # Get old and new map units.
            old.map.unit <- getMapUnit(loc)
            validateMapUnit(old.map.unit)
            new.map.unit <- map.unit
            
            # Get old and new map types.
            old.map.type <- const$known.map.types[old.map.unit]
            new.map.type <- const$known.map.types[new.map.unit]
            
            # Get basic units of old and new map types.
            old.basic.unit <- const$basic.map.unit[old.map.type]
            new.basic.unit <- const$basic.map.unit[new.map.type]
            
            stopifnot( old.basic.unit %in% names(maps) )
            stopifnot( new.basic.unit %in% names(maps) )
            
            # Get old and new maps.
            oldmap <- maps[[old.basic.unit]]
            newmap <- maps[[new.basic.unit]]
            
            # Check that all sequences are in mapkey.
            loc.seqs <- unique(loc[, seqcol.index])
            missing.seqs <- loc.seqs[ ! loc.seqs %in% names(oldmap) ]
            if ( length(missing.seqs) > 0 ) {
                stop("sequences not found in mapkey - '", toString(missing.seqs), "'")
            }
            
            # Convert each locus position.
            for ( i in getRowIndices(loc) ) {
                
                # Get this locus position.
                x <- loc[i, poscol.index]
                
                # Convert locus position to basic units, if needed.
                if ( old.map.unit != old.basic.unit ) {
                    x <- x * const$map.info[[old.map.type]][old.map.unit, 'factor']
                }
                
                # Convert locus position from one map type to another, if needed.
                if ( old.basic.unit != new.basic.unit ) {
                    
                    locus.seq <- loc[i, seqcol.index]
                    old.pos <- oldmap[[locus.seq]]
                    new.pos <- newmap[[locus.seq]]
                    
                    # If locus is within the range of the old map,
                    # interpolate its position on the new map..
                    if ( inRange(x, range(old.pos)) ) {
                        
                        # Get index of lower flanking locus.
                        xld <- x - old.pos
                        lower.indices <- which( xld == min( xld[ xld >= 0 ] ) )
                        lower.index <- lower.indices[ length(lower.indices) ]
                        
                        # Get index of upper flanking locus.
                        xud <- old.pos - x
                        upper.indices <- which( xud == min( xud[ xud >= 0 ] ) )
                        upper.index <- upper.indices[1]
                        
                        # Get flanking loci in old and new maps.
                        x1 <- old.pos[lower.index]
                        x2 <- old.pos[upper.index]
                        y1 <- new.pos[lower.index]
                        y2 <- new.pos[upper.index]
                        
                        # Interpolate between flanking loci.
                        x <- y1 + (x - x1) * (y2 - y1) / (x2 - x1)
                        
                    } else { # ..otherwise extrapolate its position from map endpoints.
                        
                        # Get terminal loci in old and new maps.
                        x1 <- old.pos[1]
                        x2 <- old.pos[ length(old.pos) ]
                        y1 <- new.pos[1]
                        y2 <- new.pos[ length(new.pos) ]
                        
                        # Extrapolate from map endpoints.
                        if ( x < x1 ) {
                            x <- y1 - (x1 - x) * (y2 - y1) / (x2 - x1)
                        } else { # x > x2
                            x <- y2 + (x - x2) * (y2 - y1) / (x2 - x1)
                        }
                    }
                }
                
                # Convert locus position from basic units, if needed.
                if ( new.map.unit != new.basic.unit ) {
                    x <- x / const$map.info[[new.map.type]][new.map.unit, 'factor']
                }
                
                loc[i, poscol.index] <- unname(x)
            }
            
            attr(loc, 'map.unit') <- map.unit
            
            return(loc)
        }
        
        class(map.key) <- c('mapkey', 'function')

    } else { # ..otherwise return default mapkey object.
        
        # Get name of current genome.
        genome <- genomeOpt()
        
        # Get default mapkey of current genome.
        map.key <- const$default$mapkeys[[genome]]
    }
    
    return(map.key)
}

# mapkeyOpt --------------------------------------------------------------------
#' Set/get \pkg{shmootl} \code{mapkey} option.
#' 
#' This function returns the current \pkg{shmootl} \code{mapkey}. If a new
#' \code{mapkey} is specified, the \pkg{shmootl} \code{mapkey} option is
#' updated with the given \code{mapkey} function.
#' 
#' @param value New \code{mapkey} function.
#' 
#' @return Current \code{mapkey} function.
#' 
#' @export
#' @rdname mapkeyOpt
mapkeyOpt <- function(value) {
    
    result <- getOption('shmootl.mapkey', default=mapkey())
    stopifnot( 'mapkey' %in% class(result) )
    
    if ( ! missing(value) ) {
        stopifnot( 'mapkey' %in% class(value) )
        options(shmootl.mapkey=value)
    }
    
    return(result)
}

# matchLoci --------------------------------------------------------------------
#' Get object loci matching those specified.
#' 
#' @param x Object containing map data.
#' @param loc Locus \code{mapframe} specifying map positions.
#' 
#' @return Locus mapframe of those loci specified that were found in the target
#' object.
#' 
#' @keywords internal
#' @rdname matchLoci
matchLoci <- function(x, loc) {
    x.loci <- pullLoci(x)
    indices <- matchLocusRowIndices(x.loci, loc)
    return(x.loci[indices, ])
}

# matchLocusRowIndices ---------------------------------------------------------
#' Find matching locus row indices.
#' 
#' @param x A \code{mapframe} or \code{scanone} object.
#' @param loc Locus \code{mapframe} specifying map positions.
#' 
#' @return Vector of row indices with respect to the target \code{mapframe},
#' for which a matching locus exists in the locus mapframe.
#' 
#' @keywords internal
#' @rdname matchLocusRowIndices
matchLocusRowIndices <- function(x, loc) {
    
    stopifnot( 'mapframe' %in% class(x) || 'scanone' %in% class(x) )
    stopifnot( 'mapframe' %in% class(loc) )
    
    map.unit <- getMapUnit(x)
    
    if ( is.na(map.unit) ) {
        stop("map unit not found")
    }
    
    if ( getMapUnit(loc) != map.unit ) {
        stop("map unit mismatch")
    }
    
    if ( nrow(x) == 0 ) {
        stop("cannot match locus row indices in empty mapframe")
    }
    
    if ( nrow(loc) == 0 ) {
        stop("locus mapframe must have at least one row")
    }
    
    # Get sequences and positions of individual map loci.
    x.seqs <- pullLocusSeq(x)
    x.pos <- pullLocusPos(x)
    
    # Get sequences and positions of specified loci.
    loc.seqs <- pullLocusSeq(loc)
    loc.pos <- pullLocusPos(loc)
    
    indices <- rep.int(NA_integer_, nrow(x))
    
    # Traverse both mapframes in sync, setting indices of matching loci.
    # NB: we can assume mapframe loci will match in order, if at all.
    i <- j <- 1
    while ( i <= nrow(x) ) {
        
        while ( ( loc.seqs[j] < x.seqs[i] || ( loc.seqs[j] == x.seqs[i] &&
            loc.pos[j] < x.pos[i] ) ) && j < nrow(loc) ) {
            j <- j + 1
        }
        
        while ( ( x.seqs[i] < loc.seqs[j] || ( x.seqs[i] == loc.seqs[j] &&
            x.pos[i] < loc.pos[j] ) ) && i < nrow(x) ) {
            i <- i + 1
        }
        
        # If closest locus is within numeric tolerance, set row index.
        if ( x.seqs[i] == loc.seqs[j] &&
            abs(x.pos[i] - loc.pos[j]) <= .Machine$double.eps^0.5 ) {
            indices[i] <- i
        }
        
        i <- i + 1
    }
    
    return(indices[ ! is.na(indices) ])
}

# matchSeqRowIndices -----------------------------------------------------------
#' Find row indices of matching sequences.
#' 
#' @param x A \code{mapframe} object.
#' @param sequences Sequences for which row indices should be returned.
#' 
#' @return Vector of row indices corresponding to the specified sequences. An
#' empty vector is returned if no sequences match the input \code{mapframe}.
#' 
#' @keywords internal
#' @rdname matchSeqRowIndices
matchSeqRowIndices <- function(x, sequences, simplify=FALSE) {
    
    stopifnot( is.data.frame(x) )
    
    if ( nrow(x) > 0 ) {
        
        norm.seqs <- normSeq(sequences)
        
        x.seqs <- pullLocusSeq(x)
        
        index.list <- vector('list', length(sequences))
        for ( i in getIndices(sequences) ) {
            index.list[[i]] <- which( x.seqs == norm.seqs[i] )
        }
        
        if (simplify) {
            indices <- sort( unlist(index.list) )
        } else {
            indices <- index.list
        }
        
    } else {
        
        indices <- integer()
    }
    
    return(indices)
}

# orderMap ---------------------------------------------------------------------
#' Put object in map order.
#' 
#' @param x Object containing map data.
#' 
#' @return Input object ordered by sequence, then map position.
#' 
#' @export
#' @family map utilities
#' @rdname orderMap
orderMap <- function(x) {
    UseMethod('orderMap', x)
}

# orderMap.data.frame ----------------------------------------------------------
#' @export
#' @rdname orderMap
orderMap.data.frame <- function(x) {
    
    if ( nrow(x) > 0 ) {
        x.seqs <- pullLocusSeq(x)
        x.pos <- pullLocusPos(x)
        x <- x[ order(rankSeq(x.seqs), x.pos), ]
    }
    
    return(x)
}

# orderMap.map -----------------------------------------------------------------
#' @export
#' @rdname orderMap
orderMap.map <- function(x) {
    
    # Sort map sequences.
    x <- subsetMap( x, sortSeq( names(x) ) )
    
    # Sort map positions.
    for ( i in getIndices(x) ) {
        x[[i]] <- sort.int(x[[i]])
    }
    
    return(x)
}

# pullLoci ---------------------------------------------------------------------
#' Pull loci from object.
#' 
#' @param x Object containing map data.
#' 
#' @return Locus mapframe of individual loci.
#' 
#' @keywords internal
#' @rdname pullLoci
pullLoci <- function(x) {
    
    x.seqs <- pullLocusSeq(x)
    
    x.pos <- pullLocusPos(x)
    
    loc <- data.frame(chr=x.seqs, pos=x.pos)
    
    class(loc) <- c('mapframe', 'data.frame')
    
    attr(loc, 'map.unit') <- getMapUnit(x)
    
    return(loc)
}

# pullLocusSeq -----------------------------------------------------------------
#' Pull sequence labels for individual loci.
#' 
#' @param x Object containing map data.
#' 
#' @return Character vector of sequence labels associated with individual loci.
#' 
#' @keywords internal
#' @rdname pullLocusSeq
pullLocusSeq <- function(x) {
    UseMethod('pullLocusSeq', x)
}

# pullLocusSeq.data.frame ------------------------------------------------------
#' @method pullLocusSeq data.frame
#' @rdname pullLocusSeq
pullLocusSeq.data.frame <- function(x) {
    return( as.character(x[, getSeqColIndex(x)]) )
}

# pullLocusSeq.list ------------------------------------------------------------
#' @method pullLocusSeq list
#' @rdname pullLocusSeq
pullLocusSeq.list <- function(x) {
    validateMap(x)
    return( pullLocusSeq.map(x) )
}

# pullLocusSeq.map -------------------------------------------------------------
#' @method pullLocusSeq map
#' @rdname pullLocusSeq
pullLocusSeq.map <- function(x) {
    return( unlist( lapply( getIndices(x), function(i)
        rep_len(names(x)[i], length(x[[i]])) ) ) )
}

# pullLocusIDs -----------------------------------------------------------------
#' Pull individual locus IDs.
#' 
#' @param x Object containing map data.
#' 
#' @return Character vector of individual locus IDs.
#' 
#' @keywords internal
#' @rdname pullLocusIDs
pullLocusIDs <- function(x) {
    UseMethod('pullLocusIDs', x)
}

# pullLocusIDs.data.frame ------------------------------------------------------
#' @method pullLocusIDs data.frame
#' @rdname pullLocusIDs
pullLocusIDs.data.frame <- function(x) {
    return( if ( nrow(x) > 0 && hasRownames(x) ) { rownames(x) } else { NULL } )
}

# pullLocusIDs.list ------------------------------------------------------------
#' @method pullLocusIDs list
#' @rdname pullLocusIDs
pullLocusIDs.list <- function(x) {
    validateMap(x)
    return( pullLocusIDs.map(x) )
}
  
# pullLocusIDs.map -------------------------------------------------------------
#' @method pullLocusIDs map
#' @rdname pullLocusIDs
pullLocusIDs.map <- function(x) {
    return( unlist( lapply(x, function(loci) names(loci)), use.names=FALSE ) )
}

# pullLocusPos -----------------------------------------------------------------
#' Pull individual locus positions.
#' 
#' @param x Object containing map data.
#' 
#' @return Numeric vector of individual locus positions.
#' 
#' @keywords internal
#' @rdname pullLocusPos
pullLocusPos <- function(x) {
    UseMethod('pullLocusPos', x)
}

# pullLocusPos.data.frame ------------------------------------------------------
#' @method pullLocusPos data.frame
#' @rdname pullLocusPos
pullLocusPos.data.frame <- function(x) {
    
    poscol.index <- getPosColIndex(x)
    
    if ( is.numeric(x[, poscol.index]) ) {
        locus.pos <- x[, poscol.index]
    } else {
        locus.pos <- setPosColDataMapUnit(x[, poscol.index], NULL)
    }
    
    return(locus.pos)
}

# pullLocusPos.list ------------------------------------------------------------
#' @method pullLocusPos list
#' @rdname pullLocusPos
pullLocusPos.list <- function(x) {
    validateMap(x)
    return( pullLocusPos.map(x) )
}

# pullLocusPos.map -------------------------------------------------------------
#' @method pullLocusPos map
#' @rdname pullLocusPos
pullLocusPos.map <- function(x) {
    return( unlist( lapply(x, function(loci) unname(loci) ) ) )
}

# pullLocusPos.mapframe --------------------------------------------------------
#' @method pullLocusPos mapframe
#' @rdname pullLocusPos
pullLocusPos.mapframe <- function(x) {
    return( x[, getPosColIndex(x)] )
}

# pullMap ----------------------------------------------------------------------
#' Pull map from object.
#' 
#' @param x Object containing map data.
#' 
#' @return An \pkg{R/qtl} \code{map} extracted from the input object.
#' 
#' @keywords internal
#' @rdname pullMap
pullMap <- function(x) {
    UseMethod('pullMap', x)
}

# pullMap.geno -----------------------------------------------------------------
#' @method pullMap geno
#' @rdname pullMap
pullMap.geno <- function(x) {
    return( as.map( lapply(x, function(obj)
        obj$map), map.unit='cM' ) )
}

# pushLoci ---------------------------------------------------------------------
#' Replace loci of object.
#' 
#' @param x Object containing map data.
#' @return loc Locus mapframe of individual loci.
#' 
#' @return Input object with individual loci replaced.
#' 
#' @keywords internal
#' @rdname pushLoci
pushLoci <- function(x, loc) {
    
    stopifnot( 'mapframe' %in% class(loc) )
    stopifnot( ncol(loc) == 2 )
    
    object.mapunit <- getMapUnit(x)
    loc.mapunit <- getMapUnit(loc)
    
    if ( ! is.na(object.mapunit) && loc.mapunit != object.mapunit ) {
        stop("map unit mismatch")
    }
    
    x <- pushLocusSeq(x, pullLocusSeq(loc))
    x <- pushLocusPos(x, pullLocusPos(loc))
    
    if ( is.na(object.mapunit) ) {
        attr(x, 'map.unit') <- loc.mapunit
    }
    
    return(x)
}

# pushLocusSeq -----------------------------------------------------------------
#' Replace sequence labels of individual loci.
#' 
#' @param x Object containing map data.
#' @param value Character vector of sequence labels for individual loci.
#' 
#' @return Input object with individual locus sequence labels replaced.
#' 
#' @keywords internal
#' @rdname pushLocusSeq
pushLocusSeq <- function(x, value) {
    stopifnot( all( isValidID(value) ) )
    UseMethod('pushLocusSeq', x)
}

# pushLocusSeq.data.frame ------------------------------------------------------
#' @method pushLocusSeq data.frame
#' @rdname pushLocusSeq
pushLocusSeq.data.frame <- function(x, value) {
    
    seqcol.index <- getSeqColIndex(x)
    
    # Get key info about object.
    x.runs <- rle( as.character(x[, seqcol.index]) ) # run-length encoding
    x.nmar <- x.runs$lengths         # number of markers by sequence
    x.totmar <- nrow(x)              # total number of markers
    
    # Get key info about value.
    value.runs <- rle(value)         # run-length encoding
    value.nmar <- value.runs$lengths # number of markers by sequence
    value.totmar <- length(value)    # total number of markers
    
    if ( value.totmar != x.totmar ) {
        stop("cannot push ", value.totmar, " locus sequences into ", x.totmar, " loci")
    }
    
    if ( ! identical(value.nmar, x.nmar) ) { # NB: checks equal number of sequences
        stop("cannot push locus sequences - sequence boundaries not preserved")
    }
    
    x[, seqcol.index] <- value
    
    return(x)
}

# pushLocusSeq.list ------------------------------------------------------------
#' @method pushLocusSeq list
#' @rdname pushLocusSeq
pushLocusSeq.list <- function(x, value) {
    validateMap(x)
    return( pushLocusSeq.map(x, value) )
}

# pushLocusSeq.map -------------------------------------------------------------
#' @method pushLocusSeq map
#' @rdname pushLocusSeq
pushLocusSeq.map <- function(x, value) {
    
    # Get key info about object.
    x.nmar <- lengths(x)             # number of markers by sequence
    x.totmar <- sum(x.nmar)          # total number of markers
    
    # Get key info about value.
    value.runs <- rle(value)         # run-length encoding
    value.nmar <- value.runs$lengths # number of markers by sequence
    value.totmar <- length(value)    # total number of markers
    
    if ( value.totmar != x.totmar ) {
        stop("cannot push ", value.totmar, " locus sequences into ", x.totmar, " loci")
    }
    
    if ( ! identical(value.nmar, x.nmar) ) { # NB: checks equal number of sequences
        stop("cannot push locus sequences - sequence boundaries not preserved")
    }
    
    names(x) <- value.runs$values
    
    return(x)
}

# pushLocusIDs -----------------------------------------------------------------
#' Replace individual locus IDs.
#' 
#' @param x Object containing map data.
#' @param value Character vector of individual locus IDs.
#' 
#' @return Input object with individual locus IDs replaced.
#' 
#' @keywords internal
#' @rdname pushLocusIDs
pushLocusIDs <- function(x, value) {
    stopifnot( all( isValidID(value) ) )
    stopif( anyDuplicated(value) )
    UseMethod('pushLocusIDs', x)
}

# pushLocusIDs.data.frame ------------------------------------------------------
#' @method pushLocusIDs data.frame
#' @rdname pushLocusIDs
pushLocusIDs.data.frame <- function(x, value) {
    
    # Get total number of markers in object.
    x.totmar <- nrow(x)

    # Get total number of markers in value.
    value.totmar <- length(value)
    
    if ( value.totmar != x.totmar ) {
        stop("cannot push ", value.totmar, " locus IDs into ", x.totmar, " loci")
    }
    
    rownames(x) <- value
    
    return(x)
}

# pushLocusIDs.list ------------------------------------------------------------
#' @method pushLocusIDs list
#' @rdname pushLocusIDs
pushLocusIDs.list <- function(x, value) {
    validateMap(x)
    return( pushLocusIDs.map(x, value) )
}

# pushLocusIDs.map -------------------------------------------------------------
#' @method pushLocusIDs map
#' @rdname pushLocusIDs
pushLocusIDs.map <- function(x, value) {
    
    # Get total number of markers in object.
    x.totmar <- sum( lengths(x) )
    
    # Get total number of markers in value.
    value.totmar <- length(value)
    
    if ( value.totmar != x.totmar ) {
        stop("cannot push ", value.totmar, " locus IDs into ", x.totmar, " loci")
    }
    
    map.seqs <- names(x)
    
    seq.index.list <- getRunIndexList( unlist( lapply( seq_along(x), function(i)
        rep_len(map.seqs[i], length(x[[i]])) ) ) )
    
    for ( map.seq in map.seqs ) {
        names(x[[map.seq]]) <- value[ seq.index.list[[map.seq]] ]
    }
    
    return(x)
}

# pushLocusPos -----------------------------------------------------------------
#' Replace individual locus positions.
#' 
#' @param x Object containing map data.
#' @param value Numeric vector of individual locus positions.
#' 
#' @return Input object with individual locus positions replaced.
#' 
#' @keywords internal
#' @rdname pushLocusPos
pushLocusPos <- function(x, value) {
    UseMethod('pushLocusPos', x)
}

# pushLocusPos.data.frame ------------------------------------------------------
#' @method pushLocusPos data.frame
#' @rdname pushLocusPos
pushLocusPos.data.frame <- function(x, value) {
    
    object.mapunit <- getMapUnit(x)
    
    value.mapunit <- getPosColDataMapUnit(value)
    
    if ( ! is.na(object.mapunit) && ! is.na(value.mapunit) &&
        value.mapunit != object.mapunit ) {
        stop("map unit mismatch")
    }
    
    value <- setPosColDataMapUnit(value, NULL)
    stopifnot( is.numeric(value) )
    
    # Get total number of markers in object.
    x.totmar <- nrow(x)
    
    # Get total number of markers in value.
    value.totmar <- length(value)
    
    if ( value.totmar != x.totmar ) {
        stop("cannot push ", value.totmar, " locus positions into ", x.totmar, " loci")
    }
    
    poscol.index <- getPosColIndex(x)
    
    x[, poscol.index] <- value
    
    stopifnot( inMapOrder(x) )
    
    if ( is.na(object.mapunit) ) {
        attr(x, 'map.unit') <- value.mapunit
    }
    
    return(x)
}

# pushLocusPos.list ------------------------------------------------------------
#' @method pushLocusPos list
#' @rdname pushLocusPos
pushLocusPos.list <- function(x, value) {
    validateMap(x)
    return( pushLocusPos.map(x, value) )
}

# pushLocusPos.map -------------------------------------------------------------
#' @method pushLocusPos map
#' @rdname pushLocusPos
pushLocusPos.map <- function(x, value) {
    
    object.mapunit <- getMapUnit(x)
    
    value.mapunit <- getPosColDataMapUnit(value)
    
    if ( ! is.na(object.mapunit) && ! is.na(value.mapunit) &&
        value.mapunit != object.mapunit ) {
        stop("map unit mismatch")
    }
    
    value <- setPosColDataMapUnit(value, NULL)
    stopifnot( is.numeric(value) )
    
    # Get total number of markers in object.
    x.totmar <- sum( lengths(x) )
    
    # Get total number of markers in value.
    value.totmar <- length(value)
    
    if ( value.totmar != x.totmar ) {
        stop("cannot push ", value.totmar, " locus positions into ", x.totmar, " loci")
    }
    
    map.seqs <- names(x)
    
    seq.index.list <- getRunIndexList( unlist( lapply( seq_along(x), function(i)
        rep_len(map.seqs[i], length(x[[i]])) ) ) )
    
    for ( map.seq in map.seqs ) {
        
        seq.pos <- value[ seq.index.list[[map.seq]] ]
        
        attributes(seq.pos) <- attributes(x[[map.seq]])
        
        x[[map.seq]] <- seq.pos
    }
    
    stopifnot( inMapOrder(x) )
    
    if ( is.na(object.mapunit) ) {
        attr(x, 'map.unit') <- value.mapunit
    }
    
    return(x)
}

# pushMap ----------------------------------------------------------------------
#' Push map into object.
#' 
#' @param x Object that can contain map data.
#' @param map An \pkg{R/qtl} \code{map} object.
#' 
#' @return Input object incorporating the given map data.
#' 
#' @keywords internal
#' @rdname pushMap
pushMap <- function(x, map) {
    UseMethod('pushMap', x)
}

# pushMap.geno -----------------------------------------------------------------
#' @method pushMap geno
#' @rdname pushMap
pushMap.geno <- function(x, map) {
    
    stopifnot( 'map' %in% class(map) )
    
    map.unit <- 'cM'
    
    if ( getMapUnit(map) != map.unit ) {
        stop("cross map positions must be in centiMorgans (e.g. '47 cM')")
    }
    
    norm.map.seqs <- normSeq( names(map) )
    
    geno.map <- as.map( lapply(x, function(obj) obj$map), map.unit=map.unit )
    
    geno.map.seqs <- names(geno.map)
    
    if ( length(norm.map.seqs) != length(geno.map.seqs) ||
         any( norm.map.seqs != geno.map.seqs ) ) {
        stop("map sequence mismatch")
    }
    
    for ( map.seq in norm.map.seqs ) {
        
        geno.seq.ids <- names(geno.map[[map.seq]])
        map.seq.ids <- names(map[[map.seq]])
        if ( length(map.seq.ids) != length(geno.seq.ids) ||
             any( map.seq.ids != geno.seq.ids ) ) {
            stop("marker mismatch")
        }
        
        seq.map <- map[[map.seq]]
        class(seq.map) <- 'numeric'
        
        x[[map.seq]]$map <- seq.map
    }
    
    return(x)
}

# setMapUnit -------------------------------------------------------------------
#' Set map unit.
#' 
#' Set the map unit of the object. For a \code{data.frame}, the position column
#' may be updated to reflect the new map unit.
#' 
#' @template mapunits
#' 
#' @param x Object containing map data.
#' @param map.unit Map unit.
#' 
#' @return Input object with the specified map unit.
#' 
#' @export
#' @family map utilities
#' @rdname setMapUnit
setMapUnit <- function(x, map.unit) {
    validateMapUnit(map.unit)
    UseMethod('setMapUnit', x)
}

# setMapUnit.list --------------------------------------------------------------
#' @export
#' @rdname setMapUnit
setMapUnit.list <- function(x, map.unit) {
    attr(x, 'map.unit') <- map.unit
    return(x)
}

# setMapUnit.map ---------------------------------------------------------------
#' @export
#' @rdname setMapUnit
setMapUnit.map <- function(x, map.unit) {
    
    existing.mapunit <- getMapUnit(x)
    
    if ( ! is.na(existing.mapunit) ) {
        x <- convertMapUnit(x, map.unit)
    } else {
        attr(x, 'map.unit') <- map.unit
    }
    
    return(x)
}

# setMapUnit.data.frame --------------------------------------------------------
#' @export
#' @rdname setMapUnit
setMapUnit.data.frame <- function(x, map.unit) {
    
    existing.mapunit <- getMapUnit(x)
    heading.mapunit <- getPosColNameMapUnit(x)
    positions.mapunit <- getPosColDataMapUnit(x)
    
    if ( ! is.na(heading.mapunit) ) {
        x <- setPosColNameMapUnit(x, NULL)
    }
    
    if ( ! is.na(positions.mapunit) ) {
        x <- setPosColNameMapUnit(x, NULL)
    }
    
    if ( ! is.na(existing.mapunit) ) {
        x <- convertMapUnit(x, map.unit)
    } else {
        attr(x, 'map.unit') <- map.unit
    }
    
    if ( ! is.na(heading.mapunit) ) {
        x <- setPosColNameMapUnit(x, map.unit)
    }
    
    if ( ! is.na(positions.mapunit) ) {
        x <- setPosColNameMapUnit(x, map.unit)
    }
    
    return(x)
}

# setMapUnit.mapframe ----------------------------------------------------------
#' @export
#' @rdname setMapUnit
setMapUnit.mapframe <- function(x, map.unit) {
    
    existing.mapunit <- getMapUnit(x)
    
    if ( ! is.na(existing.mapunit) ) {
        x <- convertMapUnit(x, map.unit)
    } else {
        attr(x, 'map.unit') <- map.unit
    }
    
    return(x)
}

# setPosColDataMapUnit ---------------------------------------------------------
#' Set map unit of position column data.
#' 
#' @template mapunits
#' 
#' @param x A \code{mapframe} or equivalent \code{data.frame},
#' or a vector of map positions.
#' @param map.unit Map unit. Set to \code{NULL} to remove map unit information
#' from position column data.
#' 
#' @return Input object with the specified map unit appended to position column
#' data.
#' 
#' @keywords internal
#' @rdname setPosColDataMapUnit
setPosColDataMapUnit <- function(x, map.unit) {
    
    # TODO: optimise this function.
    
    if ( ! is.null(map.unit) ) {
        validateMapUnit(map.unit)
    }
    
    UseMethod('setPosColDataMapUnit', x)
}

# setPosColDataMapUnit.data.frame ----------------------------------------------
#' @rdname setPosColDataMapUnit
setPosColDataMapUnit.data.frame <- function(x, map.unit) {
    
    existing.mapunit <- getMapUnit(x)
    
    if ( ! is.null(map.unit) && ! is.na(existing.mapunit) &&
        map.unit != existing.mapunit ) {
        stop("pos column data map-unit mismatch")
    }
    
    poscol.index <- getPosColIndex(x)
    
    x[, poscol.index] <- setPosColDataMapUnit(x[, poscol.index], map.unit)
    
    return(x)
}

# setPosColDataMapUnit.character -----------------------------------------------
#' @rdname setPosColDataMapUnit
setPosColDataMapUnit.character <- function(x, map.unit) {
    
    # Remove any existing map unit suffixes.
    existing.mapunit <- getPosColDataMapUnit(x)
    if ( ! is.na(existing.mapunit) ) {
        map.type <- const$known.map.types[existing.mapunit]
        u <- which( rownames(const$map.info[[map.type]]) == existing.mapunit )
        mapunit.pattern <- const$map.info[[map.type]]$pattern[u]
        res <- suppressWarnings( as.numeric(
            sub(mapunit.pattern, '', x, ignore.case=TRUE) ) )
    } else {
        res <- suppressWarnings( as.numeric(x) )
    }
    
    invalid <- unique(x[ is.na(res) ])
    if ( length(invalid) > 0 ) {
        stop("invalid map positions - '", invalid, "'")
    }
    
    # Append map unit suffixes.
    if ( length(res) > 0 && ! is.null(map.unit) ) {
        res <- paste(as.character(res), map.unit)
    }
    
    return(res)
}

# setPosColDataMapUnit.numeric -------------------------------------------------
#' @rdname setPosColDataMapUnit
setPosColDataMapUnit.numeric <- function(x, map.unit) {
    return( setPosColDataMapUnit( as.character(x), map.unit) )
}

# setPosColNameMapUnit ---------------------------------------------------------
#' Set map unit of position column heading.
#' 
#' @template mapunits
#' 
#' @param x A \code{mapframe} or equivalent \code{data.frame},
#' or a position column heading.
#' @param map.unit Map unit. Set to \code{NULL} to remove map unit information
#' from position column heading.
#' 
#' @return Input object with the specified map unit added to position column
#' heading.
#' 
#' @keywords internal
#' @rdname setPosColNameMapUnit
setPosColNameMapUnit <- function(x, map.unit) {
    
    if ( ! is.null(map.unit) ) {
        validateMapUnit(map.unit)
    }
    
    UseMethod('setPosColNameMapUnit', x)
}

# setPosColNameMapUnit.data.frame ----------------------------------------------
#' @rdname setPosColNameMapUnit
setPosColNameMapUnit.data.frame <- function(x, map.unit) {

    existing.mapunit <- getMapUnit(x)
    
    if ( ! is.null(map.unit) && ! is.na(existing.mapunit) &&
         map.unit != existing.mapunit ) {
        stop("pos column heading map unit mismatch")
    }
    
    poscol.index <- getPosColIndex(x)
    
    colnames(x)[poscol.index] <- setPosColNameMapUnit(
        colnames(x)[poscol.index], map.unit)
    
    return(x)
}

# setPosColNameMapUnit.character -----------------------------------------------
#' @rdname setPosColNameMapUnit
setPosColNameMapUnit.character <- function(x, map.unit) {
    
    stopifnot( length(x) == 1 )
    
    existing.mapunit <- getPosColNameMapUnit(x)
    
    if ( ! is.null(map.unit) ) {
        
        # If pos column has an existing map unit, replace it..
        if ( ! is.na(existing.mapunit) ) {
            x <- sub(existing.mapunit, map.unit, x)
        } else { # ..otherwise insert map unit in parentheses.
            x <- paste0(x, ' (', map.unit, ')')
        }
        
    } else if ( ! is.na(existing.mapunit) ) {
        
        # Compose pattern allowing for text before and after map
        # unit, all possible enclosing brackets, and any whitespace.
        pattern <- paste0('^(.*?)[[:space:]]*[[({<]?[[:space:]]*',
            existing.mapunit, '[[:space:]]*[])}>]?[[:space:]]*(.*?)$')
        
        # Match pattern, get any flanking text.
        m <- regexec(pattern, x)
        flank.text <- (regmatches(x, m))[[1]][2:3]
        flank.text <- flank.text[ flank.text != '' ]
        
        # If no flanking text, reset pos column heading.
        if ( length(flank.text) == 0 ) {
            parts <- 'pos'
        }
        
        # Merge parts of pos column heading.
        x <- paste(flank.text, collapse=' ')
    }
    
    return(x)
}

# setupDefaultMapkeys ----------------------------------------------------------
#' Setup default \code{mapkey} objects.
#' 
#' Default \code{mapkey} objects are generated from genome sequence info. Simple
#' component maps are created from sequence length info: a genetic map from the
#' \code{'maplengths'} column, and a physical map from the \code{'seqlengths'}
#' column. These are used to create a basic \code{mapkey} object for each genome.
#' 
#' @include const.R
#' @keywords internal
#' @rdname setupDefaultMapkeys
setupDefaultMapkeys <- function() {
    
    stopifnot( exists('const') )
    stopif( environmentIsLocked(const) )
    
    known.map.types <- unique(const$known.map.types)
    
    default.mapkeys <- list()
    
    # Create default mapkey for each genome.
    for ( genome in names(const$seqinfo) ) {
        
        seqinfo <- const$seqinfo[[genome]]
        
        map.seqs <- seqinfo$seqids
        
        cmaps <- vector('list', length(known.map.types))
        
        # Create component map for each known map type.
        for ( i in seq_along(known.map.types) ) {
            
            map.type <- known.map.types[i]
            
            # Get basic map unit for this map type.
            map.unit <- const$basic.map.unit[map.type]
            
            # Get key map info.
            if ( map.type == 'gmap' ) {
                makeLocusIDs <- makePseudomarkerIDs
                map.starts <- rep_len(0, length(map.seqs))
                map.ends <- seqinfo$maplengths
            } else if ( map.type == 'pmap' ) {
                makeLocusIDs <- makeDefaultMarkerIDs
                map.starts <- rep_len(1, length(map.seqs))
                map.ends <- seqinfo$seqlengths
            } else {
                stop("unknown map type - '", map.type, "'")
            }
            
            # Create map data-frame from sequence start- and end-points.
            map.table <- data.frame( chr=c(map.seqs, map.seqs),
                                     pos=c(map.starts, map.ends) )
            map.table <- map.table[ order(rankSeq(map.table$chr), map.table$pos), ]
            rownames(map.table) <- makeLocusIDs(map.table)
            
            # Get locus info.
            locus.seqs <- pullLocusSeq(map.table)
            locus.pos <- pullLocusPos(map.table)
            locus.ids <- pullLocusIDs(map.table)
            
            cmap <- list()
            
            # Set map data for each sequence.
            for ( map.seq in map.seqs ) {
                
                # Get indices for this sequence.
                indices <- which( locus.seqs == map.seq )
                
                # Get map positions for this sequence.
                seq.pos <- locus.pos[indices]
                
                # Set map position names from locus IDs.
                names(seq.pos) <- locus.ids[indices]
                
                # Set class of map position vector.
                class(seq.pos) <- 'A' # NB: assumes no 'X' chromosomes.
                
                # Set map positions for sequence.
                cmap[[map.seq]] <- seq.pos
            }
            
            # Set map unit of component map.
            attr(cmap, 'map.unit') <- map.unit
            
            # Set class of component map.
            class(cmap) <- 'map'
            
            cmaps[[i]] <- cmap
        }
        
        # Create default mapkey from component maps.
        default.mapkeys[[genome]] <- do.call(mapkey, cmaps)
    }
    
    const$default[['mapkeys']] <- default.mapkeys
    
    return( invisible() )
}

# subsetMap --------------------------------------------------------------------
#' Subset \pkg{shmootl} \code{map} object.
#' 
#' \pkg{R/qtl} implements a \code{subset.map} method, but this method strips
#' some attributes from the \code{map} object. This function transfers any
#' stripped attributes to the subsetted \code{map} object.
#' 
#' @param x A \code{map} object.
#' @param ... Arguments passed to \code{subset}.
#' 
#' @return Subsetted \code{map} object.
#' 
#' @keywords internal
#' @rdname subsetMap
subsetMap <- function(x, ...) {
    stopifnot( 'map' %in% class(x) )
    others <- otherattributes(x)
    x <- subset(x, subset)
    otherattributes(x) <- others
    return(x)
}

# validateMap ------------------------------------------------------------------
#' Validate \code{map} object.
#' 
#' @template map
#' 
#' @param x A \code{map} or equivalent \code{list}.
#' 
#' @return TRUE if object is a valid \code{map} or equivalent \code{list};
#' otherwise, returns first error.
#' 
#' @keywords internal
#' @rdname validateMap
validateMap <- function(x, ...) {
    UseMethod('validateMap', x)
}

# validateMap.list -------------------------------------------------------------
#' @rdname validateMap
validateMap.list <- function(x) {
    
    if ( length(x) < const$min.spm ) {
        stop("map-like object has too few sequences (min=", const$min.spm, ")")
    }
    
    map.seqs <- names(x)
    
    if ( is.null(map.seqs) ) {
        stop("map-like object has no sequence labels")
    }
    
    tryCatch({
        norm.map.seqs <- normSeq(map.seqs) # NB: validates sequence labels.
    }, error=function(e) {
        stop("map-like object has invalid sequence labels")
    })
    
    if ( anyDuplicated(map.seqs) ) {
        stop("map-like object has duplicate sequence labels")
    }
    
    if ( anyDuplicated(norm.map.seqs) ) {
        stop("map-like object has inconsistent sequence labels")
    }
    
    empty.seqs <- map.seqs[ is.na(x) ]
    if ( length(empty.seqs) > 0 ) {
        stop("map-like object has empty sequences - '",
            toString(empty.seqs), "'")
    }
    
    nonnum.seqs <- map.seqs[ ! sapply(x, is.numeric) ]
    if ( length(nonnum.seqs) > 0 ) {
        stop("map-like object has non-numeric sequence vectors - '",
            toString(nonnum.seqs), "'")
    }
    
    short.seqs <- map.seqs[ lengths(x) < const$min.lps ]
    if ( length(short.seqs) > 0 ) {
        stop("map-like object has too few loci in sequences - '",
            toString(short.seqs), "'")
    }
    
    absent.ids <- map.seqs[ ! sapply(x, hasNames) ]
    if ( length(absent.ids) > 0 ) {
        stop("locus IDs not found for sequences - '",
            toString(absent.ids), "'")
    }
    
    x.ids <- pullLocusIDs(x)
    
    if ( ! all( isValidID(x.ids) ) ) {
        stop("map-like object has invalid locus IDs")
    }
    
    if ( anyDuplicated(x.ids) ) {
        stop("map-like object has duplicate locus IDs")
    }
    
    tryCatch({
        x.pos <- pullLocusPos(x) # NB: validates map positions.
    }, error=function(e) {
        stop("map-like object has invalid map positions")
    })
    
    validateMapUnit(x)
    
    if ( isPhysicalMap(x) ) {
        subzero <- unique( x.pos[ x.pos < 0 ] )
        if ( length(subzero) > 0 ) {
            stop("map-like object has invalid physical map positions - '",
                toString(subzero), "'")
        }
    }
    
    return(TRUE)
}

# validateMap.map --------------------------------------------------------------
#' @rdname validateMap
validateMap.map <- function(x, package=c('qtl', 'shmootl')) {
    
    package <- match.arg(package)
    
    if ( length(x) < const$min.spm ) {
        stop("map has too few sequences (min=", const$min.spm, ")")
    }
    
    map.seqs <- names(x)
    
    if ( is.null(map.seqs) ) {
        stop("map has no sequence labels")
    }
    
    if ( package == 'qtl' ) {
    
        tryCatch({
            norm.map.seqs <- normSeq(map.seqs) # NB: validates sequence labels.
        }, error=function(e) {
            stop("qtl map has invalid sequence labels")
        })
    
    } else {
        
        invalid <- map.seqs[ ! isNormSeq(map.seqs) ]
        if ( length(invalid) > 0 ) {
            stop("shmootl map has invalid sequence labels - '", toString(invalid), "'")
        }
    }
    
    dup.seqs <- map.seqs[ duplicated(map.seqs) ]
    if ( length(dup.seqs) > 0 ) {
        stop("map has duplicate sequence labels - '", toString(dup.seqs), "'")
    }
    
    empty.seqs <- map.seqs[ is.na(x) ]
    if ( length(empty.seqs) > 0 ) {
        stop("map has empty sequences - '", toString(empty.seqs), "'")
    }
    
    nonnum.seqs <- map.seqs[ ! sapply(x, is.numeric) ]
    if ( length(nonnum.seqs) > 0 ) {
        stop("map has non-numeric sequence vectors - '", toString(nonnum.seqs), "'")
    }
    
    short.seqs <- map.seqs[ lengths(x) < const$min.lps ]
    if ( length(short.seqs) > 0 ) {
        stop("map has too few loci in sequences - '", toString(short.seqs), "'")
    }
    
    absent.ids <- map.seqs[ ! sapply(x, hasNames) ]
    if ( length(absent.ids) > 0 ) {
        stop("locus IDs not found for sequences - '", toString(absent.ids), "'")
    }
    
    x.ids <- pullLocusIDs(x)
    
    if ( ! all( isValidID(x.ids) ) ) {
        stop("map has invalid locus IDs")
    }
    
    if ( anyDuplicated(x.ids) ) {
        stop("map has duplicate locus IDs")
    }
    
    x.pos <- pullLocusPos(x)
    
    if ( anyNA(x.pos) ) {
        stop("map has invalid map positions")
    }
    
    validateMapUnit(x)
    
    if ( isPhysicalMap(x) ) {
        subzero <- unique( x.pos[ x.pos < 0 ] )
        if ( length(subzero) > 0 ) {
            stop("map has invalid physical map positions - '",
                 toString(subzero), "'")
        }
    }
    
    if ( ! inMapOrder(x) ) {
        stop("map must be ordered by sequence and map position")
    }
    
    if ( package == 'qtl' ) {
        
        if ( ! all( sapply(x, class) %in% c('A', 'X') ) ) {
            stop("map sequences must be of class 'A' or 'X'")
        }
    
    } else {
        
        if ( ! all( sapply(x, class) == 'A' ) ) { # NB: assumes no 'X' chromosomes.
            stop("shmootl map sequences must be of class 'A'")
        }
        
        if ( is.na( getMapUnit(x) ) ) {
            stop("shmootl map must have map unit info")
        }
    }
    
    return(TRUE)
}

# validateMapframe -------------------------------------------------------------
#' Validate \code{mapframe} object.
#' 
#' @template mapframe
#' 
#' @param x A \code{mapframe} or equivalent \code{data.frame}.
#' 
#' @return TRUE if object is a valid \code{mapframe} or equivalent
#' \code{data.frame}; otherwise, returns first error.
#' 
#' @keywords internal
#' @rdname validateMapframe
validateMapframe <- function(x, ...) {
    UseMethod('validateMapframe', x)
}

# validateMapframe.data.frame --------------------------------------------------
#' @rdname validateMapframe
validateMapframe.data.frame <- function(x) {
    
    seqcol.index <- getSeqColIndex(x) # NB: confirms object has sequence column.
    poscol.index <- getPosColIndex(x) # NB: confirms object has position column.

    if ( nrow(x) > 0 ) {
        
        x.ids <- pullLocusIDs(x)

        if ( ! is.null(x.ids) ) {
            
            if ( ! all( isValidID(x.ids) ) ) {
                stop("mapframe-like object has invalid locus IDs")
            }

            if ( anyDuplicated(x.ids) ) {
                stop("mapframe-like object has duplicate locus IDs")
            }
        }

        tryCatch({
            norm.map.seqs <- normSeq( unique( pullLocusSeq(x) ) ) # NB: validates sequence labels.
        }, error=function(e) {
            stop("mapframe-like object has invalid sequence labels")
        })

        if ( anyDuplicated(norm.map.seqs) ) {
            stop("mapframe-like object has inconsistent sequence labels")
        }

        tryCatch({
            x.pos <- pullLocusPos(x) # NB: validates map positions.
        }, error=function(e) {
            stop("mapframe-like object has invalid map positions")
        })

        validateMapUnit(x)

        if ( isPhysicalMapframe(x) ) {
            subzero <- unique( x.pos[ x.pos < 0 ] )
            if ( length(subzero) > 0 ) {
                stop("mapframe-like object has invalid physical map positions - '",
                     toString(subzero), "'")
            }
        }
    }
    
    return(TRUE)
}

# validateMapframe.mapframe ----------------------------------------------------
#' @rdname validateMapframe
validateMapframe.mapframe <- function(x) {
    
    mapframe.colnames <- colnames(x)
    
    if ( mapframe.colnames[1] != 'chr' ) {
        stop("name of first mapframe column must be 'chr', not '",
            mapframe.colnames[1], "'")
    }
    
    if ( mapframe.colnames[2] != 'pos' ) {
        stop("name of second mapframe column must be 'pos', not '",
             mapframe.colnames[2], "'")
    }
    
    if ( nrow(x) > 0 ) {
        
        x.ids <- pullLocusIDs(x)
        
        if ( ! is.null(x.ids) ) {
            
            if ( ! all( isValidID(x.ids) ) ) {
                stop("mapframe has invalid locus IDs")
            }
            
            if ( anyDuplicated(x.ids) ) {
                stop("mapframe has duplicate locus IDs")
            }
        }
        
        map.seqs <- unique( pullLocusSeq(x) )
        invalid <- map.seqs[ ! isNormSeq(map.seqs) ]
        if ( length(invalid) > 0 ) {
            stop("mapframe has invalid sequence labels - '",
                toString(invalid), "'")
        }
        
        x.pos <- pullLocusPos(x)
        
        if ( ! is.numeric(x.pos) ) {
            stop("mapframe has non-numeric map positions")
        }
        
        if ( anyNA(x.pos) ) {
            stop("mapframe has invalid map positions")
        }
        
        validateMapUnit(x)
        
        if ( isPhysicalMapframe(x) ) {
            subzero <- unique( x.pos[ x.pos < 0 ] )
            if ( length(subzero) > 0 ) {
                stop("mapframe has invalid physical map positions - '",
                     toString(subzero), "'")
            }
        }
        
        if ( ! inMapOrder(x) ) {
            stop("mapframe must be ordered by sequence and position")
        }
    }
   
    if ( is.na( getMapUnit(x) ) ) {
        stop("mapframe must have map unit info")
    }
    
    return(TRUE)
}

# validateMapUnit --------------------------------------------------------------
#' Validate map unit.
#' 
#' @template mapunits
#' 
#' @param x Map unit, or object with a \code{'map.unit'} attribute.
#' @param map.type Map type for which map unit should be validated. If not
#' specified, map unit is validated for all map types.
#' 
#' @return TRUE if the specified map unit is a valid map unit for the given
#' map type(s); otherwise, returns error.
#' 
#' @keywords internal
#' @rdname validateMapUnit
validateMapUnit <- function(x, map.type=NULL) {
    UseMethod('validateMapUnit', x)
}

# validateMapUnit.character ----------------------------------------------------
#' @rdname validateMapUnit
validateMapUnit.character <- function(x, map.type=NULL) {
    
    if ( ! is.null(map.type) ) {
        
        if ( ! map.type %in% names(const$map.info) ) {
            stop("unknown map type - '", map.type, "'")
        }
        
        known.map.units <- rownames(const$map.info[[map.type]])
        
    } else {
        
        known.map.units <- const$known.map.units
    }

    if ( length(x) != 1 || is.na(x) || ! x %in% known.map.units ) {
        stop("invalid map unit - '", toString(x), "'")
    }
    
    return(TRUE)
}

# validateMapUnit.default ------------------------------------------------------
#' @rdname validateMapUnit
validateMapUnit.default <- function(x, map.type=NULL) {
    
    map.unit <- getMapUnit(x)
    
    if ( length(map.unit) == 1 && is.na(map.unit) ) {
        stop("no map unit found")
    }
    
    validateMapUnit(map.unit, map.type=map.type)
    
    return(TRUE)
}

# validateGeneticMapUnit -------------------------------------------------------
#' Validate genetic map unit.
#' 
#' @template mapunits
#' 
#' @param x Map unit, or object with a \code{'map.unit'} attribute.
#' 
#' @return TRUE if map unit is a valid genetic map unit;
#' otherwise, returns error.
#' 
#' @keywords internal
#' @rdname validateGeneticMapUnit
validateGeneticMapUnit <- function(x) {
    return( validateMapUnit(x, map.type='gmap') )
}

# validatePhysicalMapUnit ------------------------------------------------------
#' Validate physical map unit.
#' 
#' @template mapunits
#' 
#' @param x Map unit, or object with a \code{'map.unit'} attribute.
#' 
#' @return TRUE if map unit is a valid physical map unit;
#' otherwise, returns error.
#' 
#' @keywords internal
#' @rdname validatePhysicalMapUnit
validatePhysicalMapUnit <- function(x) {
    return( validateMapUnit(x, map.type='pmap') )
}

# End of map.R #################################################################