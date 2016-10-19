# Start of util.R ##############################################################

# allKwargs --------------------------------------------------------------------
#' Test if ellipsis arguments are all keyword arguments.
#' 
#' @param ... Ellipsis arguments.
#'     
#' @return \code{TRUE} if all ellipsis arguments are keyword arguments;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname allKwargs
allKwargs <- function(...) {
    args <- list(...)
    return( length(args) == 0 || hasNames(args) )
}

# anyKwargs --------------------------------------------------------------------
#' Test if any ellipsis arguments are keyword arguments.
#' 
#' @param ... Ellipsis arguments.
#' 
#' @return \code{TRUE} if any ellipsis arguments are keyword arguments;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname anyKwargs
anyKwargs <- function(...) {
    args <- list(...)
    return( length(args) > 0 && ! is.null( names(args) ) &&
        any( ! is.na(names(args)) & names(args) != '' ) )
}

# allNA ------------------------------------------------------------------------
#' Test if all elements are \code{NA} values.
#' 
#' @param x Test vector.
#'     
#' @return \code{TRUE} if vector is of length zero or
#' contains only \code{NA} values; \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname allNA
allNA <- function(x) {
    stopifnot( is.vector(x) )
    return( all( is.na(x) ) )
}

# allWhite ---------------------------------------------------------------------
#' Test if vector is all whitespace.
#' 
#' @param x Character vector.
#' 
#' @return \code{TRUE} if character vector is of length zero or
#' contains only whitespace characters; \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname allWhite
allWhite <- function(x) {
    stopifnot( is.character(x) )
    return( all( grepl( "^[[:space:]]*$", x)  ) )
}

# bstripBlankRows --------------------------------------------------------------
#' Strip blank rows from bottom of \code{data.frame}.
#' 
#' @param x A \code{data.frame} with columns of type \code{character}.
#' 
#' @return Input object in which bottommost blank rows (i.e. rows
#' containing no non-whitespace characters) have been stripped.
#' 
#' @keywords internal
#' @rdname bstripBlankRows
bstripBlankRows <- function(x) {
    stopifnot( is.data.frame(x) )
    stopifnot( all( sapply(x, class) == 'character' ) )
    while( allWhite( as.character( x[nrow(x), ]) ) ) {
        x <- x[-nrow(x), , drop=FALSE]
    }
    return(x)
}

# clamp ------------------------------------------------------------------------
#' Clamp numbers within a range.
#' 
#' @param n Numeric vector.
#' @param interval Numeric vector containing the minimum and maximum values of
#' the range, respectively. 
#'     
#' @return Input vector with values clamped within the specified range.
#' 
#' @keywords internal
#' @rdname clamp
clamp <- function(n, interval) {
    
    stopifnot( is.numeric(n) )
    stopifnot( is.numeric(interval) )
    stopifnot( length(interval) == 2 )
    stopif( anyNA(interval) )
    stopif( any( is.nan(interval) ) )
    stopifnot( diff(interval) >= 0 )
    
    n[ ! is.na(n) & n < interval[1] ] <- interval[1]
    n[ ! is.na(n) & n > interval[2] ] <- interval[2]
    
    return(n)
}

# coerceDataFrame --------------------------------------------------------------
#' Coerce \code{data.frame} columns to the specified classes.
#' 
#' @param x A \code{data.frame}.
#' @param classes Classes to which columns of \code{data.frame} are to be 
#' coerced. This should be a character vector containing class names, with 
#' the same number of elements as there are columns in the \code{data.frame}.
#'     
#' @return Coerced \code{data.frame}.
#' 
#' @keywords internal
#' @rdname coerceDataFrame
coerceDataFrame <- function(x, classes) {
    
    stopifnot( is.data.frame(x) )
    stopifnot( nrow(x) > 0 )
    stopifnot( ncol(x) > 0 )
    stopifnot( is.character(classes) )
    stopifnot( length(classes) == ncol(x) )
    
    coercions <- lapply(classes, getCoercionFromClassS3)
    
    stopif( any( is.null(coercions) ) )
    
    for ( i in getColIndices(x) ) {
        x[[i]] <- coercions[[i]](x[[i]])
    }
    
    return(x)
}

# deleteColumn -----------------------------------------------------------------
#' Delete column from object.
#' 
#' @param x A \code{data.frame} or \code{matrix}.
#' @param col.index Column index of column to remove.
#' @param col.name Column name of column to remove.
#'     
#' @return Input object with column removed.
#' 
#' @keywords internal
#' @rdname deleteColumn
deleteColumn <- function(x, col.index=NULL, col.name=NULL) {
    
    stopifnot( is.data.frame(x) || is.matrix(x) )
    
    if ( ! xor( is.null(col.index), is.null(col.name) ) ) {
        stop("deleteColumn takes either a column index or column name, but not both")
    }
    
    if ( ! is.null(col.index) ) {
        stopifnot( isSingleWholeNumber(col.index) )
        stopifnot( inRange( col.index, c( 1, ncol(x) ) ) )
    } else {
        stopifnot( isSingleString(col.name) )
        stopif( is.null( colnames(x) ) )
        stopif( anyDuplicated( colnames(x) ) )
        col.index <- which( colnames(x) == col.name )
    } 
    
    others <- otherattributes(x)
    
    # Delete column.
    x <- x[, -col.index, drop=FALSE]
    
    otherattributes(x) <- others
    
    return(x)
}

# dispatchFromClassS3 ----------------------------------------------------------
#' Dispatch method with respect to the given class vector.
#' 
#' @param generic Generic function name.
#' @param class.vector Class vector of an R object.
#' @param package Package from which available methods are taken.
#'     
#' @return Method suited to the specified class.
#' 
#' @importFrom utils lsf.str
#' @keywords internal
#' @rdname dispatchFromClassS3
dispatchFromClassS3 <- function(generic, class.vector, package) {
    
    stopifnot( isSingleString(generic) )
    stopifnot( is.character(class.vector) )
    stopifnot( length(class.vector) > 0 )
    stopifnot( isSingleString(package) )
    
    pkg <- paste0('package:', package)
    pattern <- paste0('^', generic, '[.]([^[:space:]]+)')
    x <- utils::lsf.str(pkg, pattern=pattern)
    
    m <- regexec(pattern, x)
    matches <- regmatches(x, m)

    functions <- sapply(matches, getElement, 1)
    suffixes <- sapply(matches, getElement, 2)
    
    default.index <- which( suffixes == 'default' )
    class.indices <- which( suffixes != 'default' )
    
    func.name <- NULL
    
    if ( length(class.indices) > 0 ) {
        
        classes <- suffixes[ class.indices ]
        
        class.match <- match(classes, class.vector)
        
        if ( length(class.match) > 0 && ! allNA(class.match) ) {
            
            first.match <- min(class.match, na.rm = TRUE)
            cls <- classes[ ! is.na(class.match) & class.match == first.match ]
            func.name <- functions[ suffixes == cls ]
        }
    }
    
    if ( is.null(func.name) ) {
        
        if ( length(default.index) > 0 ) {
            
            func.name <- functions[default.index]
            
        } else {
            
            stop("no matching class, no default")
        }
    }
    
    return( get(func.name) )
}

# emptyArgs --------------------------------------------------------------------
#' Test if ellipsis arguments are empty.
#' 
#' @param ... Ellipsis arguments.
#' 
#' @return \code{TRUE} if there are no ellipsis arguments;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname emptyArgs
emptyArgs <- function(...) {
    return( length( list(...) ) == 0 )
}

# getCoercionFromClassS3 -------------------------------------------------------
#' Get coercion function for the given class.
#' 
#' @param class.vector Class vector of an R object.
#'     
#' @return Coercion function.
#' 
#' @importFrom utils as.roman
#' @keywords internal
#' @rdname getCoercionFromClassS3
getCoercionFromClassS3 <- function(class.vector) {
    
    stopifnot( is.character(class.vector) )
    stopifnot( length(class.vector) > 0 )
    
    known.coercions <- list(
        character=as.character,
        double=as.double,
        factor=as.factor,
        integer=as.integer,
        logical=as.logical,
        numeric=as.numeric,
        raw=as.raw,
        roman=utils::as.roman
    )
    
    coercion <- NULL
    
    for ( class.name in class.vector ) {
        
        if ( class.name %in% names(known.coercions) ) {
            coercion <- known.coercions[[class.name]]
            break
        }
    }
    
    return(coercion)
}

# getColIndices ----------------------------------------------------------------
#' Get column indices of object.
#' 
#' @param x An object with columns.
#'     
#' @return Integer vector of all column indices. Returns an empty integer vector
#' if the object has zero columns.
#'    
#' @keywords internal
#' @rdname getColIndices
getColIndices <- function(x) {
    num.cols <- ncol(x)
    stopifnot( isSingleNonNegativeWholeNumber(num.cols) )
    return( if ( num.cols > 0 ) { 1:num.cols } else { integer() } )
}

# getIdColIndex ----------------------------------------------------------------
#' Get sample ID column index.
#' 
#' Get the index of the sample ID column in the given object. For example,
#' in a \code{cross} object, this is the column of the phenotype
#' \code{data.frame} with the heading \code{'ID'} (case-insensitive). It can
#' also be used to get the ID column index directly from a \code{data.frame},
#' in which case the ID column is that which has heading \code{'ID'}
#' (also case-insensitive).
#' 
#' @param x An object that may contain a sample ID column.
#' 
#' @return Sample ID column index.
#' 
#' @export
#' @family cross object functions
#' @rdname getIdColIndex
getIdColIndex <- function(x) {
    UseMethod('getIdColIndex', x)
}

# getIdColIndex.cross ----------------------------------------------------------
#' @export
#' @rdname getIdColIndex
getIdColIndex.cross <- function(x) {
    
    id.col <- which( tolower( colnames(x$pheno) ) == 'id' )
    
    if ( length(id.col) == 0 ) {
        id.col <- NULL
    } else if ( length(id.col) > 1 ) {
        stop("multiple ID columns found")
    }
    
    return(id.col)
}

# getIdColIndex.data.frame -----------------------------------------------------
#' @export
#' @method getIdColIndex data.frame
#' @rdname getIdColIndex
getIdColIndex.data.frame <- function(x) {
    
    id.col <- which( tolower( colnames(x) ) == 'id' )
    
    if ( length(id.col) == 0 ) {
        id.col <- NULL
    } else if ( length(id.col) > 1 ) {
        stop("multiple ID columns found")
    }
    
    return(id.col)
}

# getMissingValueFromClassS3 ---------------------------------------------------
#' Get missing value for the given class.
#' 
#' @param class.vector Class vector of an R object.
#'     
#' @return Missing value.
#' 
#' @keywords internal
#' @rdname getMissingValueFromClassS3
getMissingValueFromClassS3 <- function(class.vector) {
    
    stopifnot( is.character(class.vector) )
    stopifnot( length(class.vector) > 0 )
    
    known.missing <- list(
        character = NA_character_,
        double    = NA_real_,
        factor    = NA_character_,
        integer   = NA_integer_,
        logical   = NA,
        numeric   = NA_real_
    )
    
    missing.value <- NULL
    
    for ( class.name in class.vector ) {
        
        if ( class.name %in% names(known.missing) ) {
            missing.value <- known.missing[[class.name]]
            break
        }
    }
    
    return(missing.value)
}

# getRowIndices ----------------------------------------------------------------
#' Get row indices of object.
#' 
#' @param x An object with rows.
#'     
#' @return Integer vector of all row indices. Returns an empty integer vector if
#' the object has zero rows.
#'    
#' @keywords internal
#' @rdname getRowIndices
getRowIndices <- function(x) {
    num.rows <- nrow(x)
    stopifnot( isSingleNonNegativeWholeNumber(num.rows) )
    return( if ( num.rows > 0 ) { 1:num.rows } else { integer() } )
}

# getRunIndexList --------------------------------------------------------------
#' Get index list of successive runs in a vector.
#' 
#' @param x A vector.
#'     
#' @return List of integer vectors, each containing indices for a run of
#' repeated values in the input vector. Each list element takes its name
#' from the corresponding repeated value. Returns an empty list if the
#' input vector is of length zero.
#'    
#' @keywords internal
#' @rdname getRunIndexList
getRunIndexList <- function(x, na.rm=FALSE) {
    
    stopifnot( is.vector(x) )
    stopifnot( isBOOL(na.rm) )
    
    if ( length(x) > 0 ) {
        
        # Get run-length encoding of vector.
        runs <- rle(x)
        
        # Set run names from RLE values.
        run.names <- runs$values
        
        # Get number of runs in RLE.
        num.runs <- unique( lengths(runs) )
        
        # Get last index of each run.
        J <- cumsum(runs$lengths)
        
        # Get first index of each run.
        if ( num.runs > 1 ) {
            I <- c( 1, sapply(J[1:(length(J)-1)], function(j) j + 1) )
        } else {
            I <- 1
        }
        
        # Remove NA values, if specified.
        if (na.rm) {
            mask <- ! is.na(runs$values)
            run.names <- run.names[mask]
            I <- I[mask]
            J <- J[mask]
        }
        
        # Set index list from run index ranges.
        index.list <- mapply(function(i, j) i:j, I, J, SIMPLIFY=FALSE)
        
        # Set names of index list from run values.
        names(index.list) <- run.names
        
    } else {
        
        index.list <- list()
    }
    
    return(index.list)
}

# getRunIndices ----------------------------------------------------------------
#' Get indices of successive runs in a run-length encoding.
#' 
#' @param x An \code{rle} object.
#' 
#' @return Integer vector of all run indices in the given run-length encoding,
#' which can be used to index into the \code{lengths} and \code{values} of the
#' given \code{rle} object. Returns an empty integer vector if the run-length
#' encoding has zero runs.
#' 
#' @keywords internal
#' @rdname getRunIndices
getRunIndices <- function(x) {
    stopifnot( 'rle' %in% class(x) )
    num.runs <- union( length(x$lengths), length(x$values) )
    stopifnot( isSingleNonNegativeWholeNumber(num.runs) )
    return( if ( num.runs > 0 ) { 1:num.runs } else { integer() } )
}

# getSeqinfo -------------------------------------------------------------------
#' Get \pkg{GenomeInfoDb} \code{Seqinfo} object.
#' 
#' @param genome Genome for which a \pkg{GenomeInfoDb} \code{Seqinfo} object
#' should be returned. If this is not specified, a \code{Seqinfo} object is
#' returned for the current genome.
#' 
#' @return A \code{Seqinfo} object for the given genome.
#'
#' @importFrom GenomeInfoDb Seqinfo
#' @keywords internal
#' @rdname getSeqinfo
getSeqinfo <- function(genome=NULL) {
    
    seqtab <- getSeqTable(genome)
    
    seqinfo <- GenomeInfoDb::Seqinfo(
        seqnames   = seqtab$seqids,
        seqlengths = seqtab$seqlengths,
        isCircular = seqtab$isCircular,
        genome     = seqtab$genome
    )
    
    return(seqinfo)
}

# getSeqTable ------------------------------------------------------------------
#' Get genome sequence info table.
#' 
#' @param genome Genome for which sequence info should be returned. If this
#' is not specified, sequence info is returned for the current genome.
#' 
#' @return A \code{data.frame} containing sequence info for the given genome.
#' This can be used to create, but is distinct from, a \pkg{GenomeInfoDb}
#' \code{Seqinfo} object.
#' 
#' @keywords internal
#' @rdname getSeqTable
getSeqTable <- function(genome=NULL) {
    prev.genome <- genomeOpt(genome)
    on.exit( genomeOpt(prev.genome) )
    seqtab <- const$seqtab[[ genomeOpt() ]]
    return(seqtab)
}

# getSpecialAttributeNames -----------------------------------------------------
#' Get special attributes for the given object.
#' 
#' @param x R object.
#'     
#' @return Vector of special attribute names.
#' 
#' @keywords internal
#' @rdname getSpecialAttributeNames
getSpecialAttributeNames <- function(x) {
    
    default.specials <- const$special.attributes[['default']]
    
    attrset.names <- names(const$special.attributes)
    
    class.specials <- NULL
    for ( class.name in class(x) ) {
        if ( class.name %in% attrset.names ) {
            class.specials <- const$special.attributes[[class.name]]
            break
        }
    }
    
    return( union(class.specials, default.specials) )
}

# hasNames ---------------------------------------------------------------------
#' Test if object has names.
#' 
#' @param x Test object.
#'     
#' @return \code{TRUE} if object has nonempty names with no \code{NA} values;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname hasNames
hasNames <- function(x) {
    return( ! is.null( names(x) ) && all( ! is.na(names(x)) & names(x) != '' ) )
}

# hasRownames ------------------------------------------------------------------
#' Test if object has rownames.
#' 
#' @param x A \code{data.frame} or \code{matrix}.
#'     
#' @return \code{TRUE} if input object has non-default row names;
#' \code{FALSE} otherwise.
#'    
#' @keywords internal
#' @rdname hasRownames
hasRownames <- function(x) {
    
    rowname.status <- FALSE
    
    if ( is.data.frame(x) ) {
        row.names <- attr(x, 'row.names')
        rowname.status <- ! is.null(row.names) && ! all( row.names == getRowIndices(x) )
    } else if ( ! is.null( rownames(x) ) ) {
        rowname.status <- TRUE
    } 
    
    return(rowname.status) 
}

# inferStepSize ----------------------------------------------------------------
#' Infer step size from step values.
#' 
#' @param steps Numeric vector of steps (i.e. differences between consecutive 
#' values). 
#' @param tol Tolerance for step equality.
#'     
#' @return Inferred step size. Returns \code{NULL} if step size cannot be 
#' inferred.
#' 
#' @keywords internal
#' @rdname inferStepSize
inferStepSize <- function(steps, tol=.Machine$double.eps^0.5) {
    
    stopifnot( is.numeric(steps) )
    stopifnot( all( steps >= 0 ) )
    stopifnot( isSingleNonNegativeNumber(tol) )
    
    # Get frequency table of steps.
    step.freq <- table(steps)
    
    # Get numeric step sizes.
    step.sizes <- as.numeric( names(step.freq) )
    
    # Get differences between step sizes.
    step.diff <- diff(step.sizes)
    
    # Group step-size values that are very similar.
    size.groups <- vector('list')
    i <- 1L
    while ( i <= length(step.freq) ) {
        
        j <- i
        
        while ( j < length(step.freq) && step.diff[j] < tol ) {
            j <- j + 1L
        }
        
        size.groups <- append( size.groups, list(i:j) )
        
        i <- j + 1L
    }
    
    # Merge similar step values.
    merged.freq <- integer( length=length(size.groups) )
    for ( i in seq_along(size.groups) ) {
        g <- unlist(size.groups[i])
        merged.freq[i] <- sum(step.freq[g])
        names(merged.freq)[i] <- names( sort(step.freq[g], decreasing=TRUE) )[1]
    }
    
    # Sort steps by decreasing frequency.
    sorted.freq <- sort(merged.freq, decreasing=TRUE)
    
    # Verify that most frequent step is in majority.
    if( length(sorted.freq) > 1 && sorted.freq[1] < sum(sorted.freq[2:length(sorted.freq)]) ) {
        return(NULL)
    }
    
    # Get most frequent step.
    step.size <- as.numeric( names(sorted.freq)[1] )
    
    return(step.size)
}

# inRange ----------------------------------------------------------------------
#' Test if numbers lie within a range.
#' 
#' @param n Numeric vector.
#' @param interval Numeric vector containing the minimum and maximum values of
#' the range, respectively. 
#'     
#' @return Logical vector indicating which elements of the input numeric vector
#' lie within the specified range.
#' 
#' @keywords internal
#' @rdname inRange
inRange <- function(n, interval) {
    
    stopifnot( is.numeric(n) )
    stopifnot( is.numeric(interval) )
    stopifnot( length(interval) == 2 )
    stopif( anyNA(interval) )
    stopif( any( is.nan(interval) ) )
    stopifnot( diff(interval) >= 0 )
    
    return( findInterval(n, interval, rightmost.closed=TRUE) == 1 )
}

# insertColumn -----------------------------------------------------------------
#' Insert column in an object.
#' 
#' @param x A \code{data.frame} or \code{matrix}.
#' @param col.index Column index of inserted column.
#' @param col.name Column name of inserted column. If a column name is 
#' specified, but the object has no existing column names, default names 
#' (e.g. \code{'COL05'}) will be assigned to the existing columns. If a column
#' name is not specified, but the object has existing column names, the 
#' inserted column will be assigned a default column name.
#' @param data Optional vector of data to insert in the new column. The length 
#' of this vector should be evenly divisible by the number of rows in the data 
#' frame.
#'     
#' @return Input object with column inserted.
#' 
#' @keywords internal
#' @rdname insertColumn
insertColumn <- function(x, col.index, col.name=NULL, data=NA) {
    
    stopifnot( is.data.frame(x) || is.matrix(x) )
    stopifnot( isSingleWholeNumber(col.index) )
    stopifnot( inRange(col.index, c(1, ncol(x) + 1)) )
    
    stopifnot( isWholeNumber( nrow(x) / length(data) ) )
    
    others <- otherattributes(x)
    
    # Get column names of input object.
    object.colnames <- colnames(x)
    
    # Get number of columns before and after column insertion.
    prev.ncol <- ncol(x)
    post.ncol <- prev.ncol + 1

    # If this is a matrix, insert new column to   
    # the right of the rightmost existing column.
    if ( is.matrix(x) ) {
        x <- cbind( x, rep( NA, nrow(x) ) )
    }
    
    # If new column index is within previously existing columns, 
    # nudge subsequent columns one column to the right.
    if ( col.index <= prev.ncol ) {
        x[, (col.index + 1):post.ncol] <- x[, col.index:prev.ncol]
    }
    
    # Set new column data.
    x[, col.index] <- data

    # Update column names, if appropriate.
    if ( ! is.null(col.name) && ! is.null(object.colnames) ) {
        colnames(x) <- append(object.colnames, col.name, after=(col.index - 1))
    }
    
    otherattributes(x) <- others
    
    return(x)
}

# isBOOL -----------------------------------------------------------------------
#' Test for a single logical value.
#' 
#' @param x Test object.
#'      
#' @return \code{TRUE} if object is a single logical value;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isBOOL
isBOOL <- function(x) {
    return( isTRUE(x) | isFALSE(x) )
}

# isDefaultMarkerID ------------------------------------------------------------
#' Test for default marker IDs.
#'  
#' @param ids Vector of locus IDs.
#'      
#' @return Logical vector indicating which elements are valid
#' IDs that follow the pattern of a default marker ID.
#' 
#' @template section-default-marker-ids
#' 
#' @export
#' @rdname isDefaultMarkerID
isDefaultMarkerID <- function(ids) {
    return( isValidID(ids) & grepl(const$pattern$default.marker.id, ids) )
}

# isDefaultQTLName -------------------------------------------------------------
#' Test for default QTL names.
#' 
#' @param ids Vector of item IDs.
#'      
#' @return Logical vector indicating which elements are valid
#' IDs that follow the pattern of a default QTL name.
#' 
#' @template section-qtl-names
#' 
#' @export
#' @rdname isDefaultQTLName
isDefaultQTLName <- function(ids) {
    return( isValidID(ids) & grepl(const$pattern$default.qtl.name, ids) )
}

# isEnumAllele -----------------------------------------------------------------
#' Test if symbol is a valid enumerated allele.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if symbol is a valid enumerated allele;
#' \code{FALSE} otherwise.
#' 
#' @template section-geno-enum
#' 
#' @keywords internal
#' @rdname isEnumAllele
isEnumAllele <- function(x) {
    return( x %in% const$enum.geno.charset )
}

# isEnumGenotype ---------------------------------------------------------------
#' Test if symbol is a valid enumerated genotype.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if symbol is a valid enumerated
#'  genotype; \code{FALSE} otherwise.
#' 
#' @template section-geno-enum
#' 
#' @keywords internal
#' @rdname isEnumGenotype
isEnumGenotype <- function(x) {
    return( x %in% const$enum.geno.charset )
}

# isFALSE ----------------------------------------------------------------------
#' Test for a single \code{FALSE} value.
#' 
#' @param x Test object.
#'      
#' @return \code{TRUE} if object is a single \code{FALSE} value;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isFALSE
isFALSE <- function(x) {
    return( identical(FALSE, x) )
}

# isFounderAllele --------------------------------------------------------------
#' Test if symbol is a valid founder allele.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if symbol is a valid founder allele;
#' \code{FALSE} otherwise.
#' 
#' @template section-geno-founder
#' 
#' @keywords internal
#' @rdname isFounderAllele
isFounderAllele <- function(x) {
    return( x %in% const$founder.allele.charset )
}

# isFounderGenotype ------------------------------------------------------------
#' Test if symbol is a valid founder genotype.
#' 
#' @param x Test object.
#' @param strict Return \code{TRUE} only for complete genotypes.
#' 
#' @return \code{TRUE} if symbol is a valid founder genotype;
#' \code{FALSE} otherwise.
#' 
#' @template section-geno-founder
#' 
#' @keywords internal
#' @rdname isFounderGenotype
isFounderGenotype <- function(x, strict=FALSE) {
    
    allele.list <- strsplit(x, '')
    
    strict.charset <- const$founder.allele.charset
    
    if (strict) {
        
        result <- unlist( lapply(allele.list, function(alleles)
            all( alleles %in% strict.charset ) ) )
        
    } else {
        
        extended.charset <- c(strict.charset, const$missing.value)
        result <- unlist( lapply(allele.list, function(alleles)
            all( alleles %in% extended.charset ) &&
            any( alleles %in% strict.charset ) ) )
    }
    
    return(result)
}

# isMarkerID -------------------------------------------------------------------
#' Test for marker IDs.
#'  
#' @param loc.ids Vector of locus IDs.
#'      
#' @return Logical vector indicating which elements are valid IDs that do
#' \emph{not} follow the pattern of a pseudomarker ID.
#' 
#' @template section-marker-ids
#' 
#' @export
#' @rdname isMarkerID
isMarkerID <- function(loc.ids) {
    return( isValidID(loc.ids) & ! grepl(const$pattern$pseudomarker.id, loc.ids) )
}

# isNonNegativeNumber ----------------------------------------------------------
#' Test for non-negative numbers.
#' 
#' @param n Test vector.
#'      
#' @return Logical vector indicating which elements are greater than or equal to 
#' zero.
#' 
#' @keywords internal
#' @rdname isNonNegativeNumber
isNonNegativeNumber <- function(n) {
    return( is.numeric(n) & is.finite(n) & n >= 0 )
}

# isPositiveNumber -------------------------------------------------------------
#' Test for positive numbers.
#' 
#' @param n Test vector.
#'      
#' @return Logical vector indicating which elements are positive numbers.
#' 
#' @keywords internal
#' @rdname isPositiveNumber
isPositiveNumber <- function(n) {
    return( is.numeric(n) & is.finite(n) & n > 0 )
}

# isPositiveWholeNumber --------------------------------------------------------
#' Test for positive whole numbers.
#' 
#' @param n Test vector.
#' @param tol Numeric tolerance.
#'      
#' @return Logical vector indicating which elements are positive whole numbers.
#' 
#' @keywords internal
#' @rdname isPositiveWholeNumber
isPositiveWholeNumber <- function(n, tol=.Machine$double.eps^0.5) {
    stopifnot( tol > 0 )
    return( is.numeric(n) & is.finite(n) & n > tol & abs(n - round(n)) < abs(tol) )
}

# isProbability ----------------------------------------------------------------
#' Test for valid probabilities.
#' 
#' @param n Test vector.
#'      
#' @return Logical vector indicating which elements are valid probabilities.
#' 
#' @keywords internal
#' @rdname isProbability
isProbability <- function(n) {
    return( is.numeric(n) & is.finite(n) & n >= 0 & n <= 1 )
}

# isPseudomarkerID -------------------------------------------------------------
#' Test for pseudomarker IDs.
#'  
#' @param loc.ids Vector of locus IDs.
#'      
#' @return Logical vector indicating which elements are valid IDs that follow
#' the pattern of a pseudomarker ID.
#' 
#' @template section-pseudomarker-ids
#' 
#' @export
#' @rdname isPseudomarkerID
isPseudomarkerID <- function(loc.ids) {
    return( isValidID(loc.ids) & grepl(const$pattern$pseudomarker.id, loc.ids) )
}

# isRawAllele ------------------------------------------------------------------
#' Test if symbol is a valid raw SNP allele.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if symbol is a valid raw SNP allele;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isRawAllele
isRawAllele <- function(x) {
    return( x %in% const$raw.allele.charset )
}

# isRawGenotype ----------------------------------------------------------------
#' Test if symbol is a valid raw SNP genotype.
#' 
#' @param x Test object.
#' @param strict Return \code{TRUE} only for complete genotypes.
#' 
#' @return \code{TRUE} if symbol is a valid raw SNP genotype;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isRawGenotype
isRawGenotype <- function(x, strict=FALSE) {
    
    allele.list <- strsplit(x, '')
    
    strict.charset <- const$raw.allele.charset
    
    if (strict) {
    
        result <- unlist( lapply(allele.list, function(alleles)
            all( alleles %in% strict.charset ) ) )
        
    } else {
        
        extended.charset <- c(strict.charset, const$missing.value)
        result <- unlist( lapply(allele.list, function(alleles)
            all( alleles %in% extended.charset ) &&
            any( alleles %in% strict.charset ) ) )
    }
    
    return(result)
}

# isSingleChar -----------------------------------------------------------------
#' Test for a single character.
#' 
#' @param x Test object.
#'      
#' @return \code{TRUE} if the object is a single character;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleChar
isSingleChar <- function(x) {
    return( length(x) == 1 && is.character(x) && ! is.na(x) && nchar(x) == 1 )
}

# isSingleFiniteNumber ---------------------------------------------------------
#' Test for a single finite number.
#' 
#' @param n Test object.
#' 
#' @return \code{TRUE} if the object is a single finite number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleFiniteNumber
isSingleFiniteNumber <- function(n) {
    return( length(n) == 1 && is.numeric(n) && is.finite(n) )
}

# isSingleNonNegativeNumber ----------------------------------------------------
#' Test for a single non-negative number.
#' 
#' @param n Test object.
#'      
#' @return \code{TRUE} if the object is a single non-negative number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleNonNegativeNumber
isSingleNonNegativeNumber <- function(n) {
    return( length(n) == 1 && is.numeric(n) && is.finite(n) && n >= 0 )
}

# isSingleNonNegativeWholeNumber -----------------------------------------------
#' Test for a single non-negative whole number.
#' 
#' @param n Test object.
#' @param tol Numeric tolerance.
#'      
#' @return \code{TRUE} if the object is a single non-negative whole number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleNonNegativeWholeNumber
isSingleNonNegativeWholeNumber <- function(n, tol=.Machine$double.eps^0.5) {
    return( length(n) == 1 && is.numeric(n) && is.finite(n) && 
                n >= 0 && abs(n - round(n)) < abs(tol) )
}

# isSinglePositiveNumber -------------------------------------------------------
#' Test for a single positive number.
#' 
#' @param n Test object.
#'      
#' @return \code{TRUE} if the object is a single positive number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSinglePositiveNumber
isSinglePositiveNumber <- function(n) {
    return( length(n) == 1 && is.numeric(n) && is.finite(n) && n > 0 )
}

# isSinglePositiveWholeNumber --------------------------------------------------
#' Test for a single positive whole number.
#' 
#' @param n Test object.
#' @param tol Numeric tolerance.
#'      
#' @return \code{TRUE} if the object is a single positive whole number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSinglePositiveWholeNumber
isSinglePositiveWholeNumber <- function(n, tol=.Machine$double.eps^0.5) {
    return( length(n) == 1 && is.numeric(n) && is.finite(n) && 
        n > tol && abs(n - round(n)) < abs(tol) )
}

# isSingleProbability ----------------------------------------------------------
#' Test for a single valid probability.
#' 
#' @param n Test object.
#'      
#' @return \code{TRUE} if the object is a single valid probability;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleProbability
isSingleProbability <- function(n) {
    return( length(n) == 1 && is.numeric(n) && 
        is.finite(n) && n >= 0 && n <= 1 )
}

# isSingleString ---------------------------------------------------------------
#' Test for a single character string.
#' 
#' @param x Test object.
#'      
#' @return \code{TRUE} if the object is a single string;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleString
isSingleString <- function(x) {
    return( length(x) == 1 && is.character(x) && ! is.na(x) )
}

# isSingleWholeNumber ----------------------------------------------------------
#' Test for a single whole number.
#' 
#' @param n Test object.
#' @param tol Numeric tolerance.
#'      
#' @return \code{TRUE} if the object is a single whole number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleWholeNumber
isSingleWholeNumber <- function(n, tol=.Machine$double.eps^0.5) {
    return( length(n) == 1 && is.numeric(n) && is.finite(n) && 
        abs(n - round(n)) < abs(tol) )
}

# isValidAllele ----------------------------------------------------------------
#' Test if symbol is a valid allele.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if symbol is a valid enumerated or founder allele;
#' \code{FALSE} otherwise.
#' 
#' @template section-geno-enum
#' @template section-geno-founder
#' 
#' @keywords internal
#' @rdname isValidAllele
isValidAllele <- function(x) {
    return( isEnumAllele(x) | isFounderAllele(x) )
}

# isValidGenotype --------------------------------------------------------------
#' Test if symbol is a valid genotype.
#' 
#' @param x Test object.
#' @param strict Return \code{TRUE} only for complete genotypes.
#' 
#' @return \code{TRUE} if symbol is a valid enumerated or founder genotype;
#' \code{FALSE} otherwise.
#' 
#' @template section-geno-enum
#' @template section-geno-founder
#' 
#' @keywords internal
#' @rdname isValidGenotype
isValidGenotype <- function(x, strict=FALSE) {
    return( isEnumGenotype(x) | isFounderGenotype(x, strict=strict) )
}

# isValidID --------------------------------------------------------------------
#' Test identifier validity.
#'  
#' @param x Test vector.
#'      
#' @return Logical vector indicating which elements are valid identifiers.
#' 
#' @template section-item-ids
#' 
#' @export
#' @rdname isValidID
isValidID <- function(x) {
    return( is.character(x) & nzchar(x) & grepl(const$pattern$item.id, x) )
}

# isValidName ------------------------------------------------------------------
#' Test for syntactically valid names.
#'   
#' @param x Test vector.
#'      
#' @return Logical vector indicating which elements are syntactically valid R
#' names.
#' 
#' @template section-syntactic-names
#' 
#' @export
#' @rdname isValidName
isValidName <- function(x) {
    dotvar.pattern <- '^(?:[.]{3})|(?:[.]{2}[[:digit:]]+)$'
    return( is.character(x) & x == make.names(x) & ! grepl(dotvar.pattern, x) )
}

# isWholeNumber ----------------------------------------------------------------
#' Test for whole numbers.
#' 
#' @param n Test vector.
#' @param tol Numeric tolerance.
#'      
#' @return Logical vector indicating which elements are whole numbers.
#' 
#' @keywords internal
#' @rdname isWholeNumber
isWholeNumber <- function(n, tol=.Machine$double.eps^0.5) { 
    return( is.numeric(n) & is.finite(n) & abs(n - round(n)) < abs(tol) )
}

# loadChrTable -----------------------------------------------------------------
#' Load chromosome info.
#' 
#' This function loads standard chromosome info.
#' 
#' @return A \code{data.frame} object containing chromosome info.
#' 
#' @importFrom utils read.csv
#' @keywords internal
#' @rdname loadChrTable
loadChrTable <- function() {
    
    filepath <- requestPkgDataPath('extdata', 'genomes', 'chrtab.csv')
    
    column.classes <- c(seqids='character', seqnames='character', aliases='character', 
        isCircular='logical', genome='character')
    
    chrtab <- utils::read.csv(filepath, quote='', stringsAsFactors=FALSE,
        strip.white=TRUE, na.strings='', colClasses=column.classes)
    
    return(chrtab)
}

# loadMapping ------------------------------------------------------------------
#' Load a mapping from a line or file.
#' 
#' Given a specified line of text or filename, the input text is loaded
#' as a YAML mapping, and returned as a simple \code{mapping} object.
#' 
#' @param line Character string containing the line to be loaded.
#' @param file Character string containing the name of the file to be loaded.
#' 
#' @return A \code{mapping} of unique string keys to individual values.
#' 
#' @importFrom yaml yaml.load
#' @importFrom yaml yaml.load_file
#' @keywords internal
#' @rdname loadMapping
loadMapping <- function(line=NULL, file=NULL) {
    
    if ( ! xor( is.null(line), is.null(file) ) ) {
        stop("mapping must be loaded from either a line or a file, but not both")
    }
    
    if ( ! is.null(line) ) {
        
        stopifnot( isSingleString(line) )
        
        # Ensure line is enclosed by curly braces.
        first.char <- substr(line, 1, 1)
        last.char <- substr(line, nchar(line), nchar(line))
        if ( first.char != '{' || first.char != '}' ) {
            line <- paste0('{', line, '}', collapse='')
        }
        
        tryCatch({ # Load YAML object from line.
            x <- yaml::yaml.load(line)
        }, error=function(e) {
            stop("failed to load YAML from line - '", toString(line), "'")
        })
        
    } else if ( ! is.null(file) ) {
        
        stopifnot( isSingleString(file) )
        
        tryCatch({ # Load YAML object from file.
            x <- yaml::yaml.load_file(file)
        }, error=function(e) {
            stop("failed to load YAML from file - '", toString(file), "'")
        })
    }
    
    tryCatch({ # Convert YAML object to mapping.
        x <- mapping(x)
    }, error=function(e) {
        stop("failed to load mapping")
    })
    
    return(x)
}

# loadSeqTables ----------------------------------------------------------------
#' Load genome sequence info tables.
#' 
#' This function loads the sequence info tables of genomes for which package
#' data is available.
#' 
#' @return A \code{list} of \code{data.frame} objects, each element named for a
#' given genome and containing sequence info for that genome. Each such
#' \code{data.frame} object can be used to create, but is distinct from,
#' a \pkg{GenomeInfoDb} \code{Seqinfo} object.
#' 
#' @importFrom utils read.csv
#' @keywords internal
#' @rdname loadSeqTables
loadSeqTables <- function() {
    
    column.classes <- c(seqids='character', seqnames='character', 
        seqlengths='integer', maplengths='numeric', isCircular='logical', 
        genome='character')
    
    genome.root <- requestPkgDataPath('extdata', 'genomes')
    
    genomes <- list.dirs(genome.root, full.names=FALSE, recursive=FALSE)
    
    seqtab <- vector('list', length(genomes))
    names(seqtab) <- genomes
    
    for ( genome in genomes ) {
        
        filepath <- file.path(genome.root, genome, 'seqtab.csv')
        
        seqtab[[genome]] <- utils::read.csv(filepath, quote='',
            stringsAsFactors=FALSE, strip.white=TRUE, na.strings='', 
            colClasses=column.classes)
    }
    
    return(seqtab)
}

# loadVector -------------------------------------------------------------------
#' Load a vector from a line or file.
#' 
#' Given a specified line of text or filename, the input text is loaded
#' as a YAML list, and returned as an R \code{vector} object.
#' 
#' @param line Character string containing the line to be loaded.
#' @param file Character string containing the name of the file to be loaded.
#' 
#' @return Vector of loaded data.
#' 
#' @importFrom utils read.table
#' @importFrom yaml yaml.load
#' @importFrom yaml yaml.load_file
#' @keywords internal
#' @rdname loadVector
loadVector <- function(line=NULL, file=NULL, type=NULL) {
    
    if ( ! xor( is.null(line), is.null(file) ) ) {
        stop("vector must be loaded from either a line or a file, but not both")
    }
    
    if ( ! is.null(line) ) {
        
        stopifnot( isSingleString(line) )
        
        # Ensure line is enclosed by square brackets.
        first.char <- substr(line, 1, 1)
        last.char <- substr(line, nchar(line), nchar(line))
        if ( first.char != '[' || first.char != ']' ) {
            line <- paste0('[', line, ']', collapse='')
        }
        
        tryCatch({ # Load YAML list from line.
            x <- yaml::yaml.load(line)
            stopifnot( is.vector(x) && ! hasNames(x) )
        }, error=function(e) {
            stop("failed to load vector from line - '", toString(line), "'")
        })
        
    } else if ( ! is.null(file) ) {
        
        stopifnot( isSingleString(file) )
        
        tryCatch({
        
            # Read list file as a CSV file with one column.
            x <- utils::read.table(file, sep='\n', quote='', as.is=TRUE,
                na.strings=NA, colClasses='character', strip.white=TRUE,
                blank.lines.skip=FALSE, allowEscapes=TRUE,
                stringsAsFactors=FALSE)
            stopifnot( ncol(x) == 1 )
            
            # Trim any blank rows from the bottom.
            x <- bstripBlankRows(x)
            
            # Convert to vector of YAML-loaded values.
            x <- lapply(getRowIndices(x), function(i) yaml::yaml.load(x[i, ]))
            stopifnot( all( lengths(x) == 1 ) )
            
            # Get data types of YAML-loaded elements.
            element.types <- unique( sapply(x, typeof) )
            
            if ( ! all( element.types %in%
                c('character', 'logical', 'integer', 'double', 'numeric') ) ) {
                stop("unknown vector element types - '", toString(element.types), "'")
            }
            
            # If all elements are of the same type,
            # convert list to a vector of that type.
            if ( length(element.types) == 1 ) {
                
                if ( element.types == 'character' ) {
                    x <- as.character(x)
                } else if ( element.types == 'logical' ) {
                    x <- as.logical(x)
                } else if ( element.types == 'integer' ) {
                    x <- as.integer(x)
                } else if ( element.types %in% c('double', 'numeric') ) {
                    x <- as.numeric(x)
                }
            }
            
        }, error=function(e) {
            stop("failed to load vector from file - '", toString(file), "'")
        })
    }
    
    # If relevant, check type of vector against
    # that specified, and convert if appropriate.
    if ( ! is.null(type) ) {
        
        stopifnot( isSingleString(type) )
        
        err.type <- NULL
        if ( type == 'character' ) {
            
            if ( ! is.character(x) ) {
                x <- as.character(x)
            }
            
        } else if ( type == 'logical' ) {
            
            if ( ! is.logical(x) ) {
                err.type <- 'logical'
            }
            
        } else if ( type == 'integer' ) {
            
            if ( ! is.integer(x) ) {
                if ( is.numeric(x) ) {
                    x <- as.integer(x)
                } else {
                    err.type <- 'integer'
                }
            }
            
        } else if ( type %in% c('double', 'numeric') ) {
            
            if ( ! is.numeric(x) ) {
                err.type <- 'numeric'
            }
            
        } else if ( type == 'list' ) {
            
            if ( ! is.list(x) ) {
                x <- as.list(x)
            }
            
        } else {
            
            stop("unknown vector type - '", type, "'")
        }
        
        if ( ! is.null(err.type) ) {
            stop("vector is not of type - '", type, "'")
        }
    }
    
    return(x)
}

# makeAlleleSet ----------------------------------------------------------------
#' Make allele set from the given genotypes.
#' 
#' @param genotypes Sorted character vector of genotype symbols.
#' 
#' @return Character vector of allele symbols generated from the given genotypes.
#' 
#' @keywords internal
#' @rdname makeAlleleSet
makeAlleleSet <- function(genotypes) {
    validateGenotypeSet(genotypes)
    alleles <- sort( unique( unlist( strsplit(genotypes, '') ) ) )
    alleles <- alleles[ alleles != const$missing.value ]
    return(alleles)
}

# makeDefaultMarkerIDs ---------------------------------------------------------
#' Make default marker IDs for loci.
#' 
#' @description This function creates marker IDs from locus \code{data.frame}
#' \code{'loc'}, and validates that each locus has a map position that
#' is within range for the given reference sequence/chromosome. If any
#' locus position is found to be out-of-range, then it will stop with
#' an error message.
#' 
#' @param loc Locus \code{data.frame}, with columns \code{'chr'} and
#' \code{'pos'}, specifying physical map positions.
#' @param sep Separator between the two parts of a default marker ID.
#'      
#' @return Character vector of default marker IDs.
#' 
#' @template author-thomas-walsh
#' @template author-yue-hu
#' @template section-default-marker-ids
#' 
#' @export
#' @rdname makeDefaultMarkerIDs
makeDefaultMarkerIDs <- function(loc, sep=c(':', '-')) {
    
    stopifnot( is.data.frame(loc) )
    stopifnot( nrow(loc) > 0 )
    
    sep <- match.arg(sep)
    
    # Check map unit is physical, if present.
    map.unit <- getMapUnit(loc)
    if ( ! is.na(map.unit) ) {
        
        stopifnot( isPhysicalMapUnit(map.unit) )
        
        # Scale locus positions to base-pair units, if needed.
        basic.unit <- const$basic.map.unit[['pmap']]
        if ( map.unit != basic.unit ) {
            loc <- setMapUnit(loc, basic.unit)
        }
    }
    
    # Get locus sequences and positions.
    loc.seqs <- normSeq( pullLocusSeq(loc) )
    loc.pos <- pullLocusPos(loc)
    
    # Validate positions from genome sequence lengths.
    seqtab <- getSeqTable()
    loc.seq.lengths <- seqtab$seqlengths[ match(loc.seqs,seqtab$seqids) ]
    exrange <- loc.pos[ loc.pos < 1 | loc.pos > loc.seq.lengths ]
    if ( length(exrange) > 0 ) {
        stop("cannot make default marker IDs for positions '", toString(exrange), "'")
    }
    
    # Format locus positions.
    loc.pos <- sprintf('%07d', loc.pos)

    return( paste0("c", loc.seqs, sep, loc.pos) )
}

# makeDefaultQTLNames ----------------------------------------------------------
#' Make default QTL names for the specified loci.
#'  
#' @param loc Locus \code{data.frame}, with columns \code{'chr'} and
#' \code{'pos'}, specifying genetic map positions.
#' @param step Map step size.
#'      
#' @return Character vector of default QTL names.
#' 
#' @template section-qtl-names
#' 
#' @export
#' @rdname makeDefaultQTLNames
makeDefaultQTLNames <- function(loc, step=0) {
    
    stopifnot( is.data.frame(loc) )
    stopifnot( nrow(loc) > 0 )
    stopifnot( isSingleNonNegativeNumber(step) )
    
    # Check map unit is genetic, if present.
    map.unit <- getMapUnit(loc)
    if ( ! is.na(map.unit) ) {
        stopifnot( isGeneticMapUnit(map.unit) )
    }
    
    # Adjust number of digits in map position. 
    # NB: based on code in R/qtl::makeqtl
    digits <- 1
    if ( step > 0 ) {
        digits <- max( digits, -floor( log10(step) ) )
    }
    
    # Get locus sequences and positions.
    loc.seqs <- normSeq( pullLocusSeq(loc) )
    loc.pos <- pullLocusPos(loc)
    
    # Format locus positions.
    loc.pos <- sprintf(paste0('%.', digits, 'f'), loc.pos)
    
    return( paste(loc.seqs, loc.pos, sep='@') )
}

# makeGenotypeSet --------------------------------------------------------------
#' Make genotype set from the given alleles and cross type.
#' 
#' @param alleles Sorted character vector containing a set of allele symbols.
#' @param crosstype Cross type for which genotype set should be made.
#' @param strict Make genotype set containing only complete genotypes.
#' 
#' @return Character vector of genotype symbols generated from the
#' given alleles in a way determined by the specified cross type.
#' 
#' @keywords internal
#' @rdname makeGenotypeSet
makeGenotypeSet <- function(alleles, crosstype, strict=FALSE) {
    
    # TODO: implement for other cross types.
    # TODO: implement strict option. By default (strict=FALSE), partially
    # missing genotype symbols (e.g. 'A-') are included in the returned
    # genotype set. If strict is TRUE, the returned genotype set contains
    # only complete genotype symbols.
    
    validateAlleleSet(alleles)
    
    if ( crosstype == 'haploid' ) {
        
        # NB: strict option is irrelevant for a haploid cross type,
        # because partial genotypes cannot exist by definition.
        genotypes <- alleles
        
    } else if ( crosstype %in% const$supported.crosstypes ) {
        stop("cannot get genotype set for cross type - '", crosstype, "'")
    } else {
        stop("unsupported cross type - '", crosstype, "'")
    }
    
    return(genotypes)
}

# makeNumbers ------------------------------------------------------------------
#' Make numbers from numeric names.
#' 
#' @param x Vector of syntactically valid names that represent numbers.
#'      
#' @return Numeric vector of numbers corresponding to the input names.
#' 
#' @keywords internal
#' @rdname makeNumbers
makeNumbers <- function(x) {
    
    stopifnot( isValidName(x) )
    
    m <- regexec(const$pattern$numeric.name, x)
    matches <- regmatches(x, m)
    
    unconvertible <- x[ lengths(matches) == 0 ]
    if ( length(unconvertible) > 0 ) {
        stop("cannot make numbers from non-numeric names - '", 
            toString(unconvertible), "'")
    }
    
    is.negative <- sapply(matches, getElement, 2) == '.'
    numbers <- sapply(matches, function(regmatch) as.numeric(regmatch[3]))
    
    numbers[is.negative] <- -1 * numbers[is.negative]
    
    return(numbers)
}

# makePseudomarkerIDs ----------------------------------------------------------
#' Make pseudomarker IDs for the specified loci.
#'  
#' @param loc Locus \code{data.frame}, with columns \code{'chr'} and
#' \code{'pos'}, specifying genetic map positions.
#'      
#' @return Character vector of pseudomarker IDs.
#' 
#' @template section-pseudomarker-ids
#' 
#' @export
#' @rdname makePseudomarkerIDs
makePseudomarkerIDs <- function(loc) {
    
    stopifnot( is.data.frame(loc) )
    stopifnot( nrow(loc) > 0 )
    
    # Check map unit is genetic, if present.
    map.unit <- getMapUnit(loc)
    if ( ! is.na(map.unit) ) {
        stopifnot( isGeneticMapUnit(map.unit) )
    }
    
    # Get locus sequences and positions.
    loc.seqs <- normSeq( pullLocusSeq(loc) )
    loc.pos <- pullLocusPos(loc)
    
    return( paste0('c', loc.seqs, '.loc', loc.pos) )
}

# otherattributes --------------------------------------------------------------
#' Get non-reserved object attributes. 
#' 
#' @details As with R base function \code{mostattributes}, this provides access
#' to non-reserved attributes.
#' 
#' @param x R object.
#' 
#' @return List of non-reserved object attributes.
#' 
#' @keywords internal
#' @rdname otherattributes
otherattributes <- function(x) {
    special.attributes <- getSpecialAttributeNames(x)
    return( attributes(x)[ ! names(attributes(x)) %in% special.attributes ] )
}

# `otherattributes<-` ----------------------------------------------------------
#' Set non-reserved object attributes. 
#' 
#' @details As with R base function \code{mostattributes}, this provides access
#' to non-reserved attributes.
#' 
#' @param x R object.
#' @param value List of attributes to set.
#' 
#' @return R object with non-reserved attributes set.
#' 
#' @keywords internal
#' @rdname otherattributes
`otherattributes<-` <- function(x, value) {
    
    stopifnot( hasNames(value) )
    
    special.attributes <- getSpecialAttributeNames(x)
    
    other.names <- names(value)[ ! names(value) %in% special.attributes ]
    for ( other.name in other.names ) {
        attr(x, other.name) <- value[[other.name]]
    }
    
    return(x)
}

# parseDefaultMarkerIDs --------------------------------------------------------
#' Parse default marker IDs.
#' 
#' This function parses an input vector of default marker IDs, and returns a
#' physical \code{mapframe} with the locus in each row derived from the
#' corresponding marker ID in the input vector. An error is raised if any
#' of the input values cannot be parsed as a default marker ID.
#'  
#' @param marker.ids Vector of default marker IDs.
#' 
#' @return Physical \code{mapframe} with loci
#' corresponding to the specified marker IDs.
#' 
#' @template section-default-marker-ids
#' 
#' @export
#' @rdname parseDefaultMarkerIDs
parseDefaultMarkerIDs <- function(marker.ids) {
    
    stopifnot( all( isDefaultMarkerID(marker.ids) ) )
    
    m <- regexec(const$pattern$default.marker.id, marker.ids)
    regmatch.list <- regmatches(marker.ids, m)
    
    marker.seqs <- sapply(regmatch.list, getElement, 2)
    marker.pos <- sapply(regmatch.list, function(x) as.numeric(x[3]))
    
    return( mapframe(chr=marker.seqs, pos=marker.pos,
        row.names=marker.ids, map.unit='bp') )
}

# parseDefaultQTLNames ---------------------------------------------------------
#' Parse default QTL names.
#' 
#' This function parses an input vector of default QTL names, and returns a
#' genetic \code{mapframe} with the locus in each row derived from the
#' corresponding QTL name in the input vector. An error is raised if any
#' of the input values cannot be parsed as a default QTL name.
#' 
#' @param qtl.names Vector of default QTL names.
#' 
#' @return Genetic \code{mapframe} with loci
#' corresponding to the specified QTL names.
#' 
#' @template section-qtl-names
#' 
#' @export
#' @rdname parseDefaultQTLNames
parseDefaultQTLNames <- function(qtl.names) {
    
    stopifnot( all( isDefaultQTLName(qtl.names) ) )
    
    m <- regexec(const$pattern$default.qtl.name, qtl.names)
    regmatch.list <- regmatches(qtl.names, m)
    
    qtl.seqs <- sapply(regmatch.list, getElement, 2)
    qtl.pos <- sapply(regmatch.list, function(x) as.numeric(x[3]))
    
    return( gmapframe(chr=qtl.seqs, pos=qtl.pos, row.names=qtl.names) )
}

# parseFilenames ---------------------------------------------------------------
#' Parse filenames by pattern matching.
#' 
#' @param filenames Character vector of filenames to parse.
#' @param filename.pattern Filename pattern, which must be a valid Perl regex
#' with named capture groups. Neither the capture groups nor the regex itself
#' are required to match any given filename, but all capture groups must have
#' a unique name.
#' 
#' @return Character matrix in which each row contains the values of capture
#' groups for a given filename, and each column contains values of a given
#' capture group across files. Unmatched capture groups are represented by
#' \code{NA} values.
#' 
#' @keywords internal
#' @rdname parseFilenames
parseFilenames <- function(filenames, pattern) {
    
    stopifnot( is.character(filenames) )
    stopifnot( isSingleString(pattern) )
    
    parsed <- as.list( structure(rep(NA_character_,
        length(filenames)), names=filenames) )
    
    capture.names <- NULL
    
    # Parse each filename.
    for ( filename in filenames ) {
        
        # Apply regular expression to filename.
        pattern.match <- regexpr(pattern, filename, perl=TRUE)
        
        # Get capture-group names of pattern match.
        capture.names <- attr(pattern.match, 'capture.names')
        stopifnot( length(capture.names) > 0 )
        stopif( any( capture.names == '' ) )
        
        # Get first indices of capture groups.
        capture.first <- attr(pattern.match, 'capture.start')[1, ]
        capture.first[ capture.first %in% c(-1, 0) ] <- NA_integer_
        
        # Get lengths of capture groups.
        capture.length <- attr(pattern.match, 'capture.length')[1, ]
        capture.length[ capture.length %in% c(-1, 0) ] <- NA_integer_
        
        # Get last indices of capture groups.
        capture.last <- capture.first + capture.length - 1
        
        # Get capture-group strings.
        parsed[[filename]] <- substring(filename, capture.first, capture.last)
        names(parsed[[filename]]) <- attr(pattern.match, 'capture.names')
    }
    
    parsed <- matrix( unlist(parsed), nrow=length(parsed),
        byrow=TRUE, dimnames=list(filenames, capture.names) )
    
    return(parsed)
}

# parsePseudomarkerIDs ---------------------------------------------------------
#' Parse pseudomarker IDs.
#' 
#' This function parses an input vector of pseudomarker IDs, and returns a
#' genetic \code{mapframe} with the locus in each row derived from the
#' corresponding pseudomarker ID in the input vector. An error is raised
#' if any of the input values cannot be parsed as a pseudomarker ID.
#' 
#' @param loc.ids Vector of pseudomarker IDs.
#' 
#' @return Genetic \code{mapframe} with loci
#' corresponding to the input pseudomarker IDs.
#' 
#' @template section-pseudomarker-ids
#' 
#' @export
#' @rdname parsePseudomarkerIDs
parsePseudomarkerIDs <- function(loc.ids) {
    
    stopifnot( all( isPseudomarkerID(loc.ids) ) )
    
    m <- regexec(const$pattern$pseudomarker.id, loc.ids)
    regmatch.list <- regmatches(loc.ids, m)
    
    loc.seqs <- sapply(regmatch.list, getElement, 2)
    loc.pos <- sapply(regmatch.list, function(x) as.numeric(x[3]))
    
    return( gmapframe(chr=loc.seqs, pos=loc.pos, row.names=loc.ids) )
}

# removeColsNA -----------------------------------------------------------------
#' Remove columns containing only \code{NA} values.
#'
#' @param x A \code{data.frame} or \code{matrix}.
#'
#' @return Input object without those columns that contain only \code{NA} values.
#'
#' @keywords internal
#' @rdname removeColsNA
removeColsNA <- function(x) {
    stopifnot( is.data.frame(x) || is.matrix(x) )
    stopifnot( nrow(x) > 0 )
    mask <- ! apply( x, 2, function(column) allNA( as.vector(column) ) )
    x <- x[, mask, drop=FALSE]
    return(x)
}

# removeRowsNA -----------------------------------------------------------------
#' Remove rows containing only \code{NA} values.
#'
#' @param x A \code{data.frame} or \code{matrix}.
#'
#' @return Input object without those rows that contain only \code{NA} values.
#'
#' @keywords internal
#' @rdname removeRowsNA
removeRowsNA <- function(x) {
    stopifnot( is.data.frame(x) || is.matrix(x) )
    stopifnot( nrow(x) > 0 )
    mask <- ! apply( x, 1, function(row) allNA( as.vector(row) ) )
    x <- x[mask, , drop=FALSE]
    return(x)
}

# requestNodes -----------------------------------------------------------------
#' Request nodes for parallel execution.
#' 
#' @template param-n.cluster
#'     
#' @return Vector of node names.
#' 
#' @export
#' @rdname requestNodes
requestNodes <- function(n.cluster=1) {

    stopifnot( isSinglePositiveWholeNumber(n.cluster) )
    
    # Get the PBS node file, or NA value if not found.
    node.file <- Sys.getenv("PBS_NODEFILE", NA)

    # If PBS node file specified, create node list from node file..
    if ( ! is.na(node.file) && file.info(node.file)$size > 0 ) {
      
        # Read lines from node file.
        f <- file(node.file, "r")
        node.lines <- readLines(f)
        close(f)
      
        # Create node list with N nodes, where N is the number of processes  
        # requested or the number of nodes in the file, whichever is smaller.
        node.list <- node.lines[1:n.cluster]
        node.list <- node.list[ ! is.na(node.list) ]

    } else { # ..otherwise create node list to run on local host. 
        
        # Get number of cores available on the local host.
        local.cores <- parallel::detectCores(logical=FALSE)
        
        # Cap the number of processes at the number of cores.
        n.cluster <- min(n.cluster, local.cores)
      
        # Set node list from number of processes.
        node.list <- rep('localhost', n.cluster)
    }
    
    return(node.list)
}

# requestPkgDataPath -----------------------------------------------------------
#' Request path of \pkg{shmootl} package data.
#' 
#' @param ... Path elements. These should be specified relative to the package 
#' data root. For the source package, this is the \pkg{shmootl} \code{'inst'} 
#' directory, while for the binary package, this is the same as the package 
#' root.
#' 
#' @return An absolute path to the requested \pkg{shmootl} data file/directory.
#' 
#' @keywords internal
#' @rdname requestPkgDataPath
requestPkgDataPath <- function(...) {
    
    datapath <- file.path(...)
    
    result <- system.file(datapath, package='shmootl')
    
    if ( result == '' ) {
        
        # NB: this assumes the current working directory is 'shmootl/R'.
        result <- normalizePath( file.path('..', 'inst', datapath) )
        
        if ( ! file.exists(result) ) {
            stop("shmootl package data file/directory not found - '", 
                 datapath, "'")
        }
    }
    
    return(result)
}

# resolveQtlIndices ------------------------------------------------------------
#' Resolve indices of QTL object.
#'  
#' @param x A \code{qtl} object.
#' @param qtl.indices QTL indices referring to the given QTL object.
#' 
#' @return QTL indices resolved with respect to the given QTL object.
#' 
#' @keywords internal
#' @rdname resolveQtlIndices
resolveQtlIndices <- function(x, qtl.indices=NULL) {
    
    stopifnot( 'qtl' %in% class(x) )
    
    if ( ! is.null(qtl.indices) ) {
        
        if ( length(qtl.indices) == 0 ) {
            stop("no QTL indices specified")
        }
        
        invalid <- qtl.indices[ ! isPositiveWholeNumber(qtl.indices) ]
        if ( length(invalid) > 0 ) {
            stop("invalid QTL indices - '", toString(invalid), "'")
        }
        
        exrange <- qtl.indices[ qtl.indices < 1 | qtl.indices > x$n.qtl ]
        if ( length(exrange) > 0 ) {
            stop("QTL indices out of range - '", toString(exrange), "'")
        }
        
    } else {
        
        qtl.indices <- 1:x$n.qtl
    }
    
    return(qtl.indices)
}

# rstripBlankCols --------------------------------------------------------------
#' Strip blank columns from right of \code{data.frame}.
#'
#' @param x A \code{data.frame} with columns of type \code{character}.
#'
#' @return Input object in which rightmost blank columns (i.e. columns
#' containing no non-whitespace characters) have been stripped.
#' 
#' @keywords internal
#' @rdname rstripBlankCols
rstripBlankCols <- function(x) {
    stopifnot( is.data.frame(x) )
    stopifnot( all( sapply(x, class) == 'character' ) )
    while( allWhite( as.character(x[, ncol(x)]) ) ) {
        x <- x[, -ncol(x), drop=FALSE]
    }
    return(x)
}

# setColumnFromRownames --------------------------------------------------------
#' Move object row names to the specified column.
#' 
#' @param x A \code{data.frame} or \code{matrix}.
#' @param col.name Name of the column to which row names are to be set.
#'     
#' @return Input object with row names inserted into the leftmost column, and 
#' original row names removed.
#' 
#' @keywords internal
#' @rdname setColumnFromRownames
setColumnFromRownames <- function(x, col.name='rownames') {
    
    stopifnot( is.data.frame(x) || is.matrix(x) )
    stopifnot( hasRownames(x) )
    stopif( col.name %in% colnames(x) )
    
    x <- insertColumn(x, 1, col.name=col.name, data=rownames(x))

    rownames(x) <- NULL

    return(x)
}

# setRownamesFromColumn --------------------------------------------------------
#' Move specified column to object row names.
#' 
#' @param x A \code{data.frame} or \code{matrix}.
#' @param col.name Name of the column from which row names are to be set.
#'     
#' @return Input object with the leftmost column removed, and with row names set
#' from its contents.
#' 
#' @keywords internal
#' @rdname setRownamesFromColumn
setRownamesFromColumn <- function(x, col.name='rownames') {
    
    stopifnot( is.data.frame(x) || is.matrix(x) )
    stopif( anyDuplicated( colnames(x) ) )
    stopifnot( col.name %in% colnames(x) )
    
    rownames(x) <- as.character( x[, col.name] )
    
    x <- deleteColumn(x, col.name=col.name)
    
    return(x)
}

# stopif -----------------------------------------------------------------------
#' Stop if the expression is \code{TRUE}.
#' 
#' @details This function is based on the R base function \code{stopifnot}, but  
#' is defined only for the simple case of a single expression. It should be used
#' to avoid confusing double negatives.
#'  
#' @param expression R expression to be tested.
#' 
#' @keywords internal
#' @rdname stopif
stopif <- function(expression) {
    
    if ( length(expression) > 0 ) {
        
        mc <- match.call()
        
        ch <- deparse(mc[[2]])
        
        if ( ! anyNA(expression) && all(expression) ) {
            
            if ( length(expression) > 1 ) {
                msg <- sprintf("%s are all TRUE", ch)
            } else {
                msg <- sprintf("%s is TRUE", ch)
            }
            
            stop(msg, call. = FALSE)
        }    
    }
    
    return( invisible() )
}

# stripWhite -------------------------------------------------------------------
#' Strip leading/trailing whitespace from elements of a character vector.
#' 
#' @param x Character vector whose elements are to be stripped.
#'     
#' @return Whitespace-stripped copy of the input character vector.
#' 
#' @keywords internal
#' @rdname stripWhite
stripWhite <- function(x) {
    stopifnot( is.character(x) )
    return( gsub( "^[[:space:]]+|[[:space:]]+$", "", x) )
}

# validateAlleleSet ------------------------------------------------------------
#' Validate a set of alleles.
#' 
#' @param x Vector of alleles.
#' 
#' @return \code{TRUE} if allele set is valid; otherwise raises first error.
#' 
#' @keywords internal
#' @rdname validateAlleleSet
validateAlleleSet <- function(x) {
    
    stopifnot( is.character(x) )
    
    # Decompose allele symbols into different types.
    founder.symbols <- x[ isFounderAllele(x) ]
    enum.symbols <- x[ isEnumAllele(x) ]
    invalid.values <- x[ is.na(x) | ! isValidAllele(x) ]
    
    # Check for invalid values.
    if ( length(invalid.values) > 0 ) {
        if ( '' %in% invalid.values ) {
            stop("blank allele values")
        } else {
            stop("invalid allele symbols - '", toString(invalid.values), "'")
        }
    }
    
    # Check for duplicate allele symbols.
    dup.symbols <- unique( x[ duplicated(x) ] )
    if ( length(dup.symbols) > 0 ) {
        stop("duplicate allele symbols - '", toString(dup.symbols), "'")
    }
    
    # Check allele set order.
    if ( is.unsorted(x) ) {
        stop("allele set is not sorted - '", toString(x), "'")
    }
    
    # Check that symbols are either founder or enumerated alleles.
    if ( ! xor(length(founder.symbols) > 0, length(enum.symbols) > 0) ) {
        stop("alleles must be of enumerated or founder type, but not both")
    }
    
    # Verify that there are exactly two alleles.
    # TODO: handle more than two alleles.
    if ( length(x) != 2 ) {
        stop("unsupported number of alleles - '", length(x), "'")
    }
    
    return(TRUE)
}

# validateGenotypeSet ----------------------------------------------------------
#' Validate a set of genotypes.
#' 
#' @param x Vector of genotypes.
#' @param strict Restrict validity to complete genotypes only.
#' 
#' @return \code{TRUE} if genotype set is valid; otherwise raises first error.
#' 
#' @keywords internal
#' @rdname validateGenotypeSet
validateGenotypeSet <- function(x, strict=FALSE) {
    
    stopifnot( is.character(x) )
    
    # Decompose genotype symbols into different types.
    founder.symbols <- x[ isFounderGenotype(x, strict=strict) ]
    enum.symbols <- x[ isEnumGenotype(x) ]
    invalid.values <- x[ is.na(x) | ! isValidGenotype(x, strict=strict) ]
    
    # Check for invalid values.
    if ( length(invalid.values) > 0 ) {
        if ( '' %in% invalid.values ) {
            stop("blank genotype values")
        } else {
            stop("invalid genotype symbols - '", toString(invalid.values), "'")
        }
    }
    
    # Check that symbols are either founder or enumerated genotypes.
    if ( ! xor(length(founder.symbols) > 0, length(enum.symbols) > 0) ) {
        stop("genotypes must be of enumerated or founder type, but not both")
    }
    
    # Check for duplicate genotype symbols.
    dup.symbols <- unique( x[ duplicated(x) ] )
    if ( length(dup.symbols) > 0 ) {
        stop("duplicate genotype symbols - '", toString(dup.symbols), "'")
    }
    
    # Check genotype set order.
    if ( is.unsorted(x) ) {
        stop("genotype set is not sorted - '", toString(x), "'")
    }
    
    # Verify that genotypes are haploid.
    # TODO: handle other ploidies.
    if ( any( nchar(x) > 1 ) ) {
        stop("unsupported genotype ploidy")
    }
    
    # Verify that genotyes have consistent ploidy.
    if ( length( unique( nchar(x) ) ) != 1 ) {
        stop("inconsistent genotype ploidy")
    }
    
    # Verify that there are exactly two genotypes.
    # TODO: handle more than two genotypes.
    if ( length(x) != 2 ) {
        stop("unsupported number of genotypes - '", length(x), "'")
    }
    
    return(TRUE)
}

# End of util.R ################################################################
