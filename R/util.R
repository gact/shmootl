# Start of util.R ##############################################################

# allKwargs --------------------------------------------------------------------
#' Test if ellipsis arguments are all keyword arguments.
#' 
#' @param ... Ellipsis arguments.
#'     
#' @return TRUE if all ellipsis arguments are keyword arguments; FALSE otherwise.
#' 
#' @keywords internal
#' @rdname allKwargs
allKwargs <- function(...) {
    args <- list(...)
    return( length(args) == 0 || hasNames(args) )
}

# allNA ------------------------------------------------------------------------
#' Test if all elements are NA values.
#' 
#' @param x Test object.
#'     
#' @return TRUE if object has nonzero length and all elements are NA values; 
#' FALSE otherwise.
#' 
#' @keywords internal
#' @rdname allNA
allNA <- function(x) {
    na.vals <- is.na(x)
    return( length(na.vals) > 0 && all(na.vals) )
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
    
    # Delete column.
    x <- x[, -col.index, drop=FALSE]
    
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
#' @keywords internal
#' @rdname dispatchFromClassS3
dispatchFromClassS3 <- function(generic, class.vector, package) {
    
    stopifnot( isSingleString(generic) )
    stopifnot( is.character(class.vector) )
    stopifnot( length(class.vector) > 0 )
    stopifnot( isSingleString(package) )
    
    pkg <- paste0('package:', package)
    pattern <- paste0('^', generic, '[.]([^[:space:]]+)')
    x <- lsf.str(pkg, pattern=pattern)
    
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
        
        if ( ! allNA(class.match) ) {
            
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

# getCoercionFromClassS3 -------------------------------------------------------
#' Get coercion function for the given class.
#' 
#' @param class.vector Class vector of an R object.
#'     
#' @return Coercion function.
#' 
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
        roman=as.roman
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

# getFilePrefix ----------------------------------------------------------------
#' Get file prefixes.
#' 
#' @param x Character vector of file paths.
#'     
#' @return Vector in which each element contains the file prefix of the 
#' corresponding element of the input vector, or the input file path if 
#' no file extension was found.
#' 
#' @keywords internal
#' @rdname getFilePrefix
getFilePrefix <- function(x) {
    
    stopifnot( is.character(x) )
    stopifnot( length(x) > 0 )
    
    m <- regexec(const$pattern$file.with.ext, x)
    matches <- regmatches(x, m)
    
    prefixes <- sapply(getIndices(x), function(i) if ( length(matches[[i]]) > 0 ) 
    { matches[[i]][2] } else { x[i] } )
    
    return(prefixes)
}

# getFileExtension -------------------------------------------------------------
#' Get file extensions.
#' 
#' @param x Character vector of file paths.
#'     
#' @return Vector in which each element contains the file extension of the 
#' corresponding element of the input vector, or an empty string if no file
#' extension was found.
#' 
#' @keywords internal
#' @rdname getFileExtension
getFileExtension <- function(x) {
    
    stopifnot( is.character(x) )
    stopifnot( length(x) > 0 )
    
    m <- regexec(const$pattern$file.with.ext, x)
    matches <- regmatches(x, m)
    
    extensions <- sapply(getIndices(x), function(i) if ( length(matches[[i]]) > 0 ) 
    { matches[[i]][3] } else { '' } )
    
    return(extensions)
}

# getIndices -------------------------------------------------------------------
#' Get indices of object.
#' 
#' @param x An object with length.
#'     
#' @return Integer vector of all object indices. Returns an empty integer vector
#' if the object does not have length.
#'    
#' @keywords internal
#' @rdname getIndices
getIndices <- function(x) {
    object.length <- length(x)
    stopifnot( isSingleNonNegativeWholeNumber(object.length) )
    return( if ( object.length > 0 ) { 1:object.length } else { integer() } )
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
#' @return TRUE if object has nonempty names with no NA values; FALSE otherwise.
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
#' @return TRUE if input object has non-default row names; FALSE otherwise.
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
    for ( i in getIndices(size.groups) ) {
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
#' (e.g. 'COL05') will be assigned to the existing columns. If a column 
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
    
    return(x)
}

# isBOOL -----------------------------------------------------------------------
#' Test for a single logical value.
#' 
#' @param x Test object.
#'      
#' @return TRUE if object is a single logical value; FALSE otherwise.
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
#' @return Logical vector indicating which elements are valid IDs that follow 
#' the pattern of a default marker ID.
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
#' @return Logical vector indicating which elements are default QTL names. 
#'  
#' @export
#' @rdname isDefaultQTLName
isDefaultQTLName <- function(ids) {
    return( isValidID(ids) & grepl(const$pattern$default.qtl.name, ids) )
}

# isFALSE ----------------------------------------------------------------------
#' Test for a single FALSE value.
#' 
#' @param x Test object.
#'      
#' @return TRUE if object is a single FALSE value; FALSE otherwise.
#' 
#' @keywords internal
#' @rdname isFALSE
isFALSE <- function(x) {
    return( identical(FALSE, x) )
}

# isMarkerID -------------------------------------------------------------------
#' Test for marker IDs.
#'  
#' @param loc.ids Vector of locus IDs.
#'      
#' @return Logical vector indicating which elements are marker IDs.
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
#' @return Logical vector indicating which elements are pseudomarker IDs.
#'  
#' @export
#' @rdname isPseudomarkerID
isPseudomarkerID <- function(loc.ids) {
    return( isValidID(loc.ids) & grepl(const$pattern$pseudomarker.id, loc.ids) )
}

# isSingleChar -----------------------------------------------------------------
#' Test for a single character.
#' 
#' @param x Test object.
#'      
#' @return TRUE if the object is a single character; FALSE otherwise.
#' 
#' @keywords internal
#' @rdname isSingleChar
isSingleChar <- function(x) {
    return( length(x) == 1 && is.character(x) && ! is.na(x) && nchar(x) == 1 )
}

# isSingleNonNegativeNumber ----------------------------------------------------
#' Test for a single non-negative number.
#' 
#' @param n Test object.
#'      
#' @return TRUE if the object is a single non-negative number; 
#' FALSE otherwise.
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
#' @return TRUE if the object is a single non-negative whole number; 
#' FALSE otherwise.
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
#' @return TRUE if the object is a single positive number; 
#' FALSE otherwise.
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
#' @return TRUE if the object is a single positive whole number; 
#' FALSE otherwise.
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
#' @return TRUE if the object is a single valid probability; 
#' FALSE otherwise.
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
#' @return TRUE if the object is a single string; FALSE otherwise.
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
#' @return TRUE if the object is a single whole number; FALSE otherwise.
#' 
#' @keywords internal
#' @rdname isSingleWholeNumber
isSingleWholeNumber <- function(n, tol=.Machine$double.eps^0.5) {
    return( length(n) == 1 && is.numeric(n) && is.finite(n) && 
        abs(n - round(n)) < abs(tol) )
}

# isValidID --------------------------------------------------------------------
#' Test identifier validity.
#'  
#' @param x Test vector.
#'      
#' @return Logical vector indicating which elements are valid identifiers.
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
#' @return Logical vector indicating which elements are syntactically valid 
#' names.
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

# loadChrInfo ------------------------------------------------------------------
#' Load chromosome info.
#' 
#' This function loads standard chromosome info.
#' 
#' @return A \code{data.frame} object containing chromosome info.
#' 
#' @keywords internal
#' @rdname loadChrInfo
loadChrInfo <- function() {
    
    filepath <- requestPkgDataPath('extdata', 'genomes', 'chrinfo.csv')
    
    column.classes <- c(seqids='character', seqnames='character', aliases='character', 
        isCircular='logical', genome='character')
    
    chrinfo <- read.csv(filepath, quote='', stringsAsFactors=FALSE, 
        strip.white=TRUE, na.strings='', colClasses=column.classes)
    
    return(chrinfo)
}

# loadSeqInfo ------------------------------------------------------------------
#' Load genome sequence info.
#' 
#' This function loads the sequence info of genomes for which package data
#' is available.
#' 
#' @return A \code{list} of \code{data.frame} objects, each element named for a
#' given genome and containing sequence info for that genome.
#' 
#' @keywords internal
#' @rdname loadSeqInfo
loadSeqInfo <- function() {
    
    column.classes <- c(seqids='character', seqnames='character', 
        seqlengths='integer', maplengths='numeric', isCircular='logical', 
        genome='character')
    
    genome.root <- requestPkgDataPath('extdata', 'genomes')
    
    genomes <- list.dirs(genome.root, full.names=FALSE, recursive=FALSE)
    
    seqinfo <- vector('list', length(genomes))
    names(seqinfo) <- genomes
    
    for ( genome in genomes ) {
        
        filepath <- file.path(genome.root, genome, 'seqinfo.csv')
        
        seqinfo[[genome]] <- read.csv(filepath, quote='', 
            stringsAsFactors=FALSE, strip.white=TRUE, na.strings='', 
            colClasses=column.classes)
    }
    
    return(seqinfo)
}

# makeDefaultMarkerIDs ---------------------------------------------------------
#' Make default marker IDs for loci.
#'   
#' @param loc Locus \code{mapframe} specifying physical map positions.
#' @param sep Separator between the two parts of a default marker ID.
#'      
#' @return Character vector of default marker IDs.
#'  
#' @export
#' @rdname makeDefaultMarkerIDs
makeDefaultMarkerIDs <- function(loc, sep=c(':', '-')) {
    
    stopifnot( isPhysicalMapframe(loc) )
    stopifnot( nrow(loc) > 0 )
    
    sep <- match.arg(sep)
    
    # Scale locus positions to base-pair units. 
    if ( getMapUnit(loc) != 'bp' ) {
        loc <- setMapUnit(loc, 'bp')
    }
    
    # TODO: validate positions from genome sequence lengths.
    exrange <- loc$pos[ loc$pos < 1 | loc$pos > 9999999 ]
    if ( length(exrange) > 0 ) {
        stop("cannot make default marker IDs for positions '", toString(exrange), "'")
    }
    
    # Get locus sequences as character vector.
    loc.seqs <- as.character(loc$chr)
    
    # Format locus positions.
    loc.pos <- sprintf('%07d', loc$pos)

    return( paste0("c", loc.seqs, sep, loc.pos) )
}

# makeDefaultQTLNames ----------------------------------------------------------
#' Make default QTL names for the specified loci.
#'  
#' @param loc Locus \code{mapframe} specifying genetic map positions.
#' @param step Map step size.
#'      
#' @return Character vector of default QTL names.
#'  
#' @export
#' @rdname makeDefaultQTLNames
makeDefaultQTLNames <- function(loc, step=0) {
    
    stopifnot( isGeneticMapframe(loc) )
    stopifnot( nrow(loc) > 0 )
    stopifnot( isSingleNonNegativeNumber(step) )
    
    # Adjust number of digits in map position. 
    # NB: based on code in R/qtl::makeqtl
    digits <- 1
    if ( step > 0 ) {
        digits <- max( digits, -floor( log10(step) ) )
    }
    
    # Get locus sequences as character vector.
    loc.seqs <- as.character(loc$chr)
    
    # Format locus positions.
    loc.pos <- sprintf(paste0('%.', digits, 'f'), loc$pos)
    
    return( paste(loc.seqs, loc.pos, sep='@') )
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
#' @param loc Locus \code{mapframe} specifying map positions.
#'      
#' @return Character vector of pseudomarker IDs.
#'  
#' @export
#' @rdname makePseudomarkerIDs
makePseudomarkerIDs <- function(loc) {
    stopifnot( isGeneticMapframe(loc) )
    stopifnot( nrow(loc) > 0 )
    return( paste0('c', as.character(loc$chr), '.loc', loc$pos) )
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

# parsePseudomarkerID ----------------------------------------------------------
#' Parse a pseudomarker ID.
#'    
#' @param loc.id Locus ID to parse.
#' 
#' @return Genetic \code{mapframe} with locus corresponding to the pseudomarker 
#' ID. Returns \code{NULL} if locus ID cannot be parsed as a pseudomarker ID.
#'  
#' @keywords internal
#' @rdname parsePseudomarkerID
parsePseudomarkerID <- function(loc.id) {
    
    stopifnot( isSingleString(loc.id) )
    stopifnot( isValidID(loc.id) )
    
    m <- regexec(const$pattern$pseudomarker.id, loc.id)
    regmatch <- (regmatches(loc.id, m))[[1]]
    
    result <- NULL
    
    if ( length(regmatch) > 0 ) {
        loc.seq <- regmatch[2]
        loc.pos <- as.numeric(regmatch[3])
        result <- gmapframe(chr=loc.seq, pos=loc.pos, row.names=loc.id)
    }
    
    return(result)
}

# parseDefaultMarkerID ---------------------------------------------------------
#' Parse a default marker ID.
#'  
#' @param marker.id Marker ID to parse.
#' 
#' @return Physical \code{mapframe} with locus corresponding to the specified 
#' marker ID. Returns \code{NULL} if marker ID cannot be parsed as a default 
#' marker ID.
#' 
#' @keywords internal
#' @rdname parseDefaultMarkerID
parseDefaultMarkerID <- function(marker.id) {
    
    stopifnot( isSingleString(marker.id) )
    stopifnot( isValidID(marker.id) )
    
    m <- regexec(const$pattern$default.marker.id, marker.id)
    regmatch <- (regmatches(marker.id, m))[[1]]
    
    result <- NULL
    
    if ( length(regmatch) > 0 ) {
        marker.seq <- as.numeric(regmatch[2])
        marker.pos <- as.numeric(regmatch[3])
        result <- mapframe(chr=marker.seq, pos=marker.pos, 
            row.names=marker.id, map.unit='bp')
    }
    
    return(result)    
}

# parseDefaultQTLName ----------------------------------------------------------
#' Parse a default QTL name.
#'    
#' @param qtl.name QTL name to parse.
#' 
#' @return Genetic \code{mapframe} with locus corresponding to the specified QTL
#' name. Returns \code{NULL} if QTL name cannot be parsed as a default QTL name.
#'  
#' @keywords internal
#' @rdname parseDefaultQTLName
parseDefaultQTLName <- function(qtl.name) {

    stopifnot( isSingleString(qtl.name) ) 
    stopifnot( isValidID(qtl.name) )
    
    m <- regexec(const$pattern$default.qtl.name, qtl.name)
    regmatch <- (regmatches(qtl.name, m))[[1]]

    result <- NULL
    
    if ( length(regmatch) > 0 ) {
        qtl.seq <- regmatch[2]
        qtl.pos <- as.numeric(regmatch[3])
        result <- gmapframe(chr=qtl.seq, pos=qtl.pos, row.names=qtl.name)
    }
    
    return(result)    
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
#' Stop if the expression is TRUE.
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

# splitCSL ---------------------------------------------------------------------
#' Split a comma-separated list.
#' 
#' @param x Character vector of length one, containing items as a 
#' comma-separated list (CSL).
#'     
#' @return Character vector with each element containing one item of the input 
#' CSL.
#' 
#' @keywords internal
#' @rdname splitCSL
splitCSL <- function(x) {
    
    stopifnot( isSingleString(x) )
    
    if ( x != '' ) {
        x <- strsplit(x, ',')[[1]]
    }
    
    return(x)
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

# End of util.R ################################################################