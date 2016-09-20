# Start of mapping.R ###########################################################

# mapping ----------------------------------------------------------------------
#' Create a \code{mapping} object.
#' 
#' @param values Vector of mapping values. If the \code{keys} parameter is
#' not set, the names of this vector are used as the keys of the mapping.
#' @param keys Character vector of mapping keys.
#' 
#' @return A \code{mapping} of unique string keys to
#' individual values; equivalent to a named list in
#' which each element is a vector with one element.
#' 
#' @export
#' @family mapping object functions
#' @rdname mapping
mapping <- function(values=NULL, keys=NULL) {
    
    if ( ! is.null(keys) ) {
        
        if ( ! is.null(values) ) {
            
            if ( ! is.vector(values) ) {
                stop("cannot create mapping from values of class '",
                    toString(class(values)), "'")
            }
            
        } else {
            
            values <- as.list( rep_len(NA, length(keys)) )
        }
        
    } else {
     
        if ( ! is.null(values) ) {
            
            if ( ! is.null( names(values) ) && all(
                ! is.na(names(values)) & names(values) != '' ) ) {
                keys <- names(values)
            }
            
        } else {
            
            keys <- character()
            values <- list()
        }
    }
    
    if ( ! is.list(values) ) {
        values <- as.list(values)
    }
    
    object <- values
    class(object) <- c('mapping', 'list')
    names(object) <- keys # validates mapping object
    
    return(object)
}

# as.mapping -------------------------------------------------------------------
#' Convert to a \code{mapping} object.
#' 
#' @param x Named vector.
#' 
#' @return A \code{mapping} corresponding to the input object.
#' 
#' @export
#' @family mapping object functions
#' @rdname as.mapping
as.mapping <- function(x) {
    return( mapping(values=x, keys=names(x)) )
}

# is.mapping -------------------------------------------------------------------
#' Test if object is a \code{mapping}.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if object is a \code{mapping};
#' \code{FALSE} otherwise.
#' 
#' @export
#' @family mapping object functions
#' @rdname is.mapping
is.mapping <- function(x) {
    
    status <- FALSE
    tryCatch({ status <- validateMapping(x)
    }, error=function(e) {})
    
    return(status)
}

# validateMapping --------------------------------------------------------------
#' Validate a \code{mapping} object.
#' 
#' @param x A \code{mapping} object.
#' 
#' @return \code{TRUE} if object is a valid \code{mapping};
#' otherwise, returns first error.
#' 
#' @keywords internal
validateMapping <- function(x) {
    
    stopifnot( 'mapping' %in% class(x) )
    
    if ( length(x) > 0 ) {
        
        if ( is.null( names(x) ) || any( is.na(names(x)) | names(x) == '' ) ) {
            stop("mapping must have keys")
        }
        
        dup.keys <- names(x)[ duplicated( names(x) ) ]
        if ( length(dup.keys) > 0 ) {
            stop("mapping keys must be unique")
        }
        
        if ( ! all( unlist( lapply(x, function(value)
            is.vector(value) && length(value) == 1 ) ) ) ) {
            stop("each mapping value must be a vector of length one")
        }
    }
    
    return(TRUE)
}

# `$<-.mapping` ----------------------------------------------------------------
#' @export
#' @keywords internal
`$<-.mapping` <- function(x, i, value) {
    x <- unclass(x)
    x[[i]] <- value
    class(x) <- c('mapping', 'list')
    validateMapping(x)
    return(x)
}

# `[<-.mapping` ----------------------------------------------------------------
#' @export
#' @keywords internal
`[<-.mapping` <- function(x, i, value) {
    x <- unclass(x)
    x[i] <- value
    class(x) <- c('mapping', 'list')
    validateMapping(x)
    return(x)
}

# `[[<-.mapping` ---------------------------------------------------------------
#' @export
#' @keywords internal
`[[<-.mapping` <- function(x, i, value) {
    x <- unclass(x)
    x[[i]] <- value
    class(x) <- c('mapping', 'list')
    validateMapping(x)
    return(x)
}

# `names<-.mapping` ------------------------------------------------------------
#' @export
#' @keywords internal
`names<-.mapping` <- function(x, value) {
    
    if ( ! is.null(value) ) {
        
        if ( ! is.character(value) ) {
            stop("mapping cannot have keys of class '",
                toString(class(value)), "'")
        }
        
        if ( length(value) != length(x) ) {
            stop("mapping cannot have ", length(value), " keys and ",
                length(x), " values")
        }
    }
    
    x <- unclass(x)
    names(x) <- value
    class(x) <- c('mapping', 'list')
    validateMapping(x)
    
    return(x)
}

# values -----------------------------------------------------------------------
#' Get values of a \code{mapping} object.
#' 
#' @param x A \code{mapping} object.
#' 
#' @return Unnamed list of values in a \code{mapping} object.
#' 
#' @export
#' @family mapping object functions
#' @keywords internal
#' @rdname values
values <- function(x) {
    stopifnot( 'mapping' %in% class(x) )
    return( unname( unclass(x) ) )
}

# End of mapping.R #############################################################