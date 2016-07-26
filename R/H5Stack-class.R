# Start of H5Stack-class.R #####################################################

# TODO: validate HDF5 names, HDF5 objects.

# H5Stack ----------------------------------------------------------------------
#' An S4 class to manage a HDF5 object stack.
#' 
#' @slot file File containing the HDF5 objects to be accessed.
#' @slot ids List of HDF5 object IDs.
#'  
#' @docType class
#' @keywords internal
#' @rdname H5Stack-class
H5Stack <- setClass('H5Stack', 
                       
    slots = c( 
        file='character',
        ids='list'
    ),
                       
    prototype=list( 
        file=character(),
        ids=list()
    )
)    

# H5Stack::initialize ----------------------------------------------------------
#' Init \code{H5Stack} object.
#' 
#' @param .Object A \code{H5Stack} object.
#' @param file A HDF5 file path.
#' @param h5name Name of HDF5 object to access.
#' 
#' @return A prepared \code{H5Stack} object.
#'
#' @aliases initialize,H5Stack-method
#' @docType methods
#' @keywords internal
#' @name initialize 
#' @rdname initialize-methods
setMethod('initialize', 'H5Stack', function(.Object, file, h5name=NULL) {
    .Object <- openStack(.Object, file, h5name)
    return(.Object)
})

# closeStack -------------------------------------------------------------------
#' Close all HDF5 objects in \code{H5Stack}.
#' 
#' @param h5stack A \code{H5Stack} object.
#' 
#' @return Empty \code{H5Stack} object.
#' 
#' @docType methods
#' @importFrom rhdf5 H5Fclose
#' @importFrom rhdf5 H5Gclose
#' @importFrom rhdf5 H5Dclose
#' @importFrom rhdf5 H5Iget_type
#' @keywords internal
#' @rdname closeStack-methods
setGeneric('closeStack', function(h5stack) { standardGeneric('closeStack') })

# H5Stack::closeStack ----------------------------------------------------------
#' @aliases closeStack,H5Stack-method
#' @rdname closeStack-methods
setMethod('closeStack', signature='H5Stack', definition = function(h5stack) {
    
    while ( length(h5stack) > 0 ) {
        
        h5obj <- peek(h5stack)
        
        h5type <- as.character( rhdf5::H5Iget_type(h5obj) )
        
        if ( h5type == 'H5I_GROUP' ) {
            rhdf5::H5Gclose(h5obj)
        } else if ( h5type == 'H5I_DATASET' ) {
            rhdf5::H5Dclose(h5obj)
        } else if ( h5type == 'H5I_FILE' ) {
            rhdf5::H5Fclose(h5obj)
        } else if ( ! h5type %in% 'H5I_BADID' ) {
            stop("unknown HDF5 object type - '", h5type, "'")
        }
        
        h5stack <- pop(h5stack)
    }

    return(h5stack)
})

# fileID -----------------------------------------------------------------------
#' Get HDF5 object ID of HDF5 file.
#' 
#' @param h5stack A \code{H5Stack} object.
#' 
#' @return HDF5 object ID of HDF5 file. Returns \code{NULL} if stack is empty.
#' 
#' @docType methods
#' @keywords internal
#' @rdname fileID-methods
setGeneric('fileID', function(h5stack) { standardGeneric('fileID') })

# H5Stack::fileID --------------------------------------------------------------
#' @aliases fileID,H5Stack-method
#' @rdname fileID-methods
setMethod('fileID', signature='H5Stack', definition = function(h5stack) { 
    return( if ( length(h5stack@ids) > 0 ) { h5stack@ids[[1]] } else { NULL } )
})

# openStack --------------------------------------------------------------------
#' Open stack of HDF5 objects.
#' 
#' @param h5stack A \code{H5Stack} object.
#' @param file A HDF5 file path.
#' @param h5name Name of HDF5 object to access.
#' 
#' @return Opened \code{H5Stack} object.
#' 
#' @docType methods
#' @importFrom rhdf5 H5Fopen
#' @importFrom rhdf5 H5Fcreate
#' @importFrom rhdf5 H5Lexists
#' @importFrom rhdf5 H5Gopen
#' @importFrom rhdf5 H5Gcreate
#' @keywords internal
#' @rdname openStack-methods
setGeneric('openStack', function(h5stack, file, h5name=NULL) { 
    standardGeneric('openStack') })

# H5Stack::openStack -----------------------------------------------------------
#' @aliases openStack,H5Stack-method
#' @rdname openStack-methods
setMethod('openStack', signature='H5Stack', definition = 
    function(h5stack, file, h5name=NULL) {
    
    if ( length(h5stack) > 0 ) {
        h5stack <- closeStack(h5stack)
    }
    
    stopifnot( isSingleString(file) )
    
    h5name <- resolveH5ObjectName(h5name)
    
    h5stack@file <- file
    
    if ( file.exists(file) ) {
        file.id <- rhdf5::H5Fopen(file)
    } else {
        file.id <- rhdf5::H5Fcreate(file)
    }

    h5stack <- push(h5stack, file.id)
    
    if ( h5name != '/' ) {
        
        group.names <- splitH5ObjectName(h5name)
    
        for ( group.name in group.names ) {
            
            prev.id <- peek(h5stack)
            
            if ( rhdf5::H5Lexists(prev.id, group.name) ) {
                curr.id <- rhdf5::H5Gopen(prev.id, group.name)
            } else {
                curr.id <- rhdf5::H5Gcreate(prev.id, group.name)
            }   
            
            h5stack <- push(h5stack, curr.id)
        }
    }
    
    return(h5stack)
})

# H5Stack::length --------------------------------------------------------------
#' Get number of HDF5 objects in \code{H5Stack}.
#' 
#' @param object \code{H5Stack} object.
#' 
#' @return Number of HDF5 objects in \code{H5Stack}.
#' 
#' @aliases show,H5Stack-method
#' @docType methods
#' @keywords internal 
#' @rdname show-methods
setMethod('length', signature='H5Stack', definition = function(x) { 
    return( length(x@ids) )
})

# peek -------------------------------------------------------------------------
#' Get HDF5 object ID at top of \code{H5Stack}.
#' 
#' @param h5stack A \code{H5Stack} object.
#' 
#' @return HDF5 object ID at top of \code{H5Stack}. Returns \code{NULL} if stack
#' is empty.
#' 
#' @docType methods
#' @keywords internal
#' @rdname peek-methods
setGeneric('peek', function(h5stack) { standardGeneric('peek') })

# H5Stack::peek ----------------------------------------------------------------
#' @aliases peek,H5Stack-method
#' @rdname peek-methods
setMethod('peek', signature='H5Stack', definition = function(h5stack) { 
    return( if ( length(h5stack@ids) > 0 ) { 
        h5stack@ids[[ length(h5stack) ]] } else { NULL } )
})

# pop --------------------------------------------------------------------------
#' Remove HDF5 object ID at top of \code{H5Stack}.
#' 
#' @param h5stack A \code{H5Stack} object.
#' 
#' @return Input \code{H5Stack} object with top element removed.
#' 
#' @docType methods
#' @keywords internal
#' @rdname pop-methods
setGeneric('pop', function(h5stack) { standardGeneric('pop') })

# H5Stack::pop -----------------------------------------------------------------
#' @aliases pop,H5Stack-method
#' @rdname pop-methods
setMethod('pop', signature='H5Stack', definition = function(h5stack) { 
    
    if ( length(h5stack) == 0 ) {
        stop("cannot pop H5Stack of length zero")
    }
    
    h5stack@ids <- h5stack@ids[ -length(h5stack) ]
    
    return(h5stack)
})

# push -------------------------------------------------------------------------
#' Push HDF5 object ID onto top of \code{H5Stack}.
#' 
#' @param h5stack A \code{H5Stack} object.
#' @param h5obj A HDF5 object ID.
#' 
#' @return Input \code{H5Stack} object with new element added.
#' 
#' @docType methods
#' @keywords internal
#' @rdname push-methods
setGeneric('push', function(h5stack, h5obj) { standardGeneric('push') })

# H5Stack::push ----------------------------------------------------------------
#' @aliases push,H5Stack-method
#' @rdname push-methods
setMethod('push', signature='H5Stack', definition = function(h5stack, h5obj) { 
    
    h5stack@ids <- c(h5stack@ids, h5obj)
    
    validObject(h5stack)
    
    return(h5stack)
})

# show.H5Stack -----------------------------------------------------------------
#' Display contents of \code{H5Stack} object.
#' 
#' @param object A \code{H5Stack} object.
#' 
#' @aliases show,H5Stack-method
#' @keywords internal
#' @rdname show-methods 
show.H5Stack <- function(object) { 
    return( summary(object) )
}

# summary.H5Stack --------------------------------------------------------------
#' Display summary of \code{H5Stack} object.
#' 
#' @param object A \code{H5Stack} object.
#' 
#' @aliases summary,H5Stack-method
#' @importFrom rhdf5 H5Iget_name
#' @keywords internal
#' @rdname summary-methods 
summary.H5Stack <- function(object) { 
    object.name <- rhdf5::H5Iget_name( object@ids[[ length(object) ]] )
    return( cat("HDF5 STACK\n\n        name ", object.name, 
        "\n    filename ", object@file, "\n") )
}

# H5Stack::setValidity ---------------------------------------------------------
#' Validate \code{H5Stack} object.
#' 
#' @param object A \code{H5Stack} object.
#' 
#' @return TRUE if object is valid; otherwise, a character vector of errors.
#'
#' @aliases setValidity,H5Stack-method
#' @docType methods
#' @importFrom rhdf5 H5Iget_type
#' @importFrom rhdf5 H5Iget_name
#' @keywords internal
#' @name setValidity
#' @rdname setValidity-methods
setValidity('H5Stack', function(object) { 
    
    if ( ! isSingleString(object@file) ) {
        return("H5Stack file must be a single string")
    }
    
    if ( length(object) > 0 ) {
        
        if ( ! file.exists(object@file) ) {
            return("H5Stack file must exist, if open")
        }
        
        types <- unlist( lapply(object@ids, typeof) )
        if ( ! all(types == 'S4') ) {
            stop("H5Stack IDs must be of type S4")
        }
        
        classes <- unlist( lapply(object@ids, class) )
        if ( ! all(classes == 'H5IdComponent') ) {
            stop("H5Stack IDs must be of class H5IdComponent")
        }
        
        h5types <- as.character( unlist( lapply(object@ids, rhdf5::H5Iget_type) ) )
        
        l <- length(object)
        k <- l - 1
        j <- 2
        i <- 1
        
        if ( ( l >= 1 && h5types[i] != 'H5I_FILE' ) || 
             ( l >= 2 && ! h5types[l] %in% c('H5I_GROUP', 'H5I_DATASET') ) ||
             ( l >= 3 && any(h5types[j:k] != 'H5I_GROUP') ) ) {
            stop("H5Stack has invalid H5 types")
        }
        
        h5names <- sapply(object@ids, rhdf5::H5Iget_name)
        
        if ( h5names[1] != '/' ) {
            stop("Root group ID must have name '/'")
        }
        
        if ( l >= 3 ) {
            
            group.names <- h5names[2:length(h5names)]
            
            split.names <- lapply(group.names, splitH5ObjectName)
            
            last.split.name <- split.names[[ length(split.names) ]]
            split.names <- split.names[ -length(split.names) ]
            
            for ( split.name in split.names ) {
                
                if ( any(split.name != last.split.name[getIndices(split.name)]) ) {
                    stop("H5Stack names must match")
                }
            }
        }
    }    
    
    return(TRUE)
})

# End of H5Stack-class.R #######################################################