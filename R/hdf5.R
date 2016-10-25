# Start of hdf5.R ##############################################################

# copyObjectsHDF5 --------------------------------------------------------------
#' Copy objects between HDF5 files.
#' 
#' @param infile Input HDF5 file.
#' @param outfile Output HDF5 file.
#' @param h5names Names of HDF5 objects to copy. If none are specified,
#' every HDF5 object in the input file is copied to the output file.
#' 
#' @importFrom methods new
#' @importFrom rhdf5 H5Adelete
#' @importFrom rhdf5 H5Dclose
#' @importFrom rhdf5 H5Dopen
#' @importFrom rhdf5 H5Gcreate
#' @importFrom rhdf5 H5Gopen
#' @importFrom rhdf5 H5Lexists
#' @importFrom rhdf5 h5ls
#' @importFrom rhdf5 h5read
#' @importFrom rhdf5 h5readAttributes
#' @importFrom rhdf5 h5writeAttribute
#' @importFrom rhdf5 h5writeDataset
#' @keywords internal
#' @rdname copyObjectsHDF5
copyObjectsHDF5 <- function(infile, outfile, h5names=NULL) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    stopifnot( isSingleString(outfile) )
    
    if ( ! is.null(h5names) ) {
        stopifnot( is.character(h5names) )
        stopifnot( length(h5names) > 0 )
        h5names <- sapply(unique(h5names), resolveH5ObjectName)
    }
    
    # Open root HDF5 object of input file.
    h5stack <- methods::new('H5Stack', infile)
    
    # Ensure last opened HDF5 stack will be closed on exit.
    on.exit( closeStack(h5stack) )
    
    # Get info on all objects in HDF5 file.
    h5info <- rhdf5::h5ls(peek(h5stack), datasetinfo=FALSE, recursive=TRUE)
    
    # Close input HDF5 file.
    h5stack <- closeStack(h5stack)
    
    # Get HDF5 name and type for every object in HDF5 file.
    h5info <- data.frame(
        h5name = sapply( getRowIndices(h5info), function(i)
            joinH5ObjectNameParts( c(h5info$group[i], h5info$name[i]) ) ),
        h5type = as.character(h5info$otype), stringsAsFactors=FALSE )
    
    # Add HDF5 name and type of root group.
    h5info <- rbind(h5info, c('/', 'H5I_GROUP'))
    
    err.h5types <- h5info$h5type[ ! h5info$h5type %in% c('H5I_GROUP', 'H5I_DATASET') ]
    if ( length(err.h5types) > 0 ) {
        stop("cannot copy HDF5 objects - unsupported types (",
            toString(err.h5types), ")")
    }
    
    err.h5names <- h5names[ ! h5names %in% h5info$h5name ]
    if ( length(err.h5names) > 0 ) {
        stop( "HDF5 objects not found - ",
            paste0("'", err.h5names, "'", collapse=', ') )
    }
    
    # If HDF5 object names given, keep only specified HDF5 object names.
    if ( ! is.null(h5names) ) {
        h5info <- h5info[h5info$h5name %in% h5names, ]
    }
    
    # Get sorted vector of HDF5 group names.
    # NB: sorting ensures that ancestral groups are copied before their descendants.
    group.h5names <- sort( h5info$h5name[ h5info$h5type == 'H5I_GROUP' ] )
    
    # Copy HDF5 groups from input to output file.
    for ( group.h5name in group.h5names ) {
        
        # Get attributes of this HDF5 group in input file.
        input.attrs <- rhdf5::h5readAttributes(infile, group.h5name)
        
        # Get attributes of this HDF5 group in output file, if present.
        if ( file.exists(outfile) && hasObjectHDF5(outfile, group.h5name) ) {
            output.attrs <- rhdf5::h5readAttributes(outfile, group.h5name)
        } else {
            output.attrs <- list()
        }
        
        # Open output HDF5 group.
        h5stack <- methods::new('H5Stack', outfile, group.h5name)
        
        # Get HDF5 group object.
        h5group <- peek(h5stack)
        
        # Remove existing attributes of output group.
        for ( k in names(output.attrs) ) {
            rhdf5::H5Adelete(h5group, k)
        }
        
        # Copy input group attributes to output group.
        for ( i in seq_along(input.attrs) ) {
            rhdf5::h5writeAttribute(h5obj=h5group, name=names(input.attrs)[i],
                attr=input.attrs[[i]])
        }
        
        # Close output HDF5 group.
        h5stack <- closeStack(h5stack)
    }
    
    # Get HDF5 dataset names.
    dataset.h5names <- h5info$h5name[ h5info$h5type == 'H5I_DATASET' ]
    
    # Copy HDF5 datasets from input to output file.
    for ( dataset.h5name in dataset.h5names ) {
        
        # Read dataset from input file.
        dataset <- rhdf5::h5read(file=infile, name=dataset.h5name,
            read.attributes=TRUE)
        
        # Open output HDF5 file.
        h5stack <- methods::new('H5Stack', outfile)
        
        # Get HDF5 root object.
        h5root <- peek(h5stack)
        
        # Write dataset to output file.
        rhdf5::h5writeDataset(h5loc=h5root, name=dataset.h5name, obj=dataset)
        
        # Copy input dataset attributes to output dataset.
        h5obj <- rhdf5::H5Dopen(h5root, dataset.h5name)
        h5writeAttributes(h5obj, object=dataset)
        rhdf5::H5Dclose(h5obj)
        
        # Close output HDF5 file.
        h5stack <- closeStack(h5stack)
    }
    
    return( invisible() )
}

# getMapNamesHDF5 --------------------------------------------------------------
#' Get names of maps in HDF5 file.
#' 
#' @param infile An input HDF5 file.
#' 
#' @return Vector of map names. Returns \code{NULL}
#' if there are no maps present.
#' 
#' @keywords internal
#' @rdname getMapNamesHDF5
getMapNamesHDF5 <- function(infile) {
    return( getObjectNamesHDF5(infile, 'Maps', relative=TRUE,
        min.depth=1, max.depth=1) )
}

# getObjectClassHDF5 -----------------------------------------------------------
#' Get R class of a HDF5 object.
#' 
#' If possible, the R class of the specified object is obtained from the
#' attributes \code{'R.class'} and \code{'class'}, in that order. If neither
#' attribute is available, the object is loaded from the HDF5 file and its R
#' class determined from the resulting R object.
#' 
#' @param infile An input HDF5 file.
#' @param h5name HDF5 object name.
#' 
#' @return R class of HDF5 object.
#' 
#' @importFrom rhdf5 h5read
#' @importFrom rhdf5 h5readAttributes
#' @keywords internal
#' @rdname getObjectClassHDF5
getObjectClassHDF5 <- function(infile, h5name) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    h5name <- resolveH5ObjectName(h5name)
    
    attrs <- rhdf5::h5readAttributes(infile, h5name)
    
    if ( 'R.class' %in% names(attrs) ) {
        
        cls <- attrs[['R.class']]
        
    } else if ( 'class' %in% names(attrs) ) {
        
        cls <- attrs[['class']]
        
    } else {
        
        obj <- rhdf5::h5read(infile, h5name, read.attributes=TRUE)
        cls <- class(obj)
    }
    
    return(cls)
}

# getObjectNamesHDF5 -----------------------------------------------------------
#' Get names of HDF5 objects from HDF5 file and name.
#' 
#' @param infile An input HDF5 file.
#' @param h5name Name of HDF5 object to use as a starting-point. By
#' default, results are returned with respect to the HDF5 root object.
#' @param relative Option indicating that the returned HDF5 object names should
#' be given relative to the starting-point. By default, absolute HDF5 object
#' names are returned.
#' @param min.depth Minimum recursion depth (relative to starting-point) from
#' which to return object names. If not specified, no minimum depth constraint
#' is applied, and the results include the starting-point HDF5 object (if all
#' other constraints are satisfied).
#' @param max.depth Maximum recursion depth (relative to starting-point) from
#' which to return object names. If not specified, no maximum depth constraint
#' is applied, and the results include every HDF5 object below the
#' starting-point (for which all other constraints are satisfied).
#' 
#' @return Character vector of names of HDF5 objects that lie within the section
#' of the HDF5 file defined by the starting-point HDF5 object and any specified
#' constraints.
#' 
#' @importFrom methods new
#' @importFrom rhdf5 H5Iget_type
#' @importFrom rhdf5 H5Lexists
#' @importFrom rhdf5 h5ls
#' @keywords internal
#' @rdname getObjectNamesHDF5
getObjectNamesHDF5 <- function(infile, h5name=NULL, relative=FALSE,
    min.depth=0, max.depth=NULL) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    stopifnot( isBOOL(relative) )
    stopifnot( isSingleNonNegativeWholeNumber(min.depth) )
    
    if ( ! is.null(max.depth) ) {
        stopifnot( isSingleNonNegativeWholeNumber(max.depth) )
    }
    
    h5name <- resolveH5ObjectName(h5name)
    
    # Open HDF5 file.
    h5stack <- methods::new('H5Stack', infile, h5name)
    on.exit( closeStack(h5stack) )
    
    # Check that starting-point HDF5 object exists.
    if( ! rhdf5::H5Lexists(fileID(h5stack), h5name) ) {
        stop("HDF5 object ('", h5name, "') not found in file - '",  infile, "'")
    }
    
    h5type <- as.character( rhdf5::H5Iget_type( peek(h5stack) ) )
    
    if ( ! h5type %in% c('H5I_GROUP', 'H5I_DATASET', 'H5I_FILE') ) {
        stop("unknown HDF5 object type - '", h5type, "'")
    }
    
    # Init vector of HDF5 object names.
    object.h5names <- character()
    
    # If starting-point is not a HDF5 dataset, and maximum depth is greater
    # than zero, get the names of HDF5 objects below the starting-point.
    if ( h5type != 'H5I_DATASET' && ( is.null(max.depth) || max.depth > 0 ) ) {
        
        # Get recursive listing if maximum depth is greater than one.
        if ( is.null(max.depth) || max.depth > 1 ) {
            recursive <- TRUE
        } else {
            recursive <- FALSE
        }
        
        # Get data-frame listing HDF5 objects below starting-point.
        h5info <- rhdf5::h5ls(peek(h5stack), datasetinfo=FALSE,
            recursive=recursive)
        
        # Get relative depths of HDF5 objects from number
        # of component separators in returned group name.
        object.depths <- unlist( lapply( strsplit(h5info$group, ''),
            function(x) length( grep('/', x, fixed=TRUE) ) ) )
        
        # If maximum depth not specified, get actual maximum depth.
        if ( is.null(max.depth) ) {
            max.depth <- max(object.depths)
        }
        
        # Append each HDF5 object satisfying the specified constraints.
        for ( i in getRowIndices(h5info) ) {
            
            if ( object.depths[i] >= min.depth && object.depths[i] <= max.depth ) {
                
                h5parts <- c(h5info$group[i], h5info$name[i])
                
                if (relative) {
                    object.h5name <- joinH5ObjectNameParts(h5parts, relative=TRUE)
                } else { # absolute
                    object.h5name <- joinH5ObjectNameParts( c(h5name, h5parts) )
                }
                
                object.h5names <- c(object.h5names, object.h5name)
            }
        }
    }
    
    # If minimum depth includes starting-point, add this to HDF5 object names.
    if ( min.depth == 0 ) {
        object.h5names <- c(h5name, object.h5names)
    }
    
    return(object.h5names)
}

# getPhenotypesHDF5 ------------------------------------------------------------
#' Get names of phenotypes in HDF5 file.
#' 
#' @param infile An input HDF5 file.
#' 
#' @return Vector of phenotype names. Returns \code{NULL}
#' if there are no phenotypes present.
#' 
#' @keywords internal
#' @rdname getPhenotypesHDF5
getPhenotypesHDF5 <- function(infile) {
    
    phenotypes <- NULL
    
    result.names <- getObjectNamesHDF5(infile, 'Results', relative=TRUE,
        min.depth=1, max.depth=1)
    
    if ( ! is.null(result.names) ) {
        
        phenotypes <- result.names[ result.names != 'Overview' ]
        
        if ( length(phenotypes) == 0 ) {
            phenotypes <- NULL
        }
    }
    
    return(phenotypes)
}

# getResultNamesHDF5 -----------------------------------------------------------
#' Get names of results for a given phenotype.
#' 
#' @param infile An input HDF5 file.
#' @param phenotype Name of phenotype (or equivalent analysis unit).
#' 
#' @return Vector of result names for the given phenotype. Returns \code{NULL}
#' if there are no results present for the given phenotype.
#' 
#' @keywords internal
#' @rdname getResultNamesHDF5
getResultNamesHDF5 <- function(infile, phenotype) {
    stopifnot( phenotype %in% getPhenotypesHDF5(infile) )
    h5name <- joinH5ObjectNameParts( c('Results', phenotype) )
    phenotype.results <- getObjectNamesHDF5(infile, h5name, relative=TRUE,
        min.depth=1, max.depth=1)
    return(phenotype.results)
}

# hasMapHDF5 -------------------------------------------------------------------
#' Test if HDF5 file contains the named map.
#' 
#' @param infile An input HDF5 file.
#' @param name Map name.
#' 
#' @return \code{TRUE} if the given HDF5 file contains the named map;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname hasMapHDF5
hasMapHDF5 <- function(infile, name) {
    stopifnot( isSingleString(name) )
    mapnames <- getMapNamesHDF5(infile)
    return( name %in% mapnames )
}

# hasObjectHDF5 ----------------------------------------------------------------
#' Test if HDF5 file contains the named object.
#' 
#' @param infile An input HDF5 file.
#' @param h5name HDF5 object name.
#' 
#' @return \code{TRUE} if the given HDF5 file contains the named object;
#' \code{FALSE} otherwise.
#' 
#' @importFrom methods new
#' @importFrom rhdf5 H5Lexists
#' @keywords internal
#' @rdname hasObjectHDF5
hasObjectHDF5 <- function(infile, h5name) {
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    h5name <- resolveH5ObjectName(h5name)
    h5stack <- methods::new('H5Stack', infile)
    on.exit( closeStack(h5stack) )
    return( rhdf5::H5Lexists(peek(h5stack), h5name) )
}

# h5writeAttributes ------------------------------------------------------------
#' Write attributes to HDF5 object.
#' 
#' @param h5obj HDF5 object to which attributes will be assigned.
#' @param object An R object whose attributes will be written.
#' @param attrs List of attributes to be written.
#' 
#' @importFrom rhdf5 h5writeAttribute
#' @importFrom rhdf5 h5writeAttribute.array
#' @importFrom rhdf5 h5writeAttribute.character
#' @importFrom rhdf5 h5writeAttribute.double
#' @importFrom rhdf5 h5writeAttribute.integer
#' @importFrom rhdf5 h5writeAttribute.logical
#' @importFrom rhdf5 h5writeAttribute.matrix
#' @keywords internal
#' @rdname h5writeAttributes
h5writeAttributes <- function(h5obj, object=NULL, attrs=NULL) {
    
    if ( ! is.null(object) && ! is.null(attrs) ) {
        stop("cannot specify both object and attribute list")
    } else if ( ! is.null(object) ) {
        attrs <- attributes(object)
    } else if ( ! is.list(attrs) ) {
        stop("must specify either object or attribute list")
    }
    
    if ( is.list(attrs) && length(attrs) > 0 ) {
        
        # If 'dimnames' attribute present, split into one attribute
        # per dimension, then store dimension names separately.
        if ( 'dimnames' %in% names(attrs) ) {
            
            stopifnot( 'dim' %in% names(attrs) )
            
            exp.keys <- paste0( 'dimnames.', c('names', seq_along(attrs[['dim']]) ) )
         
            if ( any( exp.keys %in% names(attrs) ) ) {
                stop("cannot split 'dimnames' without overwriting existing attributes")
            }
            
            object.dimnames <- attrs[['dimnames']]
            
            if ( ! is.null(object.dimnames) ) {
                
                if ( ! is.null( names(object.dimnames) ) ) {
                    attrs[['dimnames.names']] <- names(object.dimnames)
                }
                
                indices <- which( ! sapply(object.dimnames, is.null) )
                
                for ( i in indices ) {
                    key <- paste0('dimnames.', i)
                    attrs[[key]] <- object.dimnames[[i]]
                }
                
                attrs[['dimnames']] <- NULL
            }
        }
        
        # If attributes include row names, remove these.
        if ( 'row.names' %in% names(attrs) ) {
            attrs[['row.names']] <- NULL
        }
        
        # Write attributes.
        for ( i in seq_along(attrs) ) {
            rhdf5::h5writeAttribute(h5obj=h5obj, name=names(attrs)[i],
                attr=attrs[[i]])
        }
    }
    
    return( invisible() )
}

# joinH5ObjectNameParts --------------------------------------------------------
#' Join components of HDF5 object name.
#'  
#' @param components HDF5 name components.
#' @param relative Option indicating that the returned HDF5 object name should
#' be relative (i.e. should not contain a leading forward slash \code{'/'}).
#'     
#' @return HDF5 object name.
#' 
#' @keywords internal
#' @rdname joinH5ObjectNameParts
joinH5ObjectNameParts <- function(components, relative=FALSE) {
    
    stopifnot( isBOOL(relative) )
    
    if ( isSingleString(components) && components == '' ) {
        
        h5name <- '/'
        
    } else {
        
        if ( any( grepl('/', components) ) ) {
            components <- unlist( lapply(components, splitH5ObjectName) )   
        }
        
        # Join HDF5 name components.
        h5name <- paste(components, collapse='/')
        
        # Collapse multiple component separators.
        h5name <- gsub('/+', '/', h5name)
        
        # Strip leading component separator.
        h5name <- gsub('^/', '', h5name)
        
        # If not returning a relative HDF5 object name,
        # prepend leading component separator.
        if ( ! relative ) {
            h5name <- paste0('/', h5name)
        }
    }
    
    return(h5name)
}

# makeDefaultElementNames ------------------------------------------------------
#' Make default HDF5 group element names.
#' 
#' Following the convention of the \pkg{rhdf5} package, default HDF5 group
#' element names are of the form \code{'ELT1'}, where 1 can be any positive
#' integer. Leading zeros are used to pad element numbers, depending on the
#' group size, so that the generated name strings will sort in numerical order.
#'  
#' @param n Number of default group element names to generate.
#'      
#' @return Character vector of default group element names.
#' 
#' @keywords internal
#' @rdname makeDefaultElementNames
makeDefaultElementNames <- function(n) {
    stopifnot( isSinglePositiveWholeNumber(n) )
    return( sprintf(paste0("ELT%0", ceiling( log10(n + 1) ), "d"), 1:n) )
}

# makeDefaultMapName -----------------------------------------------------------
#' Get default map name.
#'   
#' @param map An \pkg{R/qtl} \code{map} object.
#'      
#' @return Default map name.
#' 
#' @keywords internal
#' @rdname makeDefaultMapName
makeDefaultMapName <- function(map) {
    stopifnot( 'map' %in% class(map) )
    default.name <- ifelse(isPhysicalMap(map), 'Physical Map', 'Genetic Map')
    return(default.name)
}

# makeGroupObjectNames ---------------------------------------------------------
#' Make valid HDF5 group object names.
#' 
#' Group names are made based on the proposed group names or the group size. 
#' Default group element names will replace every group name that is either 
#' an empty string or \code{NA} value. Default group element names are of the
#' form \code{'ELT1'}, where the number indicates the position of the element
#' within the given group.
#' 
#' @param group.names HDF5 group object names.
#' @param group.size Number of elements in HDF5 group.
#'     
#' @return HDF5 group object names.
#' 
#' @keywords internal
#' @rdname makeGroupObjectNames
#' @seealso makeDefaultElementNames
makeGroupObjectNames <- function(group.names=NULL, group.size=NULL) {
    
    if ( ! is.null(group.names) ) {
        
        stopifnot( is.character(group.names) )
        stopifnot( length(group.names) > 0 )
        
        if ( ! is.null(group.size) ) {
            stopifnot( group.size == length(group.names) )
        }
        
        group.indices <- seq_along(group.names)
        
        m <- regexec(const$pattern$h5element, group.names)
        matches <- regmatches(group.names, m)
        default.element.nums <- sapply(matches, function(x) as.integer(x[2]))
        mask <- ! is.na(default.element.nums)
        
        if( any(default.element.nums[mask] == group.indices[mask]) ) {
            stop("group names have invalid default element names")
        }
        
        group.names[mask] <- NA
        
        indices <- which( is.na(group.names) | ! nzchar(group.names) )
        
        stopif( anyDuplicated(group.names[indices]) )
        
    } else if ( ! is.null(group.size) ) {
        
        stopifnot( isSinglePositiveWholeNumber(group.size) )
        
        indices <- 1:group.size
        
    } else {
        
        stop("cannot validate group object names - not enough information")
    }
    
    if ( length(indices) > 0 ) {
        group.names[indices] <- makeDefaultElementNames(group.size)[indices]
    }
    
    return(group.names)
}

# readDatasetHDF5 --------------------------------------------------------------
#' Read HDF5 dataset.
#'    
#' The named HDF5 dataset is read from the specified HDF5 file using the
#' method for the dataset type, as determined from the \code{'R.class'}
#' or \code{'class'} attribute of the HDF5 object.
#'    
#' @param infile An input HDF5 file.
#' @param h5name HDF5 dataset name.
#' @param ... Further arguments (see below).
#' @param rownames.column For a \code{data.frame}, if any column has a name that
#' matches this parameter, that column is extracted from the \code{data.frame}
#' and the rownames are set from its contents.
#' 
#' @return R object corresponding to the named HDF5 object.
#' 
#' @export
#' @family HDF5 functions
#' @importFrom rhdf5 h5read
#' @keywords internal
#' @rdname readDatasetHDF5
readDatasetHDF5 <- function(infile, h5name, ...) {
    
    class.vector <- getObjectClassHDF5(infile, h5name)
    
    readDataset <- dispatchFromClassS3('readDatasetHDF5', class.vector, 'shmootl')
    
    dataset <- readDataset(infile, h5name)
    
    return(dataset)
}

# readDatasetHDF5.array --------------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.array <- function(infile, h5name, ...) {
    dataset <- readDatasetHDF5.default(infile, h5name)
    dataset <- aperm(dataset) # NB: R is column-major, HDF5 is row-major
    return(dataset)
}

# readDatasetHDF5.data.frame ---------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.data.frame <- function(infile, h5name, rownames.column='rownames') {
    
    dataset <- readDatasetHDF5.default(infile, h5name)
    
    stopifnot( 'R.colClasses' %in% names( attributes(dataset) ) )
    
    colClasses <- attr(dataset, 'R.colClasses')
    dataset <- coerceDataFrame(dataset, colClasses)
    attr(dataset, 'R.colClasses') <- NULL
    
    if ( rownames.column %in% colnames(dataset) && ! hasRownames(dataset) ) {
        dataset <- setRownamesFromColumn(dataset, col.name=rownames.column)
    } 
    
    return(dataset)
}

# readDatasetHDF5.default ------------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.default <- function(infile, h5name, ...) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    h5name <- resolveH5ObjectName(h5name)
    
    dataset <- rhdf5::h5read(infile, h5name, read.attributes=TRUE)
    
    attr.keys <- names( attributes(dataset) )
    
    dimnames.keys <- attr.keys[ grepl('^dimnames.([[:digit:]]+|names)$', attr.keys) ]
    
    dimnames.attrs <- attributes(dataset)[dimnames.keys]
    
    # Reassemble 'dimnames' attribute.
    if ( length(dimnames.attrs) > 0 ) {
        
        stopifnot( 'dim' %in% names( attributes(dataset) ) )
        
        dims <- attr(dataset, 'dim')
        
        dataset.dimnames <- vector( 'list', length(dims) )
        
        for ( i in seq_along(dims) ) {
            
            attr.key <- paste0('dimnames.', i)
            
            if ( attr.key %in% names(dimnames.attrs) ) {
                
                stopifnot( length(dimnames.attrs[[attr.key]]) == dims[i] )
                
                dataset.dimnames[[i]] <- dimnames.attrs[[attr.key]]
                
                attr(dataset, attr.key) <- NULL
                
            } else {
                
                dataset.dimnames[[i]] <- NULL
            }
        }
        
        if ( 'dimnames.names' %in% dimnames.keys ) {
            
            stopifnot( length(dimnames.attrs[['dimnames.names']]) == length(dims) )
            names(dataset.dimnames) <- dimnames.attrs[['dimnames.names']]
            attr(dataset, 'dimnames.names') <- NULL
        }
        
        attr(dataset, 'dimnames') <- dataset.dimnames
    }
    
    return(dataset)
}

# readDatasetHDF5.list ---------------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.list <- function(infile, h5name, ...) {
    
    # Get names from HDF5 attributes first, because once the list dataset has 
    # been loaded, the attribute names will be in alphabetical order.
    dataset.attrs <- readObjectAttributesHDF5(infile, h5name)
    original.names <- dataset.attrs[['names']]
    original.names[ original.names == 'NA' ] <- NA
    
    child.names <- getObjectNamesHDF5(infile, h5name, relative=TRUE,
        min.depth=1, max.depth=1)
    
    ordered.names <- makeGroupObjectNames(group.names=original.names, 
        group.size=length(child.names) )
    
    stopifnot( all( sort(child.names) == sort(ordered.names) ) )
    
    dataset <- vector( 'list', length(child.names) )
    names(dataset) <- ordered.names
    
    for ( child.name in child.names ) {
        
        child.h5name <- joinH5ObjectNameParts( c(h5name, child.name) )
        
        dataset[[child.name]] <- readDatasetHDF5(infile, child.h5name)
    }
    
    names(dataset) <- original.names
    
    attributes(dataset) <- dataset.attrs
    
    return(dataset)
}

# readDatasetHDF5.map ----------------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.map <- function(infile, h5name, ...) {
    
    h5attrs <- readObjectAttributesHDF5(infile, h5name)
    
    dataset <- readDatasetHDF5.data.frame(infile, h5name, rownames.column='id')
    
    stopifnot( 'R.class' %in% names(h5attrs) )
    stopifnot( 'map' %in% h5attrs[['R.class']] )
    special.attributes <- getSpecialAttributeNames(dataset)
    h5attrs <- h5attrs[ ! names(h5attrs) %in% special.attributes ]
    
    dataset <- as.map(dataset)
    
    for ( name in names(h5attrs) ) {
        attr(dataset, name) <- h5attrs[[name]]
    }
    
    class(dataset) <- attr(dataset, 'R.class')
    attr(dataset, 'R.class') <- NULL
    
    return(dataset)
}

# readDatasetHDF5.mapframe -----------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.mapframe <- function(infile, h5name, ...) {
    
    dataset <- readDatasetHDF5.data.frame(infile, h5name, rownames.column='id')
    
    stopifnot( 'R.class' %in% names( attributes(dataset) ) )
    
    class(dataset) <- attr(dataset, 'R.class')
    attr(dataset, 'R.class') <- NULL
    
    return(dataset)
}

# readDatasetHDF5.matrix -------------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.matrix <- function(infile, h5name, ...) {
    dataset <- readDatasetHDF5.default(infile, h5name)
    dataset <- t(dataset) # NB: R is column-major, HDF5 is row-major
    return(dataset)
}

# readDatasetHDF5.scanone ------------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.scanone <- function(infile, h5name, ...) {
    return( readDatasetHDF5.mapframe(infile, h5name) )
}

# readDatasetHDF5.scanonebins --------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.scanonebins  <- function(infile, h5name, ...) {
    
    dataset.attrs <- readObjectAttributesHDF5(infile, h5name)
    
    dataset <- readDatasetHDF5.data.frame(infile, h5name)
    
    attr.names <- names(dataset.attrs)
    
    stopifnot( 'R.class' %in% attr.names )
    
    trans.attrs <- dataset.attrs[ ! attr.names %in% c('class', 'names', 'row.names') ]
    
    perm.indices <- as.character(dataset$perm)
    
    dataset <- deleteColumn(dataset, col.name='perm')
    
    bin.labels <- colnames(dataset)
    
    dim.names <- list(perm.indices, bin.labels, 'lod')
    
    dataset <- array(unlist(dataset), dim=lengths(dim.names), dimnames=dim.names)
    
    for ( attr.name in names(trans.attrs) ) {
        attr(dataset, attr.name) <- trans.attrs[[attr.name]]
    }
    
    class(dataset) <- attr(dataset, 'R.class')
    attr(dataset, 'R.colClasses') <- NULL # currently unused
    attr(dataset, 'R.class') <- NULL
    
    return(dataset)
}

# readDatasetHDF5.scanoneperm --------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.scanoneperm  <- function(infile, h5name, ...) {
    
    dataset.attrs <- readObjectAttributesHDF5(infile, h5name)
    
    dataset <- readDatasetHDF5.data.frame(infile, h5name)
    
    attr.names <- names(dataset.attrs)
    
    stopifnot( 'R.class' %in% attr.names )
    
    trans.attrs <- dataset.attrs[ ! attr.names %in% c('class', 'names', 'row.names') ]
    
    perm.indices <- dataset$perm
    
    dataset <- matrix(dataset$lod, dimnames=list(perm.indices, 'lod') )
    
    for ( attr.name in names(trans.attrs) ) {
        attr(dataset, attr.name) <- trans.attrs[[attr.name]]
    }
    
    class(dataset) <- attr(dataset, 'R.class')
    attr(dataset, 'R.colClasses') <- NULL # currently unused
    attr(dataset, 'R.class') <- NULL
    
    return(dataset)
}

# readDatasetHDF5.scantwo ------------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.scantwo <- function(infile, h5name, ...) {
    
    dataset <- readDatasetHDF5.list(infile, h5name)
    
    stopifnot( 'R.class' %in% names( attributes(dataset) ) )
    
    dataset$map <- setRownamesFromColumn(dataset$map, col.name='id')
    dataset$map$eq.spacing <- as.numeric( dataset$map$eq.spacing )
    dataset$map$xchr <- FALSE
    
    attr(dataset, 'fullmap') <- as.map(dataset$fullmap)
    dataset[['fullmap']] <- NULL
    
    dataset['scanoneX'] <- list(NULL)
    
    class(dataset) <- attr(dataset, 'R.class')
    attr(dataset, 'R.class') <- NULL
 
    return(dataset)
}

# readDatasetHDF5.scantwoperm --------------------------------------------------
#' @export
#' @rdname readDatasetHDF5
readDatasetHDF5.scantwoperm <- function(infile, h5name, ...) {
    
    dataset.attrs <- readObjectAttributesHDF5(infile, h5name)
    
    dataset <- readDatasetHDF5.data.frame(infile, h5name)
    
    attr.names <- names(dataset.attrs)
    
    stopifnot( 'R.class' %in% attr.names )
    
    trans.attrs <- dataset.attrs[ ! attr.names %in% c('class', 'names', 'row.names') ]
    
    indices <- which( colnames(dataset) != 'perm' )
    
    lod.types <- colnames(dataset)[indices]
    
    dataset <- lapply(indices, function (i) matrix(dataset[[i]]))
    names(dataset) <- lod.types
    
    for ( attr.name in names(trans.attrs) ) {
        attr(dataset, attr.name) <- trans.attrs[[attr.name]]
    }
    
    class(dataset) <- attr(dataset, 'R.class')
    attr(dataset, 'R.colClasses') <- NULL # currently unused
    attr(dataset, 'R.class') <- NULL
    
    for ( i in seq_along(dataset) ) {
        
        colnames(dataset[[i]]) <- attr(dataset, 'phenotypes')
    }
    attr(dataset, 'phenotypes') <- NULL
    
    return(dataset)
}

# readDatasetHDF5.summary.scanonebins ------------------------------------------
#' @export readDatasetHDF5.summary.scanonebins
#' @method readDatasetHDF5 summary.scanonebins
#' @rdname readDatasetHDF5
readDatasetHDF5.summary.scanonebins <- function(infile, h5name, ...) {
    
    dataset <- readDatasetHDF5.data.frame(infile, h5name)
    
    dataset.attrs <- attributes(dataset)
    
    attr.names <- names( attributes(dataset) )
    
    stopifnot( 'R.class' %in% attr.names )
    
    trans.attrs <- dataset.attrs[ ! attr.names %in% c('class', 'names', 'row.names') ]
    
    dataset[, 'FDR'] <- paste0(as.character(100 * dataset[, 'FDR']), '%')
    
    dataset <- setRownamesFromColumn(dataset, col.name='FDR')
    
    dataset <- data.matrix(dataset)
    
    for ( attr.name in names(trans.attrs) ) {
        attr(dataset, attr.name) <- trans.attrs[[attr.name]]
    }
    
    class(dataset) <- attr(dataset, 'R.class')
    attr(dataset, 'R.class') <- NULL
    
    return(dataset)
}

# readDatasetHDF5.summary.scanoneperm ------------------------------------------
#' @export readDatasetHDF5.summary.scanoneperm
#' @method readDatasetHDF5 summary.scanoneperm
#' @rdname readDatasetHDF5
readDatasetHDF5.summary.scanoneperm <- function(infile, h5name, ...) {
    
    dataset <- readDatasetHDF5.data.frame(infile, h5name)
    
    dataset.attrs <- attributes(dataset)
    
    attr.names <- names( attributes(dataset) )
    
    stopifnot( 'R.class' %in% attr.names )
    
    trans.attrs <- dataset.attrs[ ! attr.names %in% c('class', 'names', 'row.names') ]
    
    dataset[, 'alpha'] <- paste0(as.character(100 * dataset[, 'alpha']), '%')
    
    dataset <- setRownamesFromColumn(dataset, col.name='alpha')
    
    dataset <- data.matrix(dataset)
    
    for ( attr.name in names(trans.attrs) ) {
        attr(dataset, attr.name) <- trans.attrs[[attr.name]]
    }
    
    class(dataset) <- attr(dataset, 'R.class')
    attr(dataset, 'R.class') <- NULL
    
    return(dataset)
}

# readMapHDF5 ------------------------------------------------------------------
#' Read \code{map} from a HDF5 file.
#'   
#' @details A \code{mapframe} dataset is read from the \code{'Maps'} group at
#' the root of the HDF5 file, and converted to an \pkg{R/qtl} \code{map} object.
#' 
#' @param infile An input HDF5 file.
#' @param name Map name.
#'      
#' @return An \pkg{R/qtl} \code{map} object. Returns \code{NULL}
#' if no map is found with the specified name.
#' 
#' @export
#' @family HDF5 functions
#' @family map utility functions
#' @keywords internal
#' @rdname readMapHDF5
readMapHDF5 <- function(infile, name) {
    
    stopifnot( isValidID(name) )
    
    map.names <- getMapNamesHDF5(infile)
    
    if ( name %in% map.names ) {
        h5name <- joinH5ObjectNameParts( c('Maps', name) )
        result <- readDatasetHDF5(infile, h5name)
    } else {
        result <- NULL
    }
    
    return(result)
}

# readObjectAttributesHDF5 -----------------------------------------------------
#' Read HDF5 object attributes.
#'   
#' @param infile An input HDF5 file.
#' @param h5name HDF5 object name.
#' 
#' @return List of attributes of the named HDF5 object.
#'  
#' @importFrom rhdf5 h5readAttributes
#' @keywords internal
#' @rdname readObjectAttributesHDF5
readObjectAttributesHDF5 <- function(infile, h5name) {
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    h5name <- resolveH5ObjectName(h5name)
    return( rhdf5::h5readAttributes(infile, h5name) )
}

# readOverviewHDF5 -------------------------------------------------------------
#' Read QTL analysis results overview from a HDF5 file.
#'    
#' @param infile An input HDF5 file.
#' 
#' @return A \code{data.frame} containing a QTL analysis results overview.
#' 
#' @export
#' @family HDF5 functions
#' @keywords internal
#' @rdname readOverviewHDF5
readOverviewHDF5 <- function(infile) {
    
    h5name <- joinH5ObjectNameParts( c('Results', 'Overview') )
    overview <- readDatasetHDF5(infile, h5name)
    
    stopifnot( is.data.frame(overview) )
    stopifnot( all( colnames(overview) == c('Phenotype', 'Status', 'Comments') ) )
    stopifnot( is.character(overview[['Phenotype']]) )
    stopifnot( is.logical(overview[['Status']]) )
    stopifnot( is.character(overview[['Comments']]) )
    
    return(overview)
}

# readResultHDF5 ---------------------------------------------------------------
#' Read QTL analysis result from a HDF5 file.
#' 
#' @param infile An input HDF5 file.
#' @param phenotype Name of phenotype (or equivalent analysis unit).
#' @param name Name of QTL analysis result.
#' 
#' @return R object containing QTL analysis result. Returns \code{NULL}
#' if the given phenotype has no result with the specified name.
#' 
#' @export
#' @family HDF5 functions
#' @keywords internal
#' @rdname readResultHDF5
readResultHDF5 <- function(infile, phenotype, name) {
    
    stopifnot( isValidID(phenotype) )
    stopifnot( isValidID(name) )
    
    phenotype.results <- getResultNamesHDF5(infile, phenotype)
    
    if ( name %in% phenotype.results ) {
        h5name <- joinH5ObjectNameParts( c('Results', phenotype, name) )
        result <- readDatasetHDF5(infile, h5name)
    } else {
        result <- NULL
    }
    
    return(result)
}

# resolveH5ObjectName ----------------------------------------------------------
#' Resolve absolute HDF5 object name.
#' 
#' Validates and (if necessary) corrects a HDF5 object name: ensuring
#' that it has a leading slash and does not have a trailing slash.
#' 
#' @param components HDF5 name components.
#'     
#' @return HDF5 object name.
#' 
#' @keywords internal
#' @rdname resolveH5ObjectName
resolveH5ObjectName <- function(h5name) {
    
    if ( ! is.null(h5name) ) {
        
        stopifnot( isSingleString(h5name) )
        
        res <- sub('^/?', '/', unname(h5name))
        
        if ( res != '/' ) {
            res <- sub('/$', '', res)
        }
        
    } else {
        
        res <- '/'
    }
    
    return(res)
}

# resolveMapNameHDF5 -----------------------------------------------------------
#' Resolve map name in HDF5 file.
#' 
#' This function resolves the name of a map in the given HDF5 file: if the given
#' map is present in the HDF5 file, that map name is returned. If no map name is
#' specified and one map is found, the name of that map is returned by default;
#' otherwise a map name must be given.
#' 
#' @param infile Input HDF5 file.
#' @param name Map name.
#' 
#' @return Resolved map name.
#' 
#' @keywords internal
#' @rdname resolveMapNameHDF5
resolveMapNameHDF5 <- function(infile, name=NULL) {
    
    mapnames <- getMapNamesHDF5(infile)
    
    if ( ! is.null(name) ) {
        
        if ( ! name %in% mapnames ) {
            stop("map ('", name, "') not found in HDF5 file - '", infile, "'")
        }
        
    } else {
        
        if ( length(mapnames) == 1 ) {
            name <- mapnames[1]
        } else if ( length(mapnames) > 1 ) {
            stop("cannot resolve multiple maps from HDF5 file - please choose one")
        } else if ( length(mapnames) == 0 ) {
            stop("no maps found in HDF5 file - '", infile, "'")
        }
    }
    
    return(name)
}

# splitH5ObjectName ------------------------------------------------------------
#' Split HDF5 object name into component parts.
#' 
#' @param h5name Character vector representing a HDF5 object name. If this has
#' multiple parts, only the initial part can contain an absolute H5Object name.
#' 
#' @return HDF5 name components.
#' 
#' @keywords internal
#' @rdname splitH5ObjectName
splitH5ObjectName <- function(h5name) {
    
    stopifnot( is.character(h5name) )
    stopifnot( length(h5name) > 0 )
    
    if ( length(h5name) == 1 ) {
        
        if ( h5name == '/' ) {
            
            components <- ''
            
        } else {
            
            trimmed.h5name <- gsub('(^/)|(/$)', '', h5name)
            
            components <- unname( unlist( strsplit(trimmed.h5name, '/') ) )
            
            if ( ! all( isValidID(components) ) ) {
                stop("invalid H5Object name - '", toString(h5name), "'")
            }
        }
        
    } else if ( length(h5name) > 1 ) {
        
        if ( any( substring(h5name[-1], 1, 1) == '/' ) ) {
            stop("invalid multi-part H5Object name - '", toString(h5name), "'")
        }
        
        components <- unlist( lapply(h5name, splitH5ObjectName) )
    }
    
    return(components)
}

# writeDatasetHDF5 -------------------------------------------------------------
#' Write R object to a HDF5 dataset.
#'   
#' The given R object is written to the named HDF5 dataset using the method
#' for the given dataset type, as determined by the class of the R object.
#'  
#' @param dataset R object.
#' @param outfile An output HDF5 file.
#' @param h5name HDF5 dataset name.
#' @param ... Further arguments (see below).
#' @param rownames.column For a \code{data.frame}, any rownames are removed,
#' and then inserted into a column of the main table, with its column name
#' taken from this parameter.
#' @param overwrite Option indicating whether to force overwrite of an
#' existing dataset. By default, existing datasets cannot be overwritten.
#' 
#' @export
#' @family HDF5 functions
#' @importFrom methods new
#' @importFrom rhdf5 H5Dopen
#' @importFrom rhdf5 h5writeDataset
#' @importFrom rhdf5 h5writeDataset.array
#' @importFrom rhdf5 h5writeDataset.character
#' @importFrom rhdf5 h5writeDataset.data.frame
#' @importFrom rhdf5 h5writeDataset.double
#' @importFrom rhdf5 h5writeDataset.integer
#' @importFrom rhdf5 h5writeDataset.list
#' @importFrom rhdf5 h5writeDataset.logical
#' @importFrom rhdf5 h5writeDataset.matrix
#' @rdname writeDatasetHDF5
writeDatasetHDF5 <- function(dataset, outfile, h5name, ...) {
    UseMethod('writeDatasetHDF5', dataset)
}

# writeDatasetHDF5.array -------------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.array <- function(dataset, outfile, h5name, ...) {
    dataset <- aperm(dataset) # NB: R is column-major, HDF5 is row-major
    writeDatasetHDF5.default(dataset, outfile, h5name, ...)
    return( invisible() )
}

# writeDatasetHDF5.data.frame --------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.data.frame <- function(dataset, outfile, h5name,
    rownames.column='rownames', ...) {

    stopif( 'R.colClasses' %in% names( attributes(dataset) ) )
    
    stopifnot( nrow(dataset) > 0 )
    stopifnot( ncol(dataset) > 0 )
    
    if ( hasRownames(dataset) ) {
        dataset <- setColumnFromRownames(dataset, col.name=rownames.column)
    }
    
    colClasses <- sapply(dataset, class)
    cols <- which( colClasses %in% c('factor', 'logical') )
    
    for ( col in cols ) {
        dataset[[col]] <- as.character(dataset[[col]])
    }
    
    attr(dataset, 'R.colClasses') <- colClasses
    
    writeDatasetHDF5.default(dataset, outfile, h5name, ...)
    
    return( invisible() )
}

# writeDatasetHDF5.default -----------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.default <- function(dataset, outfile, h5name,
    overwrite=FALSE, ...) {

    stopifnot( isSingleString(outfile) )
    stopifnot( isBOOL(overwrite) )
    h5name <- resolveH5ObjectName(h5name)
    
    # If output file exists and overwrite option is TRUE,
    # prepare to overwrite the dataset via a temp file..
    if ( file.exists(outfile) && hasObjectHDF5(outfile, h5name) && overwrite ) {
        sinkfile <- tempfile()
        overwriting <- TRUE
    } else { # ..otherwise prepare to write directly to output file.
        sinkfile <- outfile
        overwriting <- FALSE
    }
    
    # Split dataset HDF5 object name into the name of the group
    # containing the dataset, and the name of the dataset itself.
    name.parts <- splitH5ObjectName(h5name)
    h5loc.name <- joinH5ObjectNameParts(name.parts[ -length(name.parts) ])
    dataset.name <- name.parts[ length(name.parts) ]
    
    # Open new HDF5 stack.
    h5stack <- methods::new('H5Stack', sinkfile, h5loc.name)
    on.exit({ closeStack(h5stack); if (overwriting) { file.remove(sinkfile) } })
    
    # Write dataset to HDF5.
    rhdf5::h5writeDataset(h5loc=peek(h5stack), name=dataset.name, obj=dataset)
    
    # Write dataset attributes to HDF5.
    h5stack <- push(h5stack, rhdf5::H5Dopen(peek(h5stack), dataset.name))
    h5writeAttributes(peek(h5stack), object=dataset)
    
    h5stack <- closeStack(h5stack)
    
    # If overwriting dataset via a temp file, transfer any remaining objects to
    # the temp file, then copy the complete temp file to the final output file.
    if (overwriting) {
        
        # Transfer all existing but unchanged objects
        # from existing output file to temp file.
        updated <- getObjectNamesHDF5(sinkfile, h5name)
        existing <- getObjectNamesHDF5(outfile)
        unchanged <- setdiff(existing, updated) # NB: guarantees unique
        copyObjectsHDF5(outfile, sinkfile, h5names=unchanged)
        
        # Move temp file to final output file.
        # NB: file.copy is used here instead of file.rename because the latter
        # can sometimes fail when moving files between different file systems.
        file.copy(sinkfile, outfile, overwrite=TRUE)
    }
    
    return( invisible() )
}

# writeDatasetHDF5.list --------------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.list <- function(dataset, outfile, h5name, ...) {
    
    writeObjectAttributesHDF5(outfile, h5name, object=dataset)
    
    child.names <- makeGroupObjectNames( group.names=names(dataset), 
        group.size=length(dataset) )
    
    for ( i in seq_along(dataset) ) {
        
        child.h5name <- joinH5ObjectNameParts( c(h5name, child.names[i]) )
        
        writeDatasetHDF5(dataset[[i]], outfile, child.h5name, ...)
    }
    
    return( invisible() )
}

# writeDatasetHDF5.map ---------------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.map <- function(dataset, outfile, h5name, ...) {
    
    stopif( 'R.class' %in% names( attributes(dataset) ) )
    
    dataset.Rclass <- class(dataset)
    
    dataset <- as.data.frame(dataset)
    
    attr(dataset, 'R.class') <- dataset.Rclass
    
    writeDatasetHDF5.data.frame(dataset, outfile, h5name,
        rownames.column='id', ...)
    
    return( invisible() )
}

# writeDatasetHDF5.mapframe ----------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.mapframe <- function(dataset, outfile, h5name, ...) {
    
    stopif( 'R.class' %in% names( attributes(dataset) ) )
    
    attr(dataset, 'R.class') <- class(dataset)
    
    dataset <- as.data.frame(dataset)
    
    writeDatasetHDF5.data.frame(dataset, outfile, h5name,
        rownames.column='id', ...)

    return( invisible() )
}

# writeDatasetHDF5.matrix ------------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.matrix <- function(dataset, outfile, h5name, ...) {
    dataset <- t(dataset) # NB: R is column-major, HDF5 is row-major
    writeDatasetHDF5.default(dataset, outfile, h5name, ...)
    return( invisible() )
}

# writeDatasetHDF5.scanone -----------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.scanone <- function(dataset, outfile, h5name, ...) {
    
    num.phenotypes <- length( getLodColIndices(dataset) )
    stopifnot( num.phenotypes == 1 )
    
    writeDatasetHDF5.mapframe(dataset, outfile, h5name, ...)
    
    return( invisible() )
}

# writeDatasetHDF5.scanonebins -------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.scanonebins <- function(dataset, outfile, h5name, ...) {
    
    stopif( 'R.class' %in% names( attributes(dataset) ) )
    
    num.phenotypes <- dim(dataset)[3]
    stopifnot( num.phenotypes == 1 )
    
    others <- otherattributes(dataset)
    
    bin.labels <- dimnames(dataset)[[2]]
    
    perm.indices <- as.integer( dimnames(dataset)[[1]] )
    
    num.bins <- dim(dataset)[2]
    
    dataset <- sapply(1:num.bins, function(i) as.vector(dataset[, i, 1]))
    
    dataset <- data.frame(dataset, check.names=FALSE)
    colnames(dataset) <- bin.labels
    
    dataset <- insertColumn(dataset, col.index=1,
        col.name='perm', data=perm.indices)
    
    otherattributes(dataset) <- others
    
    attr(dataset, 'R.class') <- c('scanonebins', 'array')
    
    writeDatasetHDF5.data.frame(dataset, outfile, h5name, ...)
    
    return( invisible() )
}

# writeDatasetHDF5.scanoneperm -------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.scanoneperm <- function(dataset, outfile, h5name, ...) {
    
    stopif( 'R.class' %in% names( attributes(dataset) ) )
    
    num.phenotypes <- ncol(dataset)
    stopifnot( num.phenotypes == 1 )
    
    dataset.attrs <- attributes(dataset)
    
    attr.keys <- names(dataset.attrs)
    
    trans.attrs <- dataset.attrs[ ! attr.keys %in% c('class', 'dim', 'dimnames') ]
    
    perm.indices <- as.integer( rownames(dataset) )
    class(dataset) <- c('scanoneperm', 'matrix')
    
    scanoneperm.dataset <- data.frame(perm=perm.indices, lod=unname(dataset))
    
    for ( attr.key in names(trans.attrs) ) {
        attr(scanoneperm.dataset, attr.key) <- trans.attrs[[attr.key]]
    }
    
    attr(scanoneperm.dataset, 'R.class') <- 'scanoneperm'
    class(scanoneperm.dataset) <- 'data.frame'
    
    writeDatasetHDF5.data.frame(scanoneperm.dataset, outfile, h5name, ...)
    
    return( invisible() )
}

# writeDatasetHDF5.scantwo -----------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.scantwo <- function(dataset, outfile, h5name, ...) {
    
    stopif( 'R.class' %in% names( attributes(dataset) ) )
    
    stopifnot( is.null(dataset$scanoneX) )
    stopif( any( isTRUE(dataset$map$xchr) ) )
    
    pheno.names <- attr(dataset, 'phenotypes')
    stopifnot( length(pheno.names) == 1 )
    
    dataset$scanoneX <- NULL
    attr(dataset$lod, 'class') <- 'matrix'
    
    scantwo.map <- dataset[['map']]
    scantwo.map <- setColumnFromRownames(scantwo.map, col.name='id')
    scantwo.map$eq.spacing <- as.logical(scantwo.map$eq.spacing)
    scantwo.map[['xchr']] <- NULL
    dataset[['map']] <- scantwo.map
    
    dataset[['fullmap']] <- as.mapframe(attr(dataset, 'fullmap'), map.unit='cM')
    attr(dataset, 'fullmap') <- NULL
    
    attr(dataset, 'R.class') <- class(dataset)
    class(dataset) <- 'list'
    
    writeDatasetHDF5.list(dataset, outfile, h5name, ...)
    
    return( invisible() )
}

# writeDatasetHDF5.scantwoperm -------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.scantwoperm <- function(dataset, outfile, h5name, ...) {
    
    stopif( 'R.class' %in% names( attributes(dataset) ) )
    
    phenames.vectors <- lapply(unname(dataset), colnames)
    stopifnot( length( unique(phenames.vectors) ) == 1 )
    pheno.names <- phenames.vectors[[1]]
    stopifnot( length(pheno.names) == 1 )
    
    dataset.attrs <- attributes(dataset)
    
    attr.keys <- names(dataset.attrs)
    
    trans.attrs <- dataset.attrs[ ! attr.keys %in% c('class', 'names') ]
    
    # The scantwoperm object should be a list of six matrices, each 
    # with one column of length equal to number of permutations.
    lod.types <- c('full', 'fv1', 'int', 'add', 'av1', 'one')
    
    num.perms <- unique( sapply(dataset, nrow) )
    
    if (  length(num.perms) != 1 ) {
        stop("inconsistent number of permutations in 'scantwoperm' object")
    }
    
    perm.indices <- 1:num.perms
    
    lod.types <- names(dataset)
    
    scantwoperm.dataset <- data.frame(perm=perm.indices)
    
    for ( lod.type in lod.types ) {
        scantwoperm.dataset[[lod.type]] <- dataset[[lod.type]][, 1]
    }
    
    for ( attr.name in names(trans.attrs) ) {
        attr(scantwoperm.dataset, attr.name) <- attr(dataset, attr.name)
    }
    
    attr(scantwoperm.dataset, 'phenotypes') <- pheno.names
    
    attr(scantwoperm.dataset, 'R.class') <- c('scantwoperm', 'list')
    class(scantwoperm.dataset) <- 'data.frame'
    
    writeDatasetHDF5.data.frame(scantwoperm.dataset, outfile, h5name, ...)
    
    return( invisible() )
}

# writeDatasetHDF5.summary.scanonebins -----------------------------------------
#' @export
#' @method writeDatasetHDF5 summary.scanonebins
#' @rdname writeDatasetHDF5
writeDatasetHDF5.summary.scanonebins <- function(dataset, outfile, h5name, ...) {
    
    stopif( 'R.class' %in% names( attributes(dataset) ) )
    num.lodcolumns <- ncol(dataset)
    stopifnot( num.lodcolumns == 1 )
    
    attr(dataset, 'R.class') <- class(dataset)
    
    dataset.attrs <- attributes(dataset)
    
    attr.keys <- names(dataset.attrs)
    
    trans.attrs <- dataset.attrs[ ! attr.keys %in% c('class', 'dim', 'dimnames') ]
    
    dataset <- data.frame(unclass(dataset), check.names=FALSE)
    
    dataset <- setColumnFromRownames(dataset, col.name='FDR')
    
    dataset[, 'FDR'] <- 0.01 * as.numeric( sub('%', '', dataset[, 'FDR']) )
    
    for ( attr.name in names(trans.attrs) ) {
        attr(dataset, attr.name) <- trans.attrs[[attr.name]]
    }
    
    writeDatasetHDF5.data.frame(dataset, outfile, h5name, ...)
    
    return( invisible() )
}

# writeDatasetHDF5.summary.scanoneperm -----------------------------------------
#' @export
#' @method writeDatasetHDF5 summary.scanoneperm
#' @rdname writeDatasetHDF5
writeDatasetHDF5.summary.scanoneperm <- function(dataset, outfile, h5name, ...) {
    
    stopif( 'R.class' %in% names( attributes(dataset) ) )
    num.lodcolumns <- ncol(dataset)
    stopifnot( num.lodcolumns == 1 )
    
    attr(dataset, 'R.class') <- class(dataset)
    
    dataset.attrs <- attributes(dataset)
    
    attr.keys <- names(dataset.attrs)
    
    trans.attrs <- dataset.attrs[ ! attr.keys %in% c('class', 'dim', 'dimnames') ]
    
    dataset <- data.frame(unclass(dataset), check.names=FALSE)
    
    dataset <- setColumnFromRownames(dataset, col.name='alpha')
    
    dataset[, 'alpha'] <- 0.01 * as.numeric( sub('%', '', dataset[, 'alpha']) )
    
    for ( attr.name in names(trans.attrs) ) {
        attr(dataset, attr.name) <- trans.attrs[[attr.name]]
    }
    
    writeDatasetHDF5.data.frame(dataset, outfile, h5name, ...)
    
    return( invisible() )
}

# writeMapHDF5 -----------------------------------------------------------------
#' Write \code{map} to a HDF5 file.
#' 
#' @details An \pkg{R/qtl} \code{map} object is written to the \code{'Maps'}
#' group at the root of the given HDF5 file. If no map name is specified, a
#' default will be assigned based on whether it is a genetic or physical map.
#' 
#' @param map An \pkg{R/qtl} \code{map} object.
#' @param outfile An output HDF5 file.
#' @param name Map name.
#' @param overwrite Option indicating whether to force overwrite of an
#' existing map. By default, existing maps cannot be overwritten.
#' 
#' @export
#' @family HDF5 functions
#' @family map utility functions
#' @keywords internal
#' @rdname writeMapHDF5
writeMapHDF5 <- function(map, outfile, name=NULL, overwrite=FALSE) {
    
    stopifnot( 'map' %in% class(map) )
    
    if ( ! is.null(name) ) {
        stopifnot( isValidID(name) )
    } else {
        name <- makeDefaultMapName(map)
    }
    
    h5name <- joinH5ObjectNameParts( c('Maps', name) )
    writeDatasetHDF5(map, outfile, h5name, overwrite=overwrite)
    
    return( invisible() )
}

# writeObjectAttributesHDF5 ----------------------------------------------------
#' Write attributes to a HDF5 object.
#' 
#' @param outfile An output HDF5 file.
#' @param h5name HDF5 object name.
#' @param object An R object whose attributes will be written.
#' @param attrs List of attributes to be written.
#' 
#' @importFrom methods new
#' @keywords internal
#' @rdname writeObjectAttributesHDF5
writeObjectAttributesHDF5 <- function(outfile, h5name, object=NULL, attrs=NULL) {
    
    stopifnot( isSingleString(outfile) )
    
    h5stack <- methods::new('H5Stack', outfile, h5name)
    
    on.exit( closeStack(h5stack) )
    
    h5writeAttributes(peek(h5stack), object=object, attrs=attrs)
    
    return( invisible() )
}

# writeOverviewHDF5 ------------------------------------------------------------
#' Write QTL analysis results overview to a HDF5 file.
#'    
#' @param overview Object containing a QTL analysis results overview.
#' @param outfile An output HDF5 file.
#' @param overwrite Option indicating whether to force overwrite of an
#' existing overview. By default, existing overviews cannot be overwritten.
#' 
#' @export
#' @family HDF5 functions
#' @rdname writeOverviewHDF5
writeOverviewHDF5 <- function(overview, outfile, overwrite=FALSE) {
    
    # TODO: review result overview format.
    
    stopifnot( is.data.frame(overview) )
    stopifnot( all( colnames(overview) == c('Phenotype', 'Status', 'Comments') ) )
    stopifnot( is.character(overview[['Phenotype']]) )
    stopifnot( is.logical(overview[['Status']]) )
    stopifnot( is.character(overview[['Comments']]) )
    
    h5name <- joinH5ObjectNameParts( c('Results', 'Overview') )
    writeDatasetHDF5(overview, outfile, h5name, overwrite=overwrite)
    return( invisible() )
}

# writeResultHDF5 --------------------------------------------------------------
#' Write QTL analysis result to a HDF5 file.
#' 
#' @param result R object containing a QTL analysis result.
#' @param outfile An output HDF5 file.
#' @param phenotype Name of phenotype (or equivalent analysis unit).
#' @param name Name of QTL analysis result.
#' @param overwrite Option indicating whether to force overwrite of an
#' existing result. By default, existing results cannot be overwritten.
#' 
#' @export
#' @family HDF5 functions
#' @rdname writeResultHDF5
writeResultHDF5 <- function(result, outfile, phenotype, name, overwrite=FALSE) {
    stopifnot( isValidID(phenotype) )
    stopifnot( isValidID(name) )
    h5name <- joinH5ObjectNameParts( c('Results', phenotype, name) )
    writeDatasetHDF5(result, outfile, h5name, overwrite=overwrite)
    return( invisible() )
}

# End of hdf5.R ################################################################
