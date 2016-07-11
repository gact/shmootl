# Start of hdf5.R ##############################################################

# HDF5 Utilities ---------------------------------------------------------------
#' HDF5 input/output utilities.
#' 
#' @description Functions for HDF5 input and output. These should only be
#' used internally. External HDF5 input/output should be done through an
#' exported package function (e.g. \code{\link{run_scanone}}).
#' 
#' @details HDF5 files are assumed to have one or more of the following root 
#' groups.
#' 
#' \subsection{Maps}{ This group contains genetic or physical maps, which can 
#' be read with \code{\link{readMapHDF5}}, or written with 
#' \code{\link{writeMapHDF5}}.
#' }
#'   
#' \subsection{Results}{ Each subgroup contains results for a single phenotype
#' (or equivalent analysis unit). Phenotype names must be unique within each
#' file. Results can be read with \code{\link{readResultHDF5}}, or written with
#' \code{\link{writeResultHDF5}}.
#' 
#' This group can also contain a results overview dataset ('/Results/Overview'),
#' which summarises the results for all phenotypes. A results overview can be
#' read with \code{\link{readOverviewHDF5}}), or written with
#' \code{\link{writeOverviewHDF5}}.
#' }
#' 
#' If needed, datasets with a specified name can be read or written using the
#' functions \code{\link{readDatasetHDF5}} or \code{\link{writeDatasetHDF5}},
#' respectively. When reading or writing a matrix or array object, the object
#' is transposed (with R function \code{t}) or permuted (with R function
#' \code{aperm}), respectively. This is because HDF5 uses row-major order
#' (like the C language) while R uses column-major order (like Fortran).
#' For more information, see the
#' \href{https://www.hdfgroup.org/HDF5/doc/}{HDF5 documentation}. 
#' 
#' Note that attribute order is not preserved when reading/writing HDF5 objects,
#' and that datasets within a file cannot currently be overwritten.
#' 
#' @docType package
#' @keywords internal
#' @name HDF5 Utilities
NULL

# getGroupMemberNamesHDF5 ------------------------------------------------------
#' Get HDF5 group member names.
#' 
#' @param infile An input HDF5 file.
#' @param h5name HDF5 group name.
#' 
#' @return Vector of group member names. Returns \code{NULL}
#' if group is not present or has no members.
#' 
#' @importFrom methods new
#' @importFrom rhdf5 H5Lexists
#' @importFrom rhdf5 h5ls
#' @keywords internal
#' @rdname getGroupMemberNamesHDF5
getGroupMemberNamesHDF5 <- function(infile, h5name) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    h5name <- resolveH5ObjectName(h5name)
    
    member.names <- NULL
    
    h5stack <- methods::new('H5Stack', infile, h5name)
    
    on.exit( closeStack(h5stack) )
    
    if( ! rhdf5::H5Lexists(fileID(h5stack), h5name) ) {
        stop("HDF5 object ('", h5name, "') not found in file - '",  infile, "'")
    }

    group.info <- rhdf5::h5ls(peek(h5stack), recursive=FALSE)
    
    if ( nrow(group.info) > 0 ) {
        member.names <- group.info$name
    }
    
    return(member.names)
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
    return( getGroupMemberNamesHDF5(infile, 'Maps') )
}

# getObjectClassHDF5 -----------------------------------------------------------
#' Get R class of a HDF5 object.
#' 
#' If possible, the R class of the specified object is obtained from the
#' attributes 'R.class' and 'class', in that order. If neither attribute
#' is available, the object is loaded from the HDF5 file and its R class
#' determined from the resulting R object.
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
    
    result.names <- getGroupMemberNamesHDF5(infile, 'Results')
    
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
    phenotype.results <- getGroupMemberNamesHDF5(infile, h5name)
    return(phenotype.results)
}

# hasObjectHDF5 ----------------------------------------------------------------
#' Test if HDF5 file contains the named object.
#' 
#' @param infile An input HDF5 file.
#' @param h5name HDF5 object name.
#' 
#' @return TRUE if the given HDF5 file contains the named object;
#' FALSE otherwise.
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
#' Write attributes to a HDF5 object.
#'   
#' @param x List of attributes, or an R object whose attributes will be written.
#' @param h5obj HDF5 object to which attributes will be assigned.
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
h5writeAttributes <- function(x, h5obj) {
    
    attr.list <- attributes(x)

    if ( ! is.null(attr.list) ) {
        
        # If R object has 'dimnames' attribute, 
        # split into one attribute per dimension,
        # store dimension names separately.
        if ( 'dimnames' %in% names(attr.list) ) {
            
            stopifnot( 'dim' %in% names(attr.list) )
            
            exp.keys <- paste0( 'dimnames.', c('names', getIndices(attr.list[['dim']]) ) )
         
            if ( any( exp.keys %in% names(attr.list) ) ) {
                stop("cannot split 'dimnames' without overwriting existing attributes")
            }
            
            object.dimnames <- attr.list[['dimnames']]
            
            if ( ! is.null(object.dimnames) ) {
                
                if ( ! is.null( names(object.dimnames) ) ) {
                    attr.list[['dimnames.names']] <- names(object.dimnames)
                }
                
                indices <- which( ! sapply(object.dimnames, is.null) )
                
                for ( i in indices ) {
                    key <- paste0('dimnames.', i)
                    attr.list[[key]] <- object.dimnames[[i]]
                }
                
                attr.list[['dimnames']] <- NULL
            }
        }
        
        # If R object has row names, remove these.
        if ( 'row.names' %in% names(attr.list) ) {
            attr.list[['row.names']] <- NULL
        }
        
        # Write attributes.
        for ( i in getIndices(attr.list) ) {
            rhdf5::h5writeAttribute(h5obj=h5obj, name=names(attr.list)[i], 
                attr=attr.list[[i]])
        }
    }
    
    return( invisible() )
}

# joinH5ObjectNameParts --------------------------------------------------------
#' Join components of HDF5 object name.
#'  
#' @param components HDF5 name components.
#' @param relative Option indicating that the returned H5Object name should be
#' relative (i.e. should not contain a leading forward slash \code{'/'}).
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
        
        h5name <- paste(components, collapse='/')
        
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
#' element names are of the form 'ELT1', where 1 can be any positive integer. 
#' Leading zeros are used to pad element numbers, depending on the group size, 
#' so that the generated name strings will sort in numerical order.
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

# makeDefaultResultName --------------------------------------------------------
#' Get default QTL analysis result name.
#'   
#' @param result Object containing a QTL analysis result.
#'      
#' @return Default QTL analysis result name.
#'  
#' @keywords internal
#' @rdname makeDefaultResultName
makeDefaultResultName <- function(result) {
    
    classes <- const$default.result.names$class
    names <- const$default.result.names$name
    
    matching <- classes %in% class(result)
    if ( length(classes[matching]) == 1 ) {
        default.name <- names[matching]
    } else if ( length(classes[matching]) > 1 ) {
        stop("cannot get unique default name - result matches multiple classes - '", 
            toString(classes[matching]), "'")
    } else {
        stop("cannot get default name - unknown result class")
    }
    
    return(default.name)
}

# makeGroupObjectNames ---------------------------------------------------------
#' Make valid HDF5 group object names.
#' 
#' Group names are made based on the proposed group names or the group size. 
#' Default group element names will replace every group name that is either 
#' an empty string or NA value. Default group element names are of the form 
#' 'ELT1', where the number indicates the position of the element within the 
#' given group. 
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
        
        group.indices <- getIndices(group.names)
        
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
#' The named HDF5 dataset is read from the specified HDF5 file using the method
#' for the dataset type, as determined from the \code{'R.class'} or \code{'class'}
#' attribute of the HDF5 object.
#' 
#' An \pkg{R/qtl} \code{map} object is not suited for HDF5 format, but can be 
#' stored as a \code{mapframe} object, which can be read from file and then 
#' converted back to a \code{map}.
#'    
#' @param infile An input HDF5 file.
#' @param h5name HDF5 dataset name.
#' @param ... Further arguments (see below).
#' @param col.name For a \code{data.frame}, if any column has a name matching
#' this parameter, that column is extracted from the \code{data.frame} and the
#' rownames are set from its contents.
#' 
#' @return R object corresponding to the named HDF5 object.
#' 
#' @export
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
readDatasetHDF5.data.frame <- function(infile, h5name, col.name='rownames') {
    
    dataset <- readDatasetHDF5.default(infile, h5name)
    
    stopifnot( 'R.colClasses' %in% names( attributes(dataset) ) )
    
    colClasses <- attr(dataset, 'R.colClasses')
    dataset <- coerceDataFrame(dataset, colClasses)
    attr(dataset, 'R.colClasses') <- NULL
    
    if ( col.name %in% colnames(dataset) && ! hasRownames(dataset) ) {
        dataset <- setRownamesFromColumn(dataset, col.name=col.name)
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
        
        for ( i in getIndices(dims) ) {
            
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
    
    child.names <- getGroupMemberNamesHDF5(infile, h5name)
    
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
    
    dataset <- readDatasetHDF5.data.frame(infile, h5name, col.name='id')
    
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
    
    dataset <- readDatasetHDF5.data.frame(infile, h5name, col.name='id')
    
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
    
    for ( i in getIndices(dataset) ) {
        
        colnames(dataset[[i]]) <- attr(dataset, 'phenotypes')
    }
    attr(dataset, 'phenotypes') <- NULL
    
    return(dataset)
}

# readMapHDF5 ------------------------------------------------------------------
#' Read \code{map} from a HDF5 file.
#'   
#' @details A \code{mapframe} dataset is read from the 'Maps' group at the root
#' of the HDF5 file, and converted to an \pkg{R/qtl} \code{map} object.
#' 
#' @param infile An input HDF5 file.
#' @param name Map name.
#'      
#' @return An \pkg{R/qtl} \code{map} object. Returns \code{NULL}
#' if no map is found with the specified name.
#' 
#' @export
#' @keywords internal
#' @rdname readMapHDF5
readMapHDF5 <- function(infile, name) {
    
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
#' @return A \code{data.frame} containing QTL analysis results overview.
#' 
#' @export
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
#' @keywords internal
#' @rdname readResultHDF5
readResultHDF5 <- function(infile, phenotype, name) {
    
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
        
        res <- sub('/$', '', res)
        
    } else {
        
        res <- '/'
    }
    
    return(res)
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
#' The given R object is written to the named HDF5 dataset using the method for 
#' the dataset type, as determined by the class of the R object. Note that 
#' although \pkg{R/qtl} \code{map} objects are not ideally suited for HDF5 
#' format, they can be written to a HDF5 file as a \code{mapframe} object.
#'  
#' @param dataset R object.
#' @param outfile An output HDF5 file.
#' @param h5name HDF5 dataset name.
#' @param ... Further arguments (see below).
#' @param col.name For a \code{data.frame}, rownames are removed, and then
#' inserted into a column of the main table, with its column name taken
#' from this parameter.
#' 
#' @export
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
#' @keywords internal
#' @rdname writeDatasetHDF5
writeDatasetHDF5 <- function(dataset, outfile, h5name, ...) {
    UseMethod('writeDatasetHDF5', dataset)
}

# writeDatasetHDF5.array -------------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.array <- function(dataset, outfile, h5name, ...) {
    dataset <- aperm(dataset) # NB: R is column-major, HDF5 is row-major
    writeDatasetHDF5.default(dataset, outfile, h5name)
}

# writeDatasetHDF5.data.frame --------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.data.frame <- function(dataset, outfile, h5name,
    col.name='rownames', ...) {

    stopif( 'R.colClasses' %in% names( attributes(dataset) ) )
    
    stopifnot( nrow(dataset) > 0 )
    stopifnot( ncol(dataset) > 0 )
    
    if ( hasRownames(dataset) ) {
        dataset <- setColumnFromRownames(dataset, col.name=col.name)
    }
    
    colClasses <- sapply(dataset, class)
    cols <- which( colClasses %in% c('factor', 'logical') )
    
    for ( col in cols ) {
        dataset[[col]] <- as.character(dataset[[col]])
    }
    
    attr(dataset, 'R.colClasses') <- colClasses
    
    writeDatasetHDF5.default(dataset, outfile, h5name)
    
    return( invisible() )
}

# writeDatasetHDF5.default -----------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.default <- function(dataset, outfile, h5name, ...) {

    stopifnot( isSingleString(outfile) )
    h5name <- resolveH5ObjectName(h5name)
    
    name.parts <- splitH5ObjectName(h5name)

    h5loc.name <- joinH5ObjectNameParts(name.parts[ -length(name.parts) ])
    dataset.name <- name.parts[ length(name.parts) ]
    
    h5stack <- methods::new('H5Stack', outfile, h5loc.name)
    
    on.exit( closeStack(h5stack) )
    
    rhdf5::h5writeDataset(h5loc=peek(h5stack), name=dataset.name, obj=dataset)

    h5stack <- push(h5stack, rhdf5::H5Dopen(peek(h5stack), dataset.name))
    h5writeAttributes(dataset, peek(h5stack))

    return( invisible() )
}

# writeDatasetHDF5.list --------------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.list <- function(dataset, outfile, h5name, ...) {
    
    writeObjectAttributesHDF5(dataset, outfile, h5name)
    
    child.names <- makeGroupObjectNames( group.names=names(dataset), 
        group.size=length(dataset) )
    
    for ( i in getIndices(dataset) ) {
        
        child.h5name <- joinH5ObjectNameParts( c(h5name, child.names[i]) )
        
        writeDatasetHDF5(dataset[[i]], outfile, child.h5name)
    }
    
    return( invisible() )
}

# writeDatasetHDF5.map ---------------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.map <- function(dataset, outfile, h5name, ...) {
    
    stopif( 'R.class' %in% names( attributes(dataset) ) )
    
    attr(dataset, 'R.class') <- class(dataset)
    
    dataset <- as.data.frame(dataset)
    
    writeDatasetHDF5.data.frame(dataset, outfile, h5name, col.name='id')
    
    return( invisible() )
}

# writeDatasetHDF5.mapframe ----------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.mapframe <- function(dataset, outfile, h5name, ...) {
    
    stopif( 'R.class' %in% names( attributes(dataset) ) )
    
    attr(dataset, 'R.class') <- class(dataset)
    
    dataset <- as.data.frame(dataset)
    
    writeDatasetHDF5.data.frame(dataset, outfile, h5name, col.name='id')

    return( invisible() )
}

# writeDatasetHDF5.matrix ------------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.matrix <- function(dataset, outfile, h5name, ...) {
    dataset <- t(dataset) # NB: R is column-major, HDF5 is row-major
    writeDatasetHDF5.default(dataset, outfile, h5name)
}

# writeDatasetHDF5.scanone -----------------------------------------------------
#' @export
#' @rdname writeDatasetHDF5
writeDatasetHDF5.scanone <- function(dataset, outfile, h5name, ...) {
    
    num.phenotypes <- length( getDatColIndices(dataset) )
    stopifnot( num.phenotypes == 1 )
    
    writeDatasetHDF5.mapframe(dataset, outfile, h5name)
    
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
    
    dataset <- insertColumn(dataset, col.index=1, col.name='perm', data=perm.indices)
    
    otherattributes(dataset) <- others
    
    attr(dataset, 'R.class') <- c('scanonebins', 'array')
    
    writeDatasetHDF5.data.frame(dataset, outfile, h5name)
    
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
    
    writeDatasetHDF5.data.frame(scanoneperm.dataset, outfile, h5name)
    
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
    
    writeDatasetHDF5.list(dataset, outfile, h5name)
    
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
    
    writeDatasetHDF5.data.frame(scantwoperm.dataset, outfile, h5name)
    
    return( invisible() )
}

# writeMapHDF5 -----------------------------------------------------------------
#' Write \code{map} to a HDF5 file.
#'   
#' @details An \pkg{R/qtl} \code{map} object is written to the 'Maps' group at
#' the root of the given HDF5 file. If no map name is specified, a default will
#' be assigned based on whether it is a genetic or physical map.
#' 
#' @param map An \pkg{R/qtl} \code{map} object.
#' @param outfile An output HDF5 file.
#' @param name Map name.
#' 
#' @export
#' @keywords internal
#' @rdname writeMapHDF5
writeMapHDF5 <- function(map, outfile, name=NULL) {
    
    stopifnot( 'map' %in% class(map) )
    
    if ( is.null(name) ) {
        name <- makeDefaultMapName(map)
    }
    
    h5name <- joinH5ObjectNameParts( c('Maps', name) )
    writeDatasetHDF5(map, outfile, h5name)
    
    return( invisible() )
}

# writeObjectAttributesHDF5 ----------------------------------------------------
#' Write attributes to a HDF5 object.
#'   
#' @param x List of attributes, or an R object whose attributes will be written.
#' @param outfile An output HDF5 file.
#' @param h5name HDF5 object name. 
#'  
#' @importFrom methods new
#' @keywords internal
#' @rdname writeObjectAttributesHDF5
writeObjectAttributesHDF5 <- function(x, outfile, h5name) {
    
    stopifnot( isSingleString(outfile) )
    
    h5stack <- methods::new('H5Stack', outfile, h5name)
    
    on.exit( closeStack(h5stack) )
    
    h5writeAttributes(x, peek(h5stack))
    
    return( invisible() )
}

# writeOverviewHDF5 ------------------------------------------------------------
#' Write QTL analysis results overview to a HDF5 file.
#'    
#' @param overview Object containing a QTL analysis results overview.
#' @param outfile An output HDF5 file.
#' 
#' @export
#' @keywords internal
#' @rdname writeOverviewHDF5
writeOverviewHDF5 <- function(overview, outfile) {
    
    # TODO: review result overview format.
    
    stopifnot( is.data.frame(overview) )
    stopifnot( all( colnames(overview) == c('Phenotype', 'Status', 'Comments') ) )
    stopifnot( is.character(overview[['Phenotype']]) )
    stopifnot( is.logical(overview[['Status']]) )
    stopifnot( is.character(overview[['Comments']]) )
    
    h5name <- joinH5ObjectNameParts( c('Results', 'Overview') )
    writeDatasetHDF5(overview, outfile, h5name)
    return( invisible() )
}

# writeResultHDF5 --------------------------------------------------------------
#' Write QTL analysis result to a HDF5 file.
#' 
#' @param result R object containing a QTL analysis result.
#' @param outfile An output HDF5 file.
#' @param phenotype Name of phenotype (or equivalent analysis unit).
#' @param name Name of QTL analysis result.
#' 
#' @export
#' @keywords internal
#' @rdname writeResultHDF5
writeResultHDF5 <- function(result, outfile, phenotype, name=NULL) {
    
    if ( is.null(name) ) {
        name <- makeDefaultResultName(result)
    }
    
    h5name <- joinH5ObjectNameParts( c('Results', phenotype, name) )
    writeDatasetHDF5(result, outfile, h5name)
    return( invisible() )
}

# End of hdf5.R ################################################################
