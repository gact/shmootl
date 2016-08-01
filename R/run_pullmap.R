# Start of run_pullmap.R #######################################################

# run_pullmap ------------------------------------------------------------------
#' Pull map from data file.
#' 
#' @description This script pulls the map from a CSV or HDF5 file, and writes it
#' to a separate map CSV file. If a HDF5 data file contains multiple maps, the
#' \code{mapname} parameter must be used to specify which map to pull.
#' 
#' @param datafile cross/geno CSV file
#' @param mapfile output map CSV file
#' @param mapname name of map (if applicable)
#' @param require.mapunit require map units in input
#' @param include.mapunit include map units in output
#' 
#' @export
#' @rdname run_pullmap
run_pullmap <- function(datafile, mapfile, mapname=NA_character_,
    require.mapunit=TRUE, include.mapunit=TRUE) {
    
    stopifnot( isSingleString(mapfile) )
    
    mapname <- if ( ! identical(mapname, NA_character_) ) { mapname } else { NULL }
    
    datafile.ext <- tools::file_ext(datafile)
    
    if ( datafile.ext %in% const$ext$csv ) {
        
        if ( ! is.null(mapname) ) {
            stop("cannot specify map name for a CSV file")
        }
        
        cross.map <- readMapCSV(datafile, require.mapunit=require.mapunit)
        
    } else if ( datafile.ext %in% const$ext$hdf5 ) {
        
        mapname <- resolveMapNameHDF5(datafile, mapname)
        
        cross.map <- readMapHDF5(datafile, mapname)
        
    } else {
        
        stop("cannot pull map - unknown extension on file '", datafile, "'")
    }
    
    # Create output temp file.
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    # Write map to temp file.
    writeMapCSV(cross.map, tmp, include.mapunit=include.mapunit)
    
    # Move temp file to final map file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, mapfile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_pullmap.R #########################################################