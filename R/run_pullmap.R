# Start of run_pullmap.R #######################################################

# run_pullmap ------------------------------------------------------------------
#' Pull map from a data file.
#' 
#' Pull the map from a CSV or HDF5 file, and write it to a separate map CSV
#' file. If a HDF5 data file contains multiple maps, the \code{mapname}
#' parameter must specify which map to pull.
#' 
#' @param datafile file containing map data [required]
#' @param mapfile output map CSV file [required]
#' @param mapname name of map (if applicable)
#' @param require.mapunit require map units in input
#' @param include.mapunit include map units in output
#' 
#' @concept shmootl:utilities
#' @export
#' @family pipeline functions
#' @rdname run_pullmap
run_pullmap <- function(datafile=NA_character_, mapfile=NA_character_,
    mapname=NA_character_, require.mapunit=TRUE, include.mapunit=TRUE) {
    
    stopifnot( isSingleString(datafile) )
    stopifnot( isSingleString(mapfile) )
    
    mapname <- if ( ! identical(mapname, NA_character_) ) { mapname } else { NULL }
    
    datafile.format <- inferFormatFromFilename(datafile)
    
    if ( datafile.format == 'CSV' ) {
        
        if ( ! is.null(mapname) ) {
            stop("cannot specify map name for a CSV file")
        }
        
        cross.map <- readMapCSV(datafile, require.mapunit=require.mapunit)
        
    } else if ( datafile.format == 'HDF5' ) {
        
        mapname <- resolveMapNameHDF5(datafile, mapname)
        
        cross.map <- readMapHDF5(datafile, mapname)
        
    } else {
        
        stop("cannot pull map from ", datafile.format, " format file")
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