# Start of run_pushmap.R #######################################################

# run_pushmap ------------------------------------------------------------------
#' Push map into a data file.
#' 
#' Add or replace the map of the specified CSV or HDF5 file with the
#' data from the given map file. If the \code{mapname} parameter is
#' given for a HDF5 data file, the map is assigned that name and when
#' being written to the HDF5 file. Otherwise, a default map name
#' (e.g. \code{'Genetic Map'}) is used.
#' 
#' @param mapfile input map CSV file [required]
#' @param datafile file with map data [required]
#' @param mapname name of map (if applicable)
#' @param require.mapunit require map units in input
#' @param include.mapunit include map units in output
#' 
#' @concept shmootl:utilities
#' @export
#' @family pipeline functions
#' @rdname run_pushmap
run_pushmap <- function(mapfile=NA_character_, datafile=NA_character_,
    mapname=NA_character_, require.mapunit=TRUE, include.mapunit=TRUE) {
    
    stopifnot( isSingleString(mapfile) )
    stopifnot( isSingleString(datafile) )
    
    mapname <- if ( ! identical(mapname, NA_character_) ) { mapname } else { NULL }
    
    datafile.format <- inferFormatFromFilename(datafile)
    
    # Create output temp file.
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    # If datafile exists, copy to temp file.
    if ( file.exists(datafile) ) {
        file.copy(datafile, tmp)
    }
    
    # Read map from mapfile.
    cross.map <- readMapCSV(mapfile, require.mapunit=require.mapunit)
    
    # Push map into temp file.
    if ( datafile.format == 'CSV' ) {
        
        if ( ! is.null(mapname) ) {
            stop("cannot specify map name for a CSV file")
        }
        
        writeMapCSV(cross.map, tmp, include.mapunit=include.mapunit)
        
    } else if ( datafile.format == 'HDF5' ) {
        
        writeMapHDF5(cross.map, tmp, name=mapname, overwrite=TRUE)
        
    } else {
        
        stop("cannot push map into ", datafile.format, " format file")
    }
    
    # Move temp file to final output file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, datafile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_pushmap.R #########################################################