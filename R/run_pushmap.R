# Start of run_pushmap.R #######################################################

# run_pushmap ------------------------------------------------------------------
#' Push map into cross/geno data file.
#' 
#' @description This script replaces the map of the specified data file with the
#' data from the given map file.
#' 
#' @param mapfile input map CSV file
#' @param datafile cross/geno CSV file
#' @param require.mapunit require map units in input
#' @param include.mapunit include map units in output
#' 
#' @export
#' @rdname run_pushmap
run_pushmap <- function(mapfile, datafile, require.mapunit=TRUE,
    include.mapunit=TRUE) {
    
    # Create output temp file.
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    # Copy datafile to temp file.
    file.copy(datafile, tmp)
    
    guess <- sniffCSV(datafile)
    
    if ( ! guess %in% c('cross', 'geno') ) {
        stop("cannot push map into file with ", guess, " data - '", datafile,"'")
    }
    
    # Read map from mapfile.
    cross.map <- readMapCSV(mapfile, require.mapunit=require.mapunit)
    
    # Push map into temp file.
    writeMapCSV(cross.map, tmp, include.mapunit=include.mapunit)
    
    # Move temp file to final output file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, datafile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_pushmap.R #########################################################