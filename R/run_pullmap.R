# Start of run_pullmap.R #######################################################

# run_pullmap ------------------------------------------------------------------
#' Pull map from cross/geno data.
#' 
#' @description This script pulls the map from an R/qtl cross or genotype file,
#' and writes it to a separate map file.
#' 
#' @param datafile cross/geno CSV file
#' @param mapfile output map CSV file
#' @param require.mapunit require map units in input
#' @param include.mapunit include map units in output
#' 
#' @export
#' @rdname run_pullmap
run_pullmap <- function(datafile, mapfile, require.mapunit=TRUE,
    include.mapunit=TRUE) {
    
    stopifnot( isSingleString(mapfile) )
    
    cross.map <- readMapCSV(datafile, require.mapunit=require.mapunit)
    
    # Create output temp file.
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    writeMapCSV(cross.map, tmp, include.mapunit=include.mapunit)
    
    # Move temp file to final map file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, mapfile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_pullmap.R #########################################################