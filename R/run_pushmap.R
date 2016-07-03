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
    
    guess <- sniffCSV(datafile)
    
    if ( is.null(guess) ) {
        stop("cannot push map - unknown input data")
    }
    
    map <- readMapCSV(mapfile)
    
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    if ( guess == 'cross' ) {
        
        cross <- readCrossCSV(datafile, require.mapunit=require.mapunit)
        
        cross <- qtl::replace.map(cross, map)
        
        writeCrossCSV(cross, tmp, include.mapunit=include.mapunit)
        
    } else if ( guess == 'geno' ) {
        
        geno <- readGenoCSV(datafile, require.mapunit=require.mapunit)
        
        geno <- pushMap(geno, map)
        
        writeGenoCSV(geno, tmp, include.mapunit=include.mapunit)
        
    } else {
        
        stop("cannot push map into ", guess," data file")
    }
    
    # Move temp file to final output file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, datafile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_pushmap.R #########################################################