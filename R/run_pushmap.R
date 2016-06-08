# Start of run_pushmap.R #######################################################

# run_pushmap ------------------------------------------------------------------
#' Push map into cross/geno data file.
#' 
#' @description This script reads cross or genotype data from the specified
#' input file, replaces its map with that of the given map file, and writes
#' the result to the specified output file.
#' 
#' @param infile input cross/geno CSV file
#' @param mapfile input map CSV file
#' @param outfile output cross/geno CSV file
#' @param require.mapunit require map units in input
#' @param include.mapunit include map units in output
#' 
#' @export
#' @rdname run_pushmap
run_pushmap <- function(infile, mapfile, outfile, require.mapunit=TRUE,
    include.mapunit=TRUE) {
    
    stopifnot( isSingleString(outfile) )
    
    guess <- sniffCSV(infile)
    
    if ( is.null(guess) ) {
        stop("cannot push map - unknown input data")
    }
    
    map <- readMapCSV(mapfile)
    
    if ( guess == 'cross' ) {
        
        cross <- readCrossCSV(infile, require.mapunit=require.mapunit)
        
        cross <- qtl::replace.map(cross, map)
        
        writeCrossCSV(cross, outfile, include.mapunit=include.mapunit)
        
    } else if ( guess == 'geno' ) {
        
        geno <- readGenoCSV(infile, require.mapunit=require.mapunit)
        
        geno <- pushMap(geno, map)
        
        writeGenoCSV(geno, outfile, include.mapunit=include.mapunit)
        
    } else {
        
        stop("cannot push map into ", guess," data file")
    }
    
    return( invisible() )
}

# End of run_pushmap.R #########################################################