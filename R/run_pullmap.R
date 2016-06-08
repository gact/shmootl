# Start of run_pullmap.R #######################################################

# run_pullmap ------------------------------------------------------------------
#' Pull map from cross/geno data.
#' 
#' @description This script pulls the map from an R/qtl cross or genotype file,
#' and writes it to a separate map file.
#' 
#' @param infile input cross/geno CSV file
#' @param mapfile output map CSV file
#' @param require.mapunit require map units in input
#' @param include.mapunit include map units in output
#' 
#' @export
#' @rdname run_pullmap
run_pullmap <- function(infile, mapfile, require.mapunit=TRUE,
    include.mapunit=TRUE) {
    
    stopifnot( isSingleString(mapfile) )
    
    guess <- sniffCSV(infile)
    
    if ( is.null(guess) ) {
        stop("cannot pull map - unknown input data")
    }
    
    if ( ! hasMapCSV(infile) ) {
        stop("no map data found in file '", infile,"'")
    }
    
    if ( guess == 'cross' ) {
        
        cross <- readCrossCSV(infile, require.mapunit=require.mapunit)
        cross.map <- qtl::pull.map(cross)
        
    } else if ( guess == 'geno' ) {
        
        geno <- readGenoCSV(infile, require.mapunit=require.mapunit)
        cross.map <- pullMap(geno)
        
    } else {
        
        stop("cannot pull map from ", guess," data")
    }
    
    writeMapCSV(cross.map, mapfile, include.mapunit=include.mapunit)
    
    return( invisible() )
}

# End of run_pullmap.R #########################################################