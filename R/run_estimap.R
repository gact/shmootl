# Start of run_estimap.R #######################################################

# run_estimap ------------------------------------------------------------------
#' Estimate map from cross/geno data.
#' 
#' @description This script reads an R/qtl cross or genotype file, and outputs a
#' map estimated from the input genotype data. Note that a genetic map cannot be
#' estimated from enumerated genotypes, as these are generated independently for
#' each marker.
#' 
#' @param datafile input cross/geno CSV file
#' @param mapfile output map CSV file
#' @param n.cluster number of threads
#' @param error.prob genotyping error rate
#' @param map.function genetic map function
#' 
#' @export
#' @rdname run_estimap
run_estimap <- function(datafile, mapfile, n.cluster=1L, error.prob=0.0001,
    map.function=c("haldane","kosambi","c-f","morgan")) {
    
    stopifnot( isSingleString(mapfile) )
    stopifnot( isSinglePositiveNumber(n.cluster) )
    
    guess <- sniffCSV(datafile)
    
    if ( guess == 'cross' ) {
        
        cross <- readCrossCSV(datafile, require.mapunit=FALSE)
        
    } else if ( guess == 'geno' ) {
        
        geno <- readGenoCSV(datafile, require.mapunit=FALSE)
        
        cross.info <- attr(geno, 'info')
        samples <- getSamples(cross.info)
        pheno <- makePlaceholderPheno(samples=samples)
        
        cross <- makeCross(geno, pheno)
        
    } else {
        
        stop("cannot estimate map from ", guess, " input data")
    }
    
    # Check for enumerated genotypes.
    if ( any( isEnumAllele( attr(cross, 'alleles') ) ) ) {
        stop("cannot estimate map from enumerated genotype data")
    }
    
    cross.map <- qtl::est.map(cross, error.prob=error.prob,
        map.function=map.function, offset=0, n.cluster=n.cluster)
    
    cross.map <- as.mapframe(cross.map)
    
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    writeMapframeCSV(cross.map, tmp)
    
    # Move temp file to final map file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, mapfile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_estimap.R #########################################################