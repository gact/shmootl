# Start of run_annoqtl.R #######################################################

# run_annoqtl ------------------------------------------------------------------
#' Annotate QTL intervals.
#' 
#' QTL intervals in a scan result HDF5 file are annotated with features from the
#' GFF annotation file, and a copy of the scan result file is output, in which
#' QTL feature annotation has been added. The GFF annotation file must relate to
#' the same reference genome as was used in QTL analysis. In addition, physical
#' positions of QTL intervals must be present, either directly in the QTL
#' intervals themselves, or indirectly by pairing coordinates in genetic and
#' physical maps for the cross data.
#' 
#' @param infile scan result file
#' @param annofile GFF annotation file
#' @param outfile annotated scan result file
#' 
#' @concept shmootl:pipelines
#' @export
#' @family pipeline functions
#' @rdname run_annoqtl
run_annoqtl <- function(infile=NA_character_, annofile=NA_character_,
    outfile=NA_character_) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    stopifnot( isSingleString(outfile) )
    
    results.sought <- 'Scanone QTL Intervals'
    results.found <- list()
    result.info <- list()
    
    if ( hasObjectHDF5(infile, 'Results') ) {
        
        # Get phenotypes from input scan result file.
        phenotypes <- getPhenotypesHDF5(infile)
        
        # Get result info for phenotypes.
        result.info <- lapply(phenotypes, function(phenotype)
            getResultNamesHDF5(infile, phenotype))
        names(result.info) <- phenotypes
        
        # Get results of interest for each phenotype.
        results.found <- lapply(result.info, function(results)
            results[ results %in% results.sought ])
    }
    
    # Get set of results of interest.
    roi <- unique( unlist(results.found) )
    
    if ( length(roi) == 0 ) {
        return( invisible() )
    }
    
    # Get genome features.
    features <- readFeaturesGFF(annofile)
    
    # Create temp output file, ensure will be removed.
    tmp <- tempfile()
    on.exit( file.remove(tmp), add=TRUE )
    
    # Init mapkey for converting genetic to physical map positions.
    map.key <- NULL
    
    # Init HDF5 object names of updated datasets.
    start.points <- character()
    
    # Annotated QTL intervals for each phenotype.
    for ( phenotype in phenotypes ) {
        
        # If Scanone QTL intervals present for this phenotype,
        # get any corresponding QTL features.
        if ( 'Scanone QTL Intervals' %in% result.info[[phenotype]] ) {
            
            # Read QTL intervals from input scan result file.
            qtl.intervals <- readResultHDF5(infile, phenotype, 'Scanone QTL Intervals')
            
            # Ensure that Scanone QTL intervals have (estimated) physical map positions.
            if ( ! hasPhysicalPositions(qtl.intervals) ) {
                
                # Ensure mapkey option has been set 
                # from map data of scan result file.
                if ( is.null(map.key) ) {
                    
                    # Get genetic and physical map data.
                    gmap <- readMapHDF5(infile, 'Genetic Map')
                    pmap <- readMapHDF5(infile, 'Physical Map')
                    
                    if ( is.null(gmap) || is.null(pmap) ) {
                        stop("cannot annotate scan results - insufficient map data in file '", infile, "'")
                    }
                    
                    # Create mapkey for converting between
                    # genetic and physical map positions.
                    map.key <- mapkey(gmap, pmap)
                    
                    # Set mapkey option.
                    prev.mapkey <- mapkeyOpt(map.key)
                    on.exit( mapkeyOpt(prev.mapkey), add=TRUE )
                }
                
                # Estimate physical map positions of QTL intervals.
                qtl.intervals <- estPhysicalPositions(qtl.intervals)
                
                # Write updated QTL intervals to temp file.
                h5name <- joinH5ObjectNameParts( c('Results', phenotype, 'Scanone QTL Intervals') )
                writeDatasetHDF5(qtl.intervals, tmp, h5name)
                start.points <- c(start.points, h5name)
            }
            
            # Get annotated features in QTL intervals.
            qtl.features <- getQTLFeatures(qtl.intervals, features)
            
            # Remove empty QTL feature data-frames.
            qtl.features <- qtl.features[ sapply(qtl.features, function(x) nrow(x) > 0) ]
            
            # Write QTL features to temp file.
            h5name <- joinH5ObjectNameParts( c('Results', phenotype, 'Scanone QTL Features') )
            writeDatasetHDF5(qtl.features, tmp, h5name)
            start.points <- c(start.points, h5name)
        }
    }
    
    # Transfer all existing but unchanged objects
    # from existing scan result file to temp file.
    updated <- unique( unlist( lapply( start.points,
        function(h5name) getObjectNamesHDF5(tmp, h5name) ) ) )
    existing <- getObjectNamesHDF5(infile)
    unchanged <- setdiff(existing, updated) # NB: guarantees unique
    copyObjectsHDF5(infile, tmp, h5names=unchanged)
    
    # Move temp file to annotated scan result file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, outfile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_annoqtl.R #########################################################