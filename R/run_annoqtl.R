# Start of run_annoqtl.R #######################################################

# run_annoqtl ------------------------------------------------------------------
#' Annotate QTL intervals.
#' 
#' QTL intervals in a HDF5 scan file are annotated with features from the GFF
#' annotation file, and a copy of the HDF5 scan file is output, in which QTL
#' feature annotation has been added. The GFF annotation file must relate to
#' the same reference genome as was used in QTL analysis. In addition, physical
#' positions of QTL intervals must be present, either directly in the QTL
#' intervals themselves, or indirectly by pairing coordinates in the cross
#' genetic map to a corresponding physical map.
#' 
#' @param infile input HDF5 scan file [required]
#' @param annofile GFF annotation file [required]
#' @param outfile annotated HDF5 scan file [required]
#' 
#' @concept shmootl:processing
#' @export
#' @family pipeline functions
#' @rdname run_annoqtl
run_annoqtl <- function(infile=NA_character_, annofile=NA_character_,
    outfile=NA_character_) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    stopifnot( isSingleString(annofile) )
    stopifnot( isSingleString(outfile) )
    
    results.sought <- list(
        'Scanone' = c('QTL Intervals')
    )
    
    result.info <- list()
    roi <- list()
    
    if ( hasObjectHDF5(infile, 'Results') ) {
        
        # Get result info for HDF5 scan file.
        for ( phenotype in getResultPhenotypesHDF5(infile) ) {
            analyses <- getResultAnalysesHDF5(infile, phenotype)
            result.info[[phenotype]] <- lapply( analyses, function(analysis)
                getResultNamesHDF5(infile, phenotype, analysis) )
            names(result.info[[phenotype]]) <- analyses
        }
        
        # Get results of interest.
        for ( phenotype in names(result.info) ) {
            for ( analysis in names(result.info[[phenotype]]) ) {
                for ( result in result.info[[phenotype]][[analysis]] ) {
                    if ( analysis %in% names(results.sought) &&
                         result %in% results.sought[[analysis]] ) {
                        roi[[analysis]] <- union(roi[[analysis]], result)
                    }
                }
            }
        }
    }
    
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
    
    # Init HDF5 object names of updated objects.
    updated.objects <- character()
    
    # Annotated QTL intervals for each phenotype.
    for ( phenotype in names(result.info) ) {
        
        # If Scanone QTL intervals present for this phenotype,
        # get any corresponding QTL features.
        if ( 'Scanone' %in% names(result.info[[phenotype]]) &&
            'QTL Intervals' %in% result.info[[phenotype]][['Scanone']] ) {
            
            # Read QTL intervals from input HDF5 scan file.
            qtl.intervals <- readResultHDF5(infile, phenotype,
                'Scanone', 'QTL Intervals')
            
            # Ensure that Scanone QTL intervals have (estimated) physical map positions.
            if ( ! hasPhysicalPositions(qtl.intervals) ) {
                
                # Ensure mapkey option has been set 
                # from map data of HDF5 scan file.
                if ( is.null(map.key) ) {
                    
                    if ( ! hasCrossHDF5(infile) ) {
                        stop("cannot annotate scan results -",
                            " no cross data in file '", infile, "'")
                    }
                    
                    if ( ! hasMapHDF5(infile, 'Physical Map') ) {
                        stop("cannot annotate scan results -",
                            " no physical map data in file '", infile, "'")
                    }
                    
                    # Get genetic and physical map data.
                    gmap <- qtl::pull.map( readCrossHDF5(infile) )
                    pmap <- readMapHDF5(infile, 'Physical Map')
                    
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
                h5name <- joinH5ObjectNameParts(
                    c('Results', phenotype, 'Scanone', 'QTL Intervals') )
                writeDatasetHDF5(qtl.intervals, tmp, h5name)
                updated.objects <- c(updated.objects, h5name)
            }
            
            # Get annotated features in QTL intervals.
            qtl.features <- getQTLFeatures(qtl.intervals, features)
            
            # Remove empty QTL feature data-frames.
            qtl.features <- qtl.features[ sapply(qtl.features,
                function(x) nrow(x) > 0) ]
            
            # Write Scanone QTL features to temp file.
            h5name <- joinH5ObjectNameParts(
                c('Results', phenotype, 'Scanone', 'QTL Features') )
            writeDatasetHDF5(qtl.features, tmp, h5name)
            updated.objects <- c(updated.objects, h5name)
        }
    }
    
    # Transfer all existing but unchanged objects
    # from existing HDF5 scan file to temp file.
    updated <- unique( unlist( lapply( updated.objects,
        function(h5name) getObjectNamesHDF5(tmp, h5name) ) ) )
    existing <- getObjectNamesHDF5(infile)
    unchanged <- setdiff(existing, updated) # NB: guarantees unique
    copyObjectsHDF5(infile, tmp, h5names=unchanged)
    
    # Move temp file to annotated HDF5 scan file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, outfile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_annoqtl.R #########################################################