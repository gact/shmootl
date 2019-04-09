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
#' genetic map to a corresponding physical map. QTL intervals can be annotated
#' for specific phenotypes and analyses with the \code{pheno} and \code{scans}
#' parameters, respectively.
#' 
#' @param infile input HDF5 scan file [required]
#' @param annofile GFF annotation file [required]
#' @param outfile annotated HDF5 scan file [required]
#' @param pheno phenotypes [default: all]
#' @param scans analyses [default: all]
#' 
#' @concept shmootl:processing
#' @export
#' @family pipeline functions
#' @rdname run_annoqtl
run_annoqtl <- function(infile=NA_character_, annofile=NA_character_,
    outfile=NA_character_, pheno=character(), scans=character()) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    stopifnot( isSingleString(annofile) )
    stopifnot( isSingleString(outfile) )
    
    pheno <- if ( ! identical(pheno, character()) ) { pheno } else { NULL }
    scans <- if ( ! identical(scans, character()) ) { scans } else { NULL }
    
    # Get result info.
    rinfo <- getResultInfoHDF5(infile, phenotypes=pheno, analyses=scans)
    
    # Get analyses of interest.
    analyses <- character()
    for ( phenotype in names(rinfo[[infile]]) ) {
        for ( analysis in names(rinfo[[infile]][[phenotype]]) ) {
            if ( 'QTL Intervals' %in% rinfo[[infile]][[phenotype]][[analysis]] ) {
                analyses <- union(analyses, analysis)
            }
        }
    }
    
    if ( length(analyses) == 0 ) {
        warning("no relevant results found for annotation in file '", infile, "'")
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
    
    # Annotate QTL intervals for each analysis..
    for ( analysis in analyses ) {
        
        # ..and phenotype.
        for ( phenotype in names(rinfo[[infile]]) ) {
            
            # If QTL intervals present for this phenotype and
            # analysis, get any corresponding QTL features.
            if ( analysis %in% names(rinfo[[infile]][[phenotype]]) &&
                'QTL Intervals' %in% rinfo[[infile]][[phenotype]][[analysis]] ) {
                
                # Read QTL intervals from input HDF5 scan file.
                qtl.intervals <- readResultHDF5(infile, phenotype,
                    analysis, 'QTL Intervals')
                
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
                    if (length(qtl.features) > 0) {
                        h5name <- joinH5ObjectNameParts( c('Results',
                            phenotype, analysis, 'QTL Intervals') )
                        writeDatasetHDF5(qtl.intervals, tmp, h5name)
                        updated.objects <- c(updated.objects, h5name)
                    }
                }
                
                # Get annotated features in QTL intervals.
                qtl.features <- getQTLFeatures(qtl.intervals, features)
                
                # Remove empty QTL feature data-frames.
                qtl.features <- qtl.features[ sapply(qtl.features,
                    function(x) nrow(x) > 0) ]
                
                # Write Scanone QTL features to temp file.
                h5name <- joinH5ObjectNameParts( c('Results',
                    phenotype, analysis, 'QTL Features') )
                writeDatasetHDF5(qtl.features, tmp, h5name)
                updated.objects <- c(updated.objects, h5name)
            }
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
