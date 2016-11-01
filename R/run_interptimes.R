# Start of run_interptimes.R ###################################################

# run_interptimes --------------------------------------------------------------
#' Interpolate time-series phenotypes in \pkg{R/qtl} CSV file.
#' 
#' Read data in an \pkg{R/qtl} cross CSV file. If the input \code{cross} has
#' time-series phenotypes, interpolate phenotype values in any gaps in the
#' time series, then write the resulting \code{cross} to the output file.
#' 
#' @param datafile cross CSV file [required]
#' @param tol time-step equality tolerance
#' 
#' @template author-yue-hu
#' 
#' @concept shmootl:utilities
#' @export
#' @family pipeline functions
#' @rdname run_interptimes
run_interptimes <- function(datafile=NA_character_, tol=1e-05) {
    
    stopifnot( isSingleString(datafile) )
    
    # Read cross input file.
    cross <- readCrossCSV(datafile)
    
    # Check whether the input phenotypes are a time-series.
    stopifnot( hasTimeSeriesPhenotypes(cross, tol=tol) )
    
    # Do linear interpolation of gaps in time-series phenotypes.
    # TODO: add a 'method' parameter to enable different interpolation methods
    cross <- interpTimeSeries(cross, tol=tol)
    
    # Create temp output file, ensure will be removed.
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    # Write cross to temp output file.
    writeCrossCSV(cross, tmp)
    
    # Move temp file to final output file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, datafile, overwrite=TRUE)
    
    return( invisible() )
} 

# End of run_interptimes.R #####################################################