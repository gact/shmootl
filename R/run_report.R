# Start of run_report.R ########################################################

# run_report -------------------------------------------------------------------
#' Create report of QTL analysis results.
#' 
#' Given a scan result HDF5 file, create a
#' report of the results of the QTL analysis.
#' 
#' @param h5file HDF5 scan file [required]
#' @param report scan report file [required]
#' 
#' @concept shmootl:processing
#' @export
#' @family pipeline functions
#' @rdname run_report
run_report <- function(h5file=NA_character_, report=NA_character_) {
    
    stopifnot( isSingleString(h5file) )
    stopifnot( isSingleString(report) )
    
    # Get digest file extension.
    report.ext <- tools::file_ext(report)
    
    # Create temp report file, ensure will be removed.
    tmp <- tempfile()
    on.exit( file.remove(tmp), add=TRUE )
    
    # Write report to temp file.
    success <- tryCatch({
        
        if ( report.ext %in% const$ext$pdf ) {
            writeReportPDF(h5file, tmp)
        } else {
            stop("cannot write report - unknown extension on file '", report, "'")
        }
        
        result <- TRUE
        
    }, error=function(e) {
        result <- FALSE
    })
    
    # If report written without error, move temp file to final report file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    if (success) {
        file.copy(tmp, report, overwrite=TRUE)
    }
    
    return( invisible() )
}

# End of run_report.R ##########################################################