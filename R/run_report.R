# Start of run_report.R ########################################################

# run_report -------------------------------------------------------------------
#' Create report of QTL analysis results.
#' 
#' Given a scan result HDF5 file, create a report of the results of the QTL
#' analysis. Results can be output for specific phenotypes and analyses.
#' These output constraints are applied with parameters \code{pheno} and
#' \code{scans}, respectively.
#' 
#' @param h5file HDF5 scan file [required]
#' @param report scan report file [required]
#' @param pheno phenotypes to output [default: all]
#' @param scans analyses to output [default: all]
#' 
#' @concept shmootl:processing
#' @export
#' @family pipeline functions
#' @rdname run_report
run_report <- function(h5file=NA_character_, report=NA_character_,
    pheno=character(), scans=character()) {
    
    stopifnot( isSingleString(h5file) )
    stopifnot( isSingleString(report) )
    
    pheno <- if ( ! identical(pheno, character()) ) { pheno } else { NULL }
    scans <- if ( ! identical(scans, character()) ) { scans } else { NULL }
    
    # Get report file extension.
    report.ext <- tools::file_ext(report)
    
    # Infer report file format.
    report.format <- inferFormatFromFilename(report)
    
    # Create temp report file, ensure will be removed.
    tmp <- tempfile()
    on.exit( file.remove(tmp), add=TRUE )
    
    # Write report to temp file.
    if ( report.format == 'PDF' ) {
        writeReportPDF(h5file, tmp, phenotypes=pheno, analyses=scans)
    } else {
        stop("report output not supported for ", report.format, " format")
    }
    
    # If report written without error, move temp file to final report file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, report, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_report.R ##########################################################