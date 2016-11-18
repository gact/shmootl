# Start of run_report.R ########################################################

# run_report -------------------------------------------------------------------
#' Create report of QTL analysis results.
#' 
#' Given a scan result HDF5 file, create a report of the results of the QTL
#' analysis. Results can be output for specific phenotypes, analyses, and,
#' in the case of Excel output, worksheets. These output constraints are
#' applied with parameters \code{pheno}, \code{scans}, and \code{sheets},
#' respectively.
#' 
#' @param h5file HDF5 scan file [required]
#' @param report scan report file [required]
#' @param pheno phenotypes to output [default: all]
#' @param scans analyses to output [default: all]
#' @param sheets worksheets to output (Excel only)
#' @param scanfile.pattern scan file name pattern (Excel only)
#' 
#' @concept shmootl:processing
#' @export
#' @family pipeline functions
#' @rdname run_report
run_report <- function(h5file=NA_character_, report=NA_character_,
    pheno=character(), scans=character(), sheets=character(),
    scanfile.pattern=NA_character_) {
    
    stopifnot( isSingleString(h5file) )
    stopifnot( isSingleString(report) )
    
    pheno <- if ( ! identical(pheno, character()) ) { pheno } else { NULL }
    scans <- if ( ! identical(scans, character()) ) { scans } else { NULL }
    sheets <- if ( ! identical(sheets, character()) ) { sheets } else { NULL }
    scanfile.pattern <- if ( ! identical(scanfile.pattern, NA_character_) ) {
        scanfile.pattern } else { NULL }
    
    # Get report file extension.
    report.ext <- tools::file_ext(report)
    
    # Infer report file format.
    report.format <- inferFormatFromFilename(report)
    
    # Create temp report file, ensure will be removed.
    # NB: must have correct extension so that Excel format can be set automatically.
    tmp <- tempfile( fileext=paste0('.', report.ext) )
    on.exit( file.remove(tmp), add=TRUE )
    
    # Write report to temp file.
    if ( report.format == 'PDF' ) {
        
        xovars <- c('sheets', 'scanfile.pattern')
        for ( xovar in xovars ) {
            if ( ! is.null( get(xovar) ) ) {
                warning("parameter '", xovar, "' ignored for ",
                    report.format, " format report")
            }
        }
        
        writeReportPDF(h5file, tmp, phenotypes=pheno, analyses=scans)
        
    } else if ( report.format == 'Excel' ) {
        
        writeReportExcel(h5file, tmp, phenotypes=pheno, analyses=scans,
            worksheets=sheets, scanfile.pattern=scanfile.pattern)
        
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