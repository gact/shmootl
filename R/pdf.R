# Start of pdf.R ###############################################################

# writeReportPDF ---------------------------------------------------------------
#' Write a PDF report of QTL scan results.
#'  
#' @param scanfile A scan result HDF5 file.
#' @param report Path of output report PDF file.
#' @param analyses Analyses for which results should be included in the report
#' file. If none are specified, results are output for all available analyses.
#' 
#' @export
#' @importFrom grDevices cairo_pdf
#' @importFrom grDevices dev.cur
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @rdname writeReportPDF
writeReportPDF <- function(scanfile, report, analyses=NULL) {
    
    stopifnot( isSingleString(scanfile) )
    stopifnot( file.exists(scanfile) )
    stopifnot( isSingleString(report) )
    
    # Set possible results to be sought in scan file.
    results.sought <- list(
        'Scanone' = c('Result')
    )
    
    # Resolve relevant analyses.
    analyses.sought <- unique( resolveAnalysisTitle(analyses) )
    
    # Set actual results sought for relevant analyses.
    results.sought <- results.sought[ names(results.sought) %in% analyses.sought ]
    
    result.info <- list()
    roi <- list()
    
    if ( hasObjectHDF5(scanfile, 'Results') ) {
        
        # Get result info for HDF5 scan file.
        for ( phenotype in getResultPhenotypesHDF5(scanfile) ) {
            analyses <- getResultAnalysesHDF5(scanfile, phenotype)
            result.info[[phenotype]] <- lapply( analyses, function(analysis)
                getResultNamesHDF5(scanfile, phenotype, analysis) )
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
        stop("cannot output report - no relevant results found in file '",
            scanfile, "'")
    }
    
    # Introduce plot device.
    plot.device <- NULL
    
    # Set figure margin defaults (values taken from R par() doc).
    fig.margin.defaults <- c(1.02, 0.82, 0.82, 0.42)
    
    # Set figure margin width in inches.
    fig.margin.width <- sum(fig.margin.defaults[c(2, 4)])
    
    # Set figure margin height in inches.
    fig.margin.height <- sum(fig.margin.defaults[c(1, 3)])
    
    # Set page margin size in inches, assuming a 1-inch margin on either side.
    page.margins <- 2.0  # 2 x 1 inch margins
    
    # Set page width in inches (landscape A4).
    page.width <- 11.693
    
    # Adjust plot width to take account of margins.   
    plot.width <- page.width - page.margins - fig.margin.width
    
    # Get plot height: divide plot width by square root of 2.
    plot.height <- plot.width * 0.7071
    
    # Set figure dimensions.
    fig.width <- plot.width + fig.margin.width
    fig.height <- plot.height + fig.margin.height
    
    # Ensure graphics device will be switched off.
    on.exit( if ( ! is.null(plot.device) ) { dev.off(plot.device) } )
    
    # Create temp report file, ensure will be removed.
    tmp <- tempfile()
    on.exit( file.remove(tmp), add=TRUE )
    
    # Load Cairo PDF graphics device.
    grDevices::cairo_pdf(tmp, width=fig.width, height=fig.height, onefile=TRUE)
    plot.device <- dev.cur()
    
    # If Cairo PDF device failed to load, fall back on regular PDF device.
    if ( names(plot.device) == 'null device' ) {
        grDevices::pdf(tmp, width=fig.width, height=fig.height, onefile=TRUE)
        plot.device <- dev.cur()
    }
    
    # Check graphics device loaded successfully.
    if ( names(plot.device) == 'null device' ) {
        plot.device <- NULL
        stop("failed to load graphics device")
    }
    
    # Write output for each phenotype.
    for ( i in seq_along(result.info) ) {
        
        phenotype <- names(result.info)[i]
        
        if ( 'Scanone' %in% names(roi) ) {
            
            if ( 'Scanone' %in% names(result.info[[phenotype]]) &&
                'Result' %in% result.info[[phenotype]][['Scanone']] ) {
                
                scanone.result <- readResultHDF5(scanfile,
                    phenotype, 'Scanone', 'Result')
                
                # Get any QTL intervals.
                qtl.intervals <- readResultHDF5(scanfile,
                    phenotype, 'Scanone', 'QTL Intervals')
                
                # If no QTL intervals, get scanone threshold for this phenotype.
                # NB: if no scanone threshold found, threshold object will be NULL.
                if ( is.null(qtl.intervals) ) {
                    scanone.threshold <- readResultHDF5(scanfile,
                        phenotype, 'Scanone', 'Threshold')
                } else {
                    scanone.threshold <- NULL
                }
                
                # Plot (zero or more) QTL intervals across all sequences.
                plotQTLScanone(scanone.result, qtl.intervals=qtl.intervals,
                    threshold=scanone.threshold, phenotype=phenotype)
                
                # If significant QTL intervals found, plot
                # all sequences with a significant QTL.
                if ( length(qtl.intervals) > 0 ) {
                    
                    interval.seqs <- unique( sapply( qtl.intervals,
                        function(x) unique(x[, 'chr']) ) )
                    
                    for ( interval.seq in interval.seqs ) {
                        plotQTLScanone(scanone.result, qtl.intervals=qtl.intervals,
                            threshold=scanone.threshold, chr=interval.seq,
                            phenotype=phenotype)
                    }
                }
            }
        }
    }
    
    # Switch off graphics device.
    dev.off(plot.device)
    plot.device <- NULL
    
    # Move temp file to final report file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, report, overwrite=TRUE)
    
    return( invisible() )
}

# End of pdf.R #################################################################