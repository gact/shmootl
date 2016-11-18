# Start of pdf.R ###############################################################

# plotReportTitlePagePDF -------------------------------------------------------
#' Plot a QTL analysis report title page.
#' 
#' @param scanfile File path of HDF5 scan file
#' from which the report is being generated.
#' 
#' @keywords internal
#' @rdname plotReportTitlePagePDF
plotReportTitlePagePDF <- function(scanfile) {
    
    opar <- graphics::par(no.readonly=TRUE)
    on.exit( graphics::par(opar), add=TRUE )
    
    stopifnot( isSingleString(scanfile) )
    stopifnot( file.exists(scanfile) )
    
    graphics::par( mar=rep(4.1, 4) )
    
    xlim <- ylim <- c(0.0, 1.0)
    
    plot(NA, bg='white', bty='n', fg='black', type='n', xaxt='n', xlab='',
        xlim=xlim, yaxt='n', ylab='', ylim=0:1)
    
    text(0.5, 0.75, 'QTL Analysis Report', adj=0.5,
        cex=2.0, family='sans', font=2)
    
    ipar <- list(cex=1, family='mono', font=1)
    
    line.spacing <- 2 * strheight(' ', cex=ipar$cex,
        family=ipar$family, font=ipar$font)
    
    midline <- 0.5
    
    usr <- graphics::par('usr')
    rect(usr[1], midline - 4*line.spacing, usr[2],
        midline + 1*line.spacing, lty=1, lwd=3)
    
    plot.width <- diff(xlim)
    char.width <- strwidth(' ', cex=ipar$cex, family=ipar$family, font=ipar$font)
    max.nchar <- floor( plot.width / char.width )
    
    max.file.nchar <- max.nchar - nchar('File: ')
    scanfile <- ellipt(scanfile, max.file.nchar, left=FALSE, right=FALSE)
    
    file.line <- paste0('File: ', scanfile)
    text(0.0, midline - 1*line.spacing, file.line, adj=0,
        cex=ipar$cex, family=ipar$family, font=ipar$font)
    
    date.line <- paste('Date:', Sys.Date())
    text(0.0, midline - 2*line.spacing, date.line, adj=0,
        cex=ipar$cex, family=ipar$family, font=ipar$font)
    
    return( invisible() )
}

# writeReportPDF ---------------------------------------------------------------
#' Write a PDF report of QTL scan results.
#'  
#' @param scanfile A scan result HDF5 file.
#' @param report Path of output report PDF file.
#' @param phenotypes Phenotypes for which results should be included in the
#' report file. If none are specified, results are output for all phenotypes.
#' @param analyses Analyses for which results should be included in the report
#' file. If none are specified, results are output for all available analyses.
#' 
#' @export
#' @importFrom grDevices cairo_pdf
#' @importFrom grDevices dev.cur
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @rdname writeReportPDF
writeReportPDF <- function(scanfile, report, phenotypes=NULL, analyses=NULL) {
    
    stopifnot( isSingleString(scanfile) )
    stopifnot( file.exists(scanfile) )
    stopifnot( isSingleString(report) )
    
    # Set possible results to be sought in scan file.
    results.sought <- supported.results <- list(
        'Scanone' = c('Result')
    )
    
    # If analyses specified, filter results sought by given analyses.
    if ( ! is.null(analyses) ) {
        
        analyses <- unique( resolveAnalysisTitle(analyses) )
        
        unsupported <- analyses[ ! analyses %in% names(supported.results) ]
        if ( length(unsupported) > 0 ) {
            stop("PDF output not supported for analyses - '",
                toString(unsupported), "'")
        }
        
        results.sought <- results.sought[ names(results.sought) %in% analyses ]
    }
    
    # Get result info.
    rinfo <- getResultInfoHDF5(scanfile, phenotypes=phenotypes,
        analyses=names(results.sought))
    
    # Get results of interest.
    roi <- list()
    for ( phenotype in names(rinfo[[scanfile]]) ) {
        for ( analysis in names(rinfo[[scanfile]][[phenotype]]) ) {
            for ( result in rinfo[[scanfile]][[phenotype]][[analysis]] ) {
                if ( analysis %in% names(results.sought) &&
                    result %in% results.sought[[analysis]] ) {
                    roi[[analysis]] <- union(roi[[analysis]], result)
                }
            }
        }
    }
    
    if ( length(roi) == 0 ) {
        stop("no relevant results found in file '", scanfile, "'")
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
        stop("failed to load PDF graphics device")
    }
    
    # Output report title page.
    plotReportTitlePagePDF(scanfile)
    
    # Write output for each phenotype.
    for ( i in seq_along(rinfo[[scanfile]]) ) {
        
        phenotype <- names(rinfo[[scanfile]])[i]
        
        if ( 'Scanone' %in% names(roi) ) {
            
            if ( 'Scanone' %in% names(rinfo[[scanfile]][[phenotype]]) &&
                'Result' %in% rinfo[[scanfile]][[phenotype]][['Scanone']] ) {
                
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