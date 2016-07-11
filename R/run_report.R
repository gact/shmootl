# Start of run_report.R ########################################################

# run_report -------------------------------------------------------------------
#' Create report of \pkg{R/qtl} analysis results.
#' 
#' @param scanfile scan result HDF5 file
#' @param report scan report file
#' 
#' @export
#' @importFrom grDevices cairo_pdf
#' @importFrom grDevices dev.cur
#' @importFrom grDevices dev.off
#' @rdname run_report
run_report <- function(scanfile, report) {
    
    stopifnot( isSingleString(scanfile) )
    stopifnot( file.exists(scanfile) )
    stopifnot( isSingleString(report) )
    
    results.sought <- 'Scanone'
    results.found <- list()
    result.info <- list()
    
    if ( hasObjectHDF5(scanfile, 'Results') ) {
        
        # Get phenotypes from scan result file.
        phenotypes <- getPhenotypesHDF5(scanfile)
        
        # Get result info for phenotypes.
        result.info <- lapply(phenotypes, function(phenotype)
            getResultNamesHDF5(scanfile, phenotype))
        names(result.info) <- phenotypes
        
        # Get results of interest for each phenotype.
        results.found <- lapply(result.info, function(results)
            results[ results %in% results.sought ])
    }
    
    # Get set of results of interest.
    roi <- unique( unlist(results.found) )
    
    if ( length(roi) == 0 ) {
        warning("cannot output report - results not found in file '", scanfile, "'")
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
    
    # Init PDF graphics device.
    grDevices::cairo_pdf(tmp, width=fig.width, height=fig.height, onefile=TRUE)
    plot.device <- dev.cur()
    
    # Write output for each phenotype.
    for ( i in getIndices(phenotypes) ) {
        
        phenotype <- phenotypes[i]
        
        if ( 'Scanone' %in% result.info[[phenotype]] ) {
            
            scanone.result <- readResultHDF5(scanfile, phenotype, 'Scanone')
            
            # Get any QTL intervals.
            qtl.intervals <- readResultHDF5(scanfile, phenotype, 'QTL Intervals')
            
            # If no QTL intervals, get scanone permutations for this phenotype,
            # and create a QTL intervals object with any threshold info.
            # NB: if no permutations found, threshold attributes will be NULL.
            if ( is.null(qtl.intervals) ) {
                scanone.perms <- readResultHDF5(scanfile, phenotype, 'Scanone Perms')
                qtl.intervals <- qtlintervals(threshold=attr(scanone.perms, 'threshold'),
                    alpha=attr(scanone.perms, 'alpha'), fdr=attr(scanone.perms, 'fdr'))
            }
            
            # Plot (zero or more) QTL intervals across all sequences.
            plotQTLScanone(scanone.result, qtl.intervals=qtl.intervals, phenotype=phenotype)
            
            # If significant QTL intervals found, plot
            # all sequences with a significant QTL.
            if ( length(qtl.intervals) > 0 ) {
                
                interval.seqs <- unique( sapply( qtl.intervals,
                    function(x) unique(x[, 'chr']) ) )
                
                for ( interval.seq in interval.seqs ) {
                    plotQTLScanone(scanone.result, qtl.intervals=qtl.intervals,
                        chr=interval.seq, phenotype=phenotype)
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

# End of run_report.R ##########################################################