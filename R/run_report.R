# Start of run_report.R ########################################################

# run_report -------------------------------------------------------------------
#' Create report of \pkg{R/qtl} analysis results.
#' 
#' @param scanfile scan result HDF5 file
#' @param report scan report PDF file
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
    
    # Get phenotypes from result file.
    results <- readGroupMemberNamesHDF5(scanfile, 'Results')
    phenotypes <- results[ results != 'Overview' ]
    
    # Write output for each phenotype.
    for ( i in getIndices(phenotypes) ) {
        
        phenotype <- phenotypes[i]
        
        # Get any QTL intervals.
        qtl.intervals <- NULL
        tryCatch({
            qtl.intervals <- readResultHDF5(scanfile, phenotype, 'QTL Intervals')
        }, error=function(e) {})
        
        if ( is.null(qtl.intervals) ) {
            
            # Get scanone permutations for this phenotype. # TODO: simplify
            h5name <- joinH5ObjectNameParts( c('Results', phenotype) )
            pheno.results <- readGroupMemberNamesHDF5(scanfile, h5name)
            index <- pmatch('Scanone Perms', pheno.results)
            stopif( is.na(index) )
            pheno.perms <- readResultHDF5(scanfile, phenotype, pheno.results[index])
            
            # Get threshold parameters from permutations.
            threshold <- attr(pheno.perms, 'threshold')
            alpha <- attr(pheno.perms, 'alpha')
            fdr <- attr(pheno.perms, 'fdr')
            
            # Create empty QTL intervals object with threshold info.
            qtl.intervals <- qtlintervals(threshold=threshold, alpha=alpha, fdr=fdr) 
        }
        
        # Get scanone result.
        pheno.result <- readResultHDF5(scanfile, phenotype, 'Scanone')
        
        # Plot (zero or more) QTL intervals across all sequences.
        plotQTLScanone(pheno.result, qtl.intervals=qtl.intervals, phenotype=phenotype)
        
        # If significant QTL intervals found, plot
        # all sequences with a significant QTL.
        if ( length(qtl.intervals) > 0 ) {
            
            interval.seqs <- unique( sapply( qtl.intervals,
                function(x) unique(x[, 'chr']) ) )
            
            for ( interval.seq in interval.seqs ) {
                plotQTLScanone(pheno.result, qtl.intervals=qtl.intervals,
                    chr=interval.seq, phenotype=phenotype)
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