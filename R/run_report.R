# Start of run_report.R ########################################################

# run_report -------------------------------------------------------------------
#' Create report of \pkg{R/qtl} output.
#' 
#' @param resultfile HDF5 result file
#' @param reportfile PDF report file
#' 
#' @export
#' @importFrom grDevices pdf
#' @rdname run_report
run_report <- function(resultfile, reportfile) {
    
    stopifnot( isSingleString(resultfile) )
    stopifnot( file.exists(resultfile) )
    stopifnot( isSingleString(reportfile) )
    
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
    on.exit({ graphics.off() })
    
    # Remove existing report file.
    if ( file.exists(reportfile) ) {
        file.remove(reportfile)
    }
    
    # Init PDF graphics device.  
    grDevices::pdf(reportfile, width=fig.width, height=fig.height)
    
    # Get phenotypes from result file.
    results <- readGroupMemberNamesHDF5(resultfile, 'Results')
    phenotypes <- results[ results != 'Overview' ]
    
    # Write output for each phenotype.
    for ( i in getIndices(phenotypes) ) {
        
        phenotype <- phenotypes[i]
        
        # Get QTL intervals.
        qtl.intervals <- readResultHDF5(resultfile, phenotype, 'QTL Intervals')
        
        # Get scanone result.
        pheno.result <- readResultHDF5(resultfile, phenotype, 'Scanone')
        
        # Plot (zero or more) QTL intervals across all sequences.
        plotQTLIntervals(qtl.intervals, pheno.result, phenotype=phenotype)
        
        # If significant QTL intervals found, plot
        # all sequences with a significant QTL.
        if ( length(qtl.intervals) > 0 ) {
            
            interval.seqs <- unique( sapply( qtl.intervals,
                function(x) unique(x[, 'chr']) ) )
            
            for ( interval.seq in interval.seqs ) {
                plotQTLIntervals(qtl.intervals, pheno.result,
                    chr=interval.seq, phenotype=phenotype)
            }
        }
    }
    
    return( invisible() )
}

# End of run_report.R ##########################################################