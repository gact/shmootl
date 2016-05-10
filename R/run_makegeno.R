# Start of run_makegeno.R ######################################################

# run_makegeno -----------------------------------------------------------------
#' Run run_makegeno pipeline.
#' 
#' @description This script will read genotype data from the specified VCF 
#' files, then write that genotype data to an \pkg{R/qtl} genotype CSV file.
#' 
#' @param samples sample VCF file
#' @param founders optional founder VCF file
#' @param genfile output genotype CSV file
#' @param digits numeric precision [default: unrounded]
#' 
#' @export
#' @rdname run_makegeno
run_makegeno <- function(samples, genfile, founders=NA, digits=NA) {
    
    stopifnot( isSingleString(genfile) )
    
    digits <- if ( ! is.na(digits) ) { strtoi(digits) } else { NULL }
    
    sample.ids <- getSamplesVCF(samples)
    
    if ( ! is.na(founders) ) {
        founder.ids <- getSamplesVCF(founders)
        infiles <- c(samples, founders)
    } else {
        founder.ids <- NULL
        infiles <- samples
    }
    
    geno <- readGenoVCF(infiles, samples=sample.ids, founders=founder.ids)
    
    writeGenoCSV(geno, genfile, digits=digits)
    
    return( invisible() )
}

# End of run_makegeno.R ########################################################