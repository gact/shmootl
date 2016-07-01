# Start of run_recode.R ########################################################

# run_recode -------------------------------------------------------------------
#' Recode data in \pkg{R/qtl} CSV file.
#' 
#' @param datafile cross/geno CSV file
#' @param geno genotype recode mapping
#' 
#' @export
#' @rdname run_recode
run_recode <- function(datafile, geno=NA) {
    
    if ( ! is.na(geno) ) {
        geno <- unlist( loadMappingFromLine(geno) )
    } else {
        geno <- NULL
    }
    
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    recodeCSV(datafile, tmp, geno=geno)
    
    file.copy(tmp, datafile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_recode.R ##########################################################