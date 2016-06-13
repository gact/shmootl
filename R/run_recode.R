# Start of run_recode.R ########################################################

# run_recode -------------------------------------------------------------------
#' Recode data in \pkg{R/qtl} CSV file.
#' 
#' @param infile input CSV file
#' @param outfile output CSV file
#' @param geno genotype recode mapping
#' 
#' @export
#' @rdname run_recode
run_recode <- function(infile, outfile, geno=NA) {
    
    if ( ! is.na(geno) ) {
        geno <- unlist( loadMappingFromLine(geno) )
    } else {
        geno <- NULL
    }
    
    recodeCSV(infile, outfile, geno=geno)
    
    return( invisible() )
}

# End of run_recode.R ##########################################################