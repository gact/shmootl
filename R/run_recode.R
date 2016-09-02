# Start of run_recode.R ########################################################

# run_recode -------------------------------------------------------------------
#' Recode data in \pkg{R/qtl} CSV file.
#' 
#' @param datafile cross/geno CSV file
#' @param geno recode genotypes from mapping
#' @param enum.geno recode to enumerated genotypes
#' 
#' @concept shmootl:utilities
#' @export
#' @rdname run_recode
run_recode <- function(datafile, geno=mapping(), enum.geno=FALSE) {
    
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    recodeCSV(datafile, tmp, geno=geno, enum.geno=enum.geno)
    
    # Move temp file to final output file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, datafile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_recode.R ##########################################################