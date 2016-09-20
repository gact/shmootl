# Start of run_makecross.R #####################################################

# run_makecross ----------------------------------------------------------------
#' Make cross from genotype and phenotype data.
#' 
#' Read separate genotype and phenotype data files,
#' and combine these into a single cross data file.
#' 
#' @param genfile input genotype CSV file
#' @param phefile input phenotype CSV file
#' @param crossfile output cross CSV file
#' 
#' @concept shmootl:utilities
#' @export
#' @family pipeline functions
#' @rdname run_makecross
run_makecross <- function(genfile, phefile, crossfile) {
    
    geno <- readGenoCSV(genfile)
    
    pheno <- readPhenoCSV(phefile)
    
    cross <- makeCross(geno, pheno)
    
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    writeCrossCSV(cross, tmp)
    
    # Move temp file to final cross file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, crossfile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_makecross.R #######################################################