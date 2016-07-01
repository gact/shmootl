# Start of run_makecross.R #####################################################

# run_makecross ----------------------------------------------------------------
#' Make cross from geno and pheno data.
#' 
#' @description This script combines separate genotype and 
#' phenotype files into a single cross data file.
#' 
#' @param genfile input genotype CSV file
#' @param phefile input phenotype CSV file
#' @param crossfile output cross CSV file
#' 
#' @export
#' @rdname run_makecross
run_makecross <- function(genfile, phefile, crossfile) {
    
    geno <- readGenoCSV(genfile)
    
    pheno <- readPhenoCSV(phefile)
    
    cross <- makeCross(geno, pheno)
    
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    writeCrossCSV(cross, tmp)
    
    file.copy(tmp, crossfile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_makecross.R #######################################################