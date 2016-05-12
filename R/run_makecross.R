# Start of run_makecross.R #####################################################

# run_makecross ----------------------------------------------------------------
#' Run makecross pipeline.
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
    
    writeCrossCSV(cross, crossfile)
    
    return( invisible() )
}

# End of run_makecross.R #######################################################