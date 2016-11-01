# Start of run_makecross.R #####################################################

# run_makecross ----------------------------------------------------------------
#' Make cross CSV file from genotype and phenotype data.
#' 
#' Read separate genotype and phenotype data files,
#' and combine these into a single cross data file.
#' 
#' @param genfile input genotype CSV file [required]
#' @param phefile input phenotype CSV file [required]
#' @param crossfile output cross CSV file [required]
#' 
#' @concept shmootl:preparation
#' @export
#' @family pipeline functions
#' @rdname run_makecross
run_makecross <- function(genfile=NA_character_, phefile=NA_character_,
    crossfile=NA_character_) {
    
    stopifnot( isSingleString(genfile) )
    stopifnot( isSingleString(phefile) )
    stopifnot( isSingleString(crossfile) )
    
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