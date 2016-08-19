# Start of run_makegeno.R ######################################################

# run_makegeno -----------------------------------------------------------------
#' Make genotype data from VCF.
#' 
#' @description This script will read genotype data from the specified VCF 
#' files, then write that genotype data to an \pkg{R/qtl} genotype CSV file.
#' 
#' @param datafile sample VCF file
#' @param fdrfile optional founder VCF file
#' @param genfile output genotype CSV file
#' @param alleles founder allele symbols
#' @param digits numeric precision [default: unrounded]
#' 
#' @concept shmootl:utilities
#' @export
#' @rdname run_makegeno
run_makegeno <- function(datafile, genfile, fdrfile=NA_character_,
    alleles=character(), digits=NA_integer_) {
    
    stopifnot( isSingleString(genfile) )
    
    alleles <- if ( ! identical(alleles, character()) ) { alleles } else { NULL }
    digits <- if ( ! identical(digits, NA_integer_) ) { digits } else { NULL }
    
    sample.ids <- getSamplesVCF(datafile)
    
    if ( ! is.na(fdrfile) ) {
        founder.ids <- getSamplesVCF(fdrfile)
        infiles <- c(datafile, fdrfile)
    } else {
        founder.ids <- NULL
        infiles <- datafile
    }
    
    geno <- readGenoVCF(infiles, samples=sample.ids,
        founders=founder.ids, alleles=alleles)
    
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    writeGenoCSV(geno, tmp, digits=digits)
    
    # Move temp file to final genotypes file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, genfile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_makegeno.R ########################################################