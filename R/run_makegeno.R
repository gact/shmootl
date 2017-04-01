# Start of run_makegeno.R ######################################################

# run_makegeno -----------------------------------------------------------------
#' Make genotype data from VCF.
#' 
#' Read raw genotype data from the sample VCF file, encode the given genotype
#' data as founder or enumerated genotypes, and write the encoded genotypes
#' to an \pkg{R/qtl} genotype CSV file.
#' 
#' If no founder VCF file is specified, markers are assigned enumerated
#' genotypes, so that each raw allele is converted to a number in the order
#' in which it is first observed at a given marker (i.e. \code{'1'},
#' \code{'2'}, etc.).
#' 
#' If a founder VCF file is given, markers are assigned a genotype symbol
#' consisting of alleles that each correspond to a specific founder.
#' 
#' If the \code{alleles} parameter is specified, this must be a mapping of
#' founder sample IDs to allele symbols. If calling this function from within
#' the \code{R} environment, this must be specified as a mapping object (e.g.
#' \code{mapping( c(DBVPG6044 = 'W', Y12 = 'S') )}). When called from the
#' command line using \code{Rscript}, the \code{alleles} parameter must be
#' specified as a YAML string (or YAML file) mapping founders to allele
#' symbols (e.g. \code{"DBVPG6044: W, Y12: S"}). If the \code{alleles}
#' parameter is not specified, allele symbols are taken from the letters
#' of the alphabet (i.e. \code{'A'}, \code{'B'} etc.).
#' 
#' @param datafile sample VCF file [required]
#' @param fdrfile optional founder VCF file
#' @param genfile output genotype CSV file [required]
#' @param alleles founder allele symbol mapping
#' @param digits numeric precision [default: unrounded]
#' 
#' @concept shmootl:preparation
#' @export
#' @family pipeline functions
#' @rdname run_makegeno
run_makegeno <- function(datafile=NA_character_, genfile=NA_character_,
    fdrfile=NA_character_, alleles=mapping(), digits=NA_integer_) {
    
    stopifnot( isSingleString(genfile) )
    
    alleles <- if ( ! identical(alleles, mapping()) ) { alleles } else { NULL }
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