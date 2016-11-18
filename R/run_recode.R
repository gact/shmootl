# Start of run_recode.R ########################################################

# run_recode -------------------------------------------------------------------
#' Recode data in \pkg{R/qtl} CSV file.
#' 
#' Recode data in an \pkg{R/qtl} cross or genotype CSV file. Genotypes can be
#' recoded by passing to the \code{geno} parameter a mapping of old to new
#' genotypes. Alternatively, genotype data can be converted to enumerated
#' genotypes with the \code{enum.geno} parameter.
#' 
#' When calling this function from within the \code{R} environment, the
#' \code{geno} parameter must be specified as a mapping object (e.g.
#' \code{mapping( c(A = 'W', B = 'S') )}). When called from the command line
#' using \code{Rscript}, the \code{geno} parameter must be specified as a a
#' YAML string (or YAML file) mapping old to new genotypes
#' (e.g. \code{"A: W, B: S"}).
#' 
#' @param datafile cross/geno CSV file [required]
#' @param geno recode genotypes from mapping
#' @param enum.geno recode to enumerated genotypes
#' 
#' @concept shmootl:utilities
#' @export
#' @family pipeline functions
#' @rdname run_recode
run_recode <- function(datafile=NA_character_, geno=mapping(),
    enum.geno=FALSE) {
    
    stopifnot( isSingleString(datafile) )
    
    geno <- if ( ! identical(geno, NA_character_) ) { geno } else { NULL }
    
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