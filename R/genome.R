# Start of genome.R ############################################################

# genomeOpt --------------------------------------------------------------------
#' Set/get genome option.
#' 
#' This function returns the current value of the \pkg{shmootl} genome option.
#' If a genome is specified, the \pkg{shmootl} genome option is updated with the
#' given value. Call the function \code{\link{getGenomes}} to see available
#' genomes.
#' 
#' @param value New value of genome option.
#' 
#' @return Current value of genome option.
#' 
#' @export
#' @family genome functions
#' @family package option functions
#' @rdname genomeOpt
genomeOpt <- function(value) {
    
    result <- getOption('shmootl.genome', default=const$default$genome)
    
    if ( ! result %in% names(const$seqtab) ) {
        stop("unknown genome - '", result, "'")
    }
    
    if ( ! missing(value) && ! is.null(value) ) {
       
        if ( ! value %in% names(const$seqtab) ) {
            stop("unknown genome - '", value, "'")
        }
        
        options(shmootl.genome=value)
    }
    
    return(result)
}

# getGenomes -------------------------------------------------------------------
#' Get \pkg{shmootl} genome names.
#' 
#' @return Character vector of \pkg{shmootl} genome names.
#' 
#' @export
#' @family genome functions
#' @rdname getGenomes
getGenomes <- function() {
    return( names(const$seqtab) )
}

# End of genome.R ##############################################################