# Start of genome.R ############################################################

# genomeOpt --------------------------------------------------------------------
#' Set/get genome option.
#' 
#' If a genome is specified, the \pkg{shmootl} genome option is updated with the 
#' given value. Otherwise, this function returns the current value of the genome
#' option.
#' 
#' @param value New value of genome option.
#' 
#' @return Current value of genome option.
#' 
#' @export
#' @rdname genomeOpt
genomeOpt <- function(value) {
    
    result <- getOption('shmootl.genome', default=const$default$genome)
    
    if ( ! result %in% names(const$seqinfo) ) {
        stop("unknown genome - '", result, "'")
    }
    
    if ( ! missing(value) && ! is.null(value) ) {
       
        if ( ! value %in% names(const$seqinfo) ) {
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
#' @keywords internal
#' @rdname getGenomes
getGenomes <- function() {
    return( names(const$seqinfo) )
}

# End of genome.R ##############################################################