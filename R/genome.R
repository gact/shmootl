# Start of genome.R ############################################################

# genomeOpt --------------------------------------------------------------------
#' Set/get genome option.
#' 
#' If a genome is specified, the \pkg{shmootl} genome option is updated with the 
#' given value. Otherwise, this function returns the current value of the genome
#' option.
#' 
#' @param genome New value of genome option.
#' 
#' @return Current value of genome option.
#' 
#' @export
#' @rdname genomeOpt
genomeOpt <- function(genome) {
    
    if ( ! missing(genome) ) {
       
        if ( ! genome %in% names(const$seqinfo) ) {
            stop("unknown genome - '", genome, "'")
        }
        
        options(shmootl.genome=genome)
        
        result <- invisible()
        
    } else {
        
        genome <- getOption('shmootl.genome', default=const$default$genome)
        
        if ( ! genome %in% names(const$seqinfo) ) {
            stop("unknown genome - '", genome, "'")
        }
        
        result <- genome
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