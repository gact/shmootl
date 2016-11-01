# Start of zzz.R ###############################################################

# .onLoad ----------------------------------------------------------------------
#' @include const.R
#' @include map.R
.onLoad <- function(libname, pkgname) {
    
    # Setup package options ----------------------------------------------------
    
    options(shmootl.genome=const$default$genome)
    
    # Setup remaining package constants ----------------------------------------
    
    setupDefaultMapkeys()
    
    const$excel$digest$Overview[['headings']] <- c('File', 'Phenotype',
        getPkgAnalysisNames())
    
    # Lock package constants ---------------------------------------------------
    
    lockEnvironment(const, bindings=TRUE)
    
    # --------------------------------------------------------------------------
}

# End of zzz.R #################################################################