# Start of run_digest.R ########################################################

# run_digest -------------------------------------------------------------------
#' Create digest of QTL analysis results.
#' 
#' Given one or more scan result HDF5 files, create
#' a digest of the results of the QTL analyses.
#' 
#' @param h5list list of HDF5 scan files [required]
#' @param digest scan digest file [required]
#' @param scanfile.pattern scan file name pattern
#' 
#' @concept shmootl:processing
#' @export
#' @family pipeline functions
#' @importFrom tools file_ext
#' @rdname run_digest
run_digest <- function(h5list=character(), digest=NA_character_,
    scanfile.pattern=NA_character_) {
    
    stopifnot( is.character(h5list) )
    stopifnot( length(h5list) > 0 )
    stopifnot( all( file.exists(h5list) ) )
    stopifnot( isSingleString(digest) )
    
    if ( identical(scanfile.pattern, NA_character_) ) {
        scanfile.pattern <- NULL
    }
    
    # Get digest file extension.
    digest.ext <- tools::file_ext(digest)
    
    # Create temp digest file, ensure will be removed.
    # NB: must have correct extension so that Excel format can be set automatically.
    tmp <- tempfile( fileext=paste0('.', digest.ext) ) 
    on.exit( file.remove(tmp), add=TRUE )
    
    # Write digest to temp file.
    success <- tryCatch({
        
        if ( digest.ext %in% const$ext$excel ) {
            writeDigestExcel(h5list, tmp, scanfile.pattern=scanfile.pattern)
        } else {
            stop("cannot create digest - unknown extension on file '", digest, "'")
        }
        
        result <- TRUE
        
    }, error=function(e) {
        result <- FALSE
    })
    
    # If digest written without error, move temp file to final digest file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    if (success) {
        file.copy(tmp, digest, overwrite=TRUE)
    }
    
    return( invisible() )
}

# End of run_digest.R ##########################################################