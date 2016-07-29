# Start of run_digest.R ########################################################

# run_digest -------------------------------------------------------------------
#' Create digest of QTL analysis results.
#' 
#' @param h5list list of scan result files
#' @param digest scan digest file
#' 
#' @export
#' @importFrom tools file_ext
#' @include const.R
#' @rdname run_digest
run_digest <- function(h5list=character(), digest=NA_character_) {
    
    stopifnot( is.character(h5list) )
    stopifnot( length(h5list) > 0 )
    stopifnot( all( file.exists(h5list) ) )
    stopifnot( isSingleString(digest) )
    
    # Get digest file extension.
    digest.ext <- tools::file_ext(digest)
    
    # Create temp digest file, ensure will be removed.
    # NB: must have correct extension so that Excel format can be set automatically.
    tmp <- tempfile( fileext=paste0('.', digest.ext) ) 
    on.exit( file.remove(tmp), add=TRUE )
    
    if ( digest.ext %in% const$ext$excel ) {
        writeDigestExcel(h5list, tmp)
    } else {
        stop("cannot create digest - unknown extension on file '", digest, "'")
    }
    
    # Move temp file to final digest file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, digest, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_digest.R ##########################################################