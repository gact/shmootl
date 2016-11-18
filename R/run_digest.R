# Start of run_digest.R ########################################################

# run_digest -------------------------------------------------------------------
#' Create digest of QTL analysis results.
#' 
#' Given one or more scan result HDF5 files, create a digest of the results of
#' the QTL analyses. Results can be output for specific phenotypes, analyses,
#' and worksheets. These output constraints are applied with parameters
#' \code{pheno}, \code{scans}, and \code{sheets}, respectively.
#' 
#' @param h5list list of HDF5 scan files [required]
#' @param digest scan digest file [required]
#' @param pheno phenotypes to output [default: all]
#' @param scans analyses to output [default: all]
#' @param sheets worksheets to output
#' @param scanfile.pattern scan file name pattern
#' 
#' @concept shmootl:processing
#' @export
#' @family pipeline functions
#' @importFrom tools file_ext
#' @rdname run_digest
run_digest <- function(h5list=character(), digest=NA_character_,
    pheno=character(), scans=character(), sheets=character(),
    scanfile.pattern=NA_character_) {
    
    stopifnot( is.character(h5list) )
    stopifnot( length(h5list) > 0 )
    stopifnot( all( file.exists(h5list) ) )
    stopifnot( isSingleString(digest) )
    
    pheno <- if ( ! identical(pheno, character()) ) { pheno } else { NULL }
    scans <- if ( ! identical(scans, character()) ) { scans } else { NULL }
    sheets <- if ( ! identical(sheets, character()) ) { sheets } else { NULL }
    scanfile.pattern <- if ( ! identical(scanfile.pattern, NA_character_) ) {
        scanfile.pattern } else { NULL }
    
    # Get digest file extension.
    digest.ext <- tools::file_ext(digest)
    
    # Infer digest file format.
    digest.format <- inferFormatFromFilename(digest)
    
    # Create temp digest file, ensure will be removed.
    # NB: must have correct extension so that Excel format can be set automatically.
    tmp <- tempfile( fileext=paste0('.', digest.ext) ) 
    on.exit( file.remove(tmp), add=TRUE )
    
    # Write digest to temp file.
    if ( digest.format == 'Excel' ) {
        writeDigestExcel(h5list, tmp, phenotypes=pheno, analyses=scans,
            worksheets=sheets, scanfile.pattern=scanfile.pattern)
    } else {
        stop("digest output not supported for ", digest.format, " format")
    }
    
    # If digest written without error, move temp file to final digest file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, digest, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_digest.R ##########################################################