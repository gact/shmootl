# Start of compareCrossInfo.R ##################################################

# compareCrossInfo -------------------------------------------------------------
#' Compare object with associated \code{CrossInfo}.
#'
#' @description Test concordance of an object with its associated 
#' \code{CrossInfo} object. If no \code{CrossInfo} object is specified,
#' it will be extracted from the test object.
#' 
#' @param x Test object to be compared.
#' @template param-CrossInfo
#' 
#' @return TRUE if \code{CrossInfo} and test object are concordant; otherwise
#' returns a character vector of mismatch errors.
#' 
#' @include CrossInfo-class.R
#' @keywords internal
#' @rdname compareCrossInfo
#' @seealso \code{\linkS4class{CrossInfo}}
compareCrossInfo <- function (x, cross.info=NULL) {
    UseMethod('compareCrossInfo', x)
}

# compareCrossInfo.cross -------------------------------------------------------
#' @rdname compareCrossInfo
compareCrossInfo.cross <- function(x, cross.info=NULL) {
    
    # If CrossInfo object not specified, 
    # extract it from the cross itself.
    if ( is.null(cross.info) ) {
        cross.info <- attr(x, 'info')
    }
    
    if ( is.null(cross.info) ) {
        stop("CrossInfo not found in cross object")
    }
    
    validObject(cross.info)
    
    errors <- vector('character')
    
    cross.map <- qtl::pull.map(x)
    cross.markers <- pullLocusIDs(cross.map)
    obj.markers <- getMarkers(cross.info)
    
    if ( any(obj.markers != cross.markers) ) {
        errors <- c(errors, "marker mismatch")
    }
    
    if ( length(obj.markers) > 0 ) {
        
        cross.locus.seqs <- pullLocusSeq(cross.map)
        obj.locus.seqs <- getMarkerSeqs(cross.info)
        
        if ( any(obj.locus.seqs != cross.locus.seqs) ) {
            errors <- c(errors, "marker sequence mismatch")
        }
    }
    
    pheno.col <- getPhenoColIndices(x)
    cross.phenames <- colnames(x$pheno)[pheno.col]
    obj.phenames <- getPhenotypeNames(cross.info)
    
    if ( any(obj.phenames != cross.phenames) ) {
        errors <- c(errors, "phenotype mismatch")
    }
    
    id.col <- getIdColIndex(x)
    
    if ( getNumSamples(cross.info) != nrow(x$pheno) ) {
        errors <- c(errors, "sample count mismatch")
    }
    
    if ( xor( hasSampleIDs(cross.info), ! is.null(id.col) ) ) {
        errors <- c(errors, "sample ID presence/absence mismatch")
    }
    
    if ( hasSampleIDs(cross.info) ) {
        
        cross.samples <- as.character(x$pheno[, id.col])
        obj.samples <- getSamples(cross.info)
        
        if ( any(obj.samples != cross.samples) ) {
            errors <- c(errors, "sample ID mismatch")
        }
    }
    
    cross.alleles <- pull.alleles(x)
    obj.alleles <- cross.info@alleles
    if ( length(obj.alleles) != length(cross.alleles) ||
         any(obj.alleles != cross.alleles) ) {
        errors <- c(errors, "genotype/allele mismatch")
    }
    
    return( if ( length(errors) == 0 ) {TRUE} else {errors} )
}

# End of compareCrossInfo.R ####################################################