# Start of compareCrossInfo.R ##################################################

# compareCrossInfo -------------------------------------------------------------
#' Compare \code{cross} and \code{CrossInfo} objects.
#'
#' @description Test concordance of an \pkg{R/qtl} \code{cross} and a CrossInfo 
#' object. If an \pkg{R/qtl} \code{cross} is the first argument, it should have 
#' attribute \code{'info'} containing a \code{CrossInfo} object. If a 
#' \code{CrossInfo} object is the first argument, the second argument should be 
#' the \pkg{R/qtl} \code{cross} object with which it is being compared.
#' 
#' @param ... Arguments (see below).
#' @template param-CrossInfo
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' 
#' @return TRUE if \code{CrossInfo} and \pkg{R/qtl} \code{cross} objects are 
#' concordant; otherwise returns a character vector of mismatch errors.
#' 
#' @include CrossInfo-class.R
#' @keywords internal
#' @rdname compareCrossInfo
#' @seealso \code{\linkS4class{CrossInfo}}
compareCrossInfo <- function (...) {
    UseMethod('compareCrossInfo', list(...)[[1]])
}

# compareCrossInfo.cross -------------------------------------------------------
#' @rdname compareCrossInfo
compareCrossInfo.cross <- function(cross) {
    cross.info <- attr(cross, 'info')
    stopifnot( 'CrossInfo' %in% class(cross.info) )
    return( compareCrossInfo(cross.info, cross) )
}

# compareCrossInfo.CrossInfo ---------------------------------------------------
#' @rdname compareCrossInfo
compareCrossInfo.CrossInfo <- function (cross.info, cross) {
    
    stopifnot( 'cross' %in% class(cross) )
    validObject(cross.info)
    
    errors <- vector('character')
    
    cross.map <- qtl::pull.map(cross)
    cross.markers <- pullLocusIDs(cross.map)
    obj.markers <- getMarkerNames(cross.info)
    
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
    
    pheno.col <- getPhenoColIndices(cross)
    cross.phenames <- colnames(cross$pheno)[pheno.col]
    obj.phenames <- getPhenotypeNames(cross.info)
    
    if ( any(obj.phenames != cross.phenames) ) {
        errors <- c(errors, "phenotype mismatch")
    }
    
    id.col <- getIdColIndex(cross)
    
    if ( getNumSamples(cross.info) != nrow(cross$pheno) ) {
        errors <- c(errors, "sample count mismatch")
    }
    
    if ( xor( hasSampleIDs(cross.info), ! is.null(id.col) ) ) {
        errors <- c(errors, "sample ID presence/absence mismatch")
    }
    
    if ( hasSampleIDs(cross.info) ) {
        
        cross.samples <- as.character(cross$pheno[, id.col])
        obj.samples <- getSamples(cross.info)
        
        if ( any(obj.samples != cross.samples) ) {
            errors <- c(errors, "sample ID mismatch")
        }
    }
    
    cross.alleles <- pull.alleles(cross)
    obj.alleles <- cross.info@alleles
    if ( length(obj.alleles) != length(cross.alleles) ||
         any(obj.alleles != cross.alleles) ) {
        errors <- c(errors, "genotype/allele mismatch")
    }

    return( if ( length(errors) == 0 ) {TRUE} else {errors} )
}

# compareCrossInfo (S4) --------------------------------------------------------
#' @rdname compareCrossInfo
setGeneric('compareCrossInfo', compareCrossInfo)

# CrossInfo::compareCrossInfo --------------------------------------------------
#' @rdname compareCrossInfo
setMethod('compareCrossInfo', signature='CrossInfo', 
    definition = compareCrossInfo.CrossInfo)

# End of compareCrossInfo.R ####################################################