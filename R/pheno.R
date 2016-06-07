# Start of pheno.R #############################################################

# as.data.frame.pheno ----------------------------------------------------------
#' Coerce \code{pheno} object to \code{data.frame}.
#' 
#' @param x A \code{pheno} object.
#' @param ... Unused arguments.
#' @param digits If specified, round numeric phenotype values to the given
#' number of digits.
#' 
#' @return A \code{data.frame} corresponding to the input \code{pheno} object.
#' 
#' @keywords internal
#' @method as.data.frame pheno
#' @rdname as.data.frame.pheno
as.data.frame.pheno <- function(x, ..., digits=NULL) {
    
    # Get CrossInfo object.
    cross.info <- attr(x, 'info')
    
    compareCrossInfo(x, cross.info)
    
    # Get mask of actual phenotype columns.
    pheno.mask <- ! tolower(colnames(x)) %in%
        const$reserved.phenotypes
    
    # If digits specified, round numeric phenotype values.
    if ( ! is.null(digits) ) {
        stopifnot( isSinglePositiveWholeNumber(digits) )
        numeric.mask <- sapply(unname(x), is.numeric)
        mask <- pheno.mask & numeric.mask
        x[, mask] <- round(x[, mask], digits=digits)
    }
    
    # Get phenotype headings from pheno object,
    # set actual phenotype IDs from CrossInfo.
    phenotypes <- colnames(x)
    phenotypes[pheno.mask] <- getPhenotypes(cross.info)
    
    # Coerce phenotype data to type character.
    x <- coerceDataFrame(x, rep('character', ncol(x)))
    
    # Bind phenotype headings and phenotype data.
    pheno.frame <- as.data.frame(rbind(phenotypes, x),
        stringsAsFactors=FALSE)
    
    return(pheno.frame)
}

# as.pheno ---------------------------------------------------------------------
#' Coerce to \code{pheno} object.
#' 
#' @param from Object containing phenotype data.
#' 
#' @return A \code{pheno} object corresponding to the input object.
#' 
#' @keywords internal
#' @rdname as.pheno
as.pheno <- function(from) {
    UseMethod('as.pheno', from)
}

# as.pheno.data.frame ----------------------------------------------------------
#' @method as.pheno data.frame
#' @rdname as.pheno
as.pheno.data.frame <- function(from) {
    
    # Get indices of initial row.
    head.row <- 1
    
    # Get index of first and last data rows.
    first.data.row <- head.row + 1
    last.data.row <- nrow(from)
    
    stopifnot( last.data.row >= first.data.row )
    
    # Get vector of data row indices.
    dat.rows <- first.data.row : last.data.row
    
    # Get phenotype headings.
    phenotypes <- as.character(from[head.row, ])
    
    # Get index of phenotype columns with R/qtl special heading 'id'. If present,
    # this contains identifiers of sampled individuals. Search is case-insensitive.
    id.col <- which( tolower(phenotypes) == 'id' )
    
    if ( length(id.col) > 0 ) {
        # If 'id' column present, set vector of
        # sample IDs and remove from phenotypes..
        stopifnot( length(id.col) == 1 )
        samples <- as.character(from[dat.rows, id.col])
        phenotypes <- phenotypes[-id.col]
        
    } else {
        # ..otherwise set vector of sample indices.
        samples <- getIndices(dat.rows)
    }
    
    # Check for blank phenotypes.
    if ( any(from[dat.rows, -id.col] == '') ) {
        stop("blank phenotype values found")
    }
    
    # Get phenotype data.
    cross.pheno <- do.call(cbind.data.frame,
        lapply(from[dat.rows, ], type.convert))
    
    # Create CrossInfo object.
    cross.info <- methods::new('CrossInfo')
    cross.info <- setPhenotypes(cross.info, phenotypes)
    cross.info <- setSamples(cross.info, samples)
    
    attr(cross.pheno, 'info') <- cross.info
    
    class(cross.pheno) <- c('pheno', 'data.frame')
    
    return(cross.pheno)
}

# End of pheno.R ###############################################################