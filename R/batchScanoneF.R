# Start of batchScanoneF.R #####################################################

# batchPermScanoneF ------------------------------------------------------------
#' Run \code{funqtl::scanoneF} on a batch of permuted \code{cross} objects.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Phenotype columns for which QTL analysis should be 
#' performed. If specified, phenotypes should contain measurements for 
#' consecutive values of the parameter of the function-valued trait 
#' (e.g. times). If no phenotypes are specified, all are used.
#' @template param-n.cluster
#' @param iseed Seed for random number generator.
#' @param n.perm Number of permutations.
#' @param perm.pheno Permute phenotype data.
#' @param perm.geno Permute genotype data.
#' @param ... Additional keyword arguments passed to \code{scanoneF}.
#'     
#' @return A \code{scanoneperm} matrix containing the results of \code{scanoneF}
#' for all permutations.
#' 
#' @export
#' @family scan utilities
#' @rdname batchPermScanoneF
batchPermScanoneF <- function(cross, pheno.col=NULL, n.cluster=1, iseed=NULL, 
    n.perm=1000, perm.pheno=TRUE, perm.geno=FALSE, ...) {
    
    stopifnot( isSinglePositiveWholeNumber(n.perm) )
    stopifnot( allKwargs(...) )

    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)
    
    # Set batch scan arguments.
    args <- list(x=1:n.perm, scanfunction=nodePermScanoneF, cross=cross, 
        pheno.col=pheno.col, n.cluster=n.cluster, iseed=iseed, 
        perm.pheno=perm.pheno, perm.geno=perm.geno)
    
    # Run permutation batch scan.
    scanone.perms <- do.call(batchScan, c(args, list(...)))
    
    # Combine permutation results.
    scanone.result <- do.call(rbind, scanone.perms)
    
    # Set class of combined result.
    class(scanone.result) <- c('scanoneperm', 'matrix')
    
    return(scanone.result)
}

# nodePermScanoneF -------------------------------------------------------------
#' Run \code{funqtl::scanoneF} on a single permutation of a \code{cross} object.
#' 
#' @param perm.id Permutation index.
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Phenotype columns for which QTL analysis should be 
#' performed. If specified, phenotypes should contain measurements for 
#' consecutive values of the parameter of the function-valued trait 
#' (e.g. times). If no phenotypes are specified, all are used.
#' @param perm.pheno Permute phenotype data.
#' @param perm.geno Permute genotype data.
#' @param ... Additional keyword arguments passed to \code{scanoneF}.
#'     
#' @return A \code{scanoneperm} matrix containing the results of \code{scanoneF}
#' for a single permutation.
#' 
#' @export
#' @family scan utilities
#' @rdname nodePermScanoneF
nodePermScanoneF <- function(perm.id, cross, pheno.col=NULL, perm.pheno=TRUE, 
    perm.geno=FALSE, ...) {
    
    stopifnot( isSinglePositiveWholeNumber(perm.id) )
    stopifnot( allKwargs(...) )
    
    kwargs <- list(...)
    
    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)
    
    # Get list of known funqtl::scanoneF arguments.
    known.args <- union(const$scan.args[['funqtl::scanoneF']], 
        const$scan.args[['qtl::scanone']])
    
    # Set vector of funqtl::scanoneF arguments that would cause problems here.
    unsupported.args <- c('batchsize', 'n.cluster', 'n.perm', 'perm.strata', 
        'perm.Xsp', 'pheno.cols', 'verbose')
    
    unknown <- names(kwargs)[ ! names(kwargs) %in% known.args ]
    if ( length(unknown) > 0 ) {
        stop("unknown funqtl::scanoneF arguments passed to nodePermScanoneF - '", toString(unknown), "'")
    }
    
    unsupported <- names(kwargs)[ names(kwargs) %in% unsupported.args ]
    if ( length(unsupported) > 0 ) {
        stop("unsupported funqtl::scanoneF arguments passed to nodePermScanoneF - '", toString(unsupported), "'")
    }
    
    # Generate permutation indices for cross object.
    perm.indices <- permIndices(cross)
    
    # Permute cross data.
    cross <- permCross(cross, perm.indices=perm.indices, 
        perm.pheno=perm.pheno, perm.geno=perm.geno)
    
    # If permuting phenotypes, permute any corresponding data in the same way.
    if (perm.pheno) {
        
        if ( ! is.null(kwargs[['addcovar']]) ) {
            kwargs[['addcovar']] <- kwargs[['addcovar']][perm.indices, ]
        }
        
        if ( ! is.null(kwargs[['intcovar']]) ) {
            kwargs[['intcovar']] <- kwargs[['intcovar']][perm.indices, ]
        }
        
        if ( ! is.null(kwargs[['weights']]) ) {
            kwargs[['weights']] <- kwargs[['weights']][perm.indices]
        }        
        
        if ( ! is.null(kwargs[['ind.noqtl']]) ) {
            kwargs[['ind.noqtl']] <- kwargs[['ind.noqtl']][perm.indices]
        }
    }
    
    # Set scan arguments.
    args <- list(cross=cross, pheno.col=pheno.col)
    
    # Run permutation scan.
    scanone.result <- do.call(qtl::scanone, c(args, kwargs))
    
    # Get LOD column indices.
    lodcol.indices = getDatColIndices(scanone.result)
    
    # Set SLOD values from the mean LOD value across phenotypes at each locus.
    scanone.slod <- rowMeans(scanone.result[, lodcol.indices, drop=FALSE])
    
    # Set MLOD values from the maximum LOD value across phenotypes at each locus.
    scanone.mlod <- apply(scanone.result[, lodcol.indices, drop=FALSE], 1, max)
    
    # Create scanoneF permutation result.
    perm.result <- c( 
        slod=max(scanone.slod), 
        mlod=max(scanone.mlod)
    )
    
    # Convert permutation result to matrix.
    perm.result <- matrix(perm.result, nrow=1, byrow=TRUE,
        dimnames=list(perm.id, names(perm.result)))
    
    # Set class of permutation result.
    class(perm.result) <- c('scanoneperm', 'matrix')
    
    return(perm.result)
}

# End of batchScanoneF.R #######################################################