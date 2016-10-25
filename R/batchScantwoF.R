# Start of batchScantwoF.R #####################################################

# batchPermScantwoF ------------------------------------------------------------
#' Run \code{funqtl::scantwoF} on batch of permuted \code{cross} objects.
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
#' @param ... Additional keyword arguments passed to \code{scantwoF}.
#'     
#' @return A \code{scantwoperm} list containing the results of \code{scantwoF}
#' for all permutations. This has three elements - \code{'one'},
#' \code{'fullvadd'}, and \code{'fv1'} - each containing a matrix
#' of permutation results for one of three models: single QTL
#' analysis, full-versus-additive QTL model, and full-versus-single
#' QTL model, respectively.
#' 
#' @template ref-broman-2003
#' @template ref-kwak-2014
#' 
#' @export
#' @family scan utility functions
#' @rdname batchPermScantwoF
batchPermScantwoF <- function(cross, pheno.col=NULL, n.cluster=1, iseed=NULL,
    n.perm=1000, perm.pheno=TRUE, perm.geno=FALSE, ...) {
    
    stopifnot( isSinglePositiveWholeNumber(n.perm) )
    stopifnot( allKwargs(...) )
    
    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)
    
    # Set batch scan arguments.
    args <- list(x=1:n.perm, scanfunction=nodePermScantwoF, cross=cross, 
        pheno.col=pheno.col, n.cluster=n.cluster, iseed=iseed, 
        perm.pheno=perm.pheno, perm.geno=perm.geno)
    
    # Run permutation batch scan.
    scantwo.perms <- do.call(batchScan, c(args, list(...)))
    
    # Get model names from first permutation. 
    perm.models <- names(scantwo.perms[[1]])
    
    # Combine permutation results by model.
    combined.result <- list()
    for ( perm.model  in perm.models ) {
        model.perms <- lapply(scantwo.perms, function(x) x[[perm.model]])
        combined.result[[perm.model]] <- do.call(rbind, model.perms)
    }
    
    # Set class of combined result.
    class(combined.result) <- c('scantwoperm', 'list')
    
    return(combined.result)
}

# nodePermScantwoF -------------------------------------------------------------
#' Run \code{funqtl::scantwoF} on a single permutation of a \code{cross} object.
#' 
#' @param perm.id Permutation index.
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Phenotype columns for which QTL analysis should be 
#' performed. If specified, phenotypes should contain measurements for 
#' consecutive values of the parameter of the function-valued trait 
#' (e.g. times). If no phenotypes are specified, all are used.
#' @param perm.pheno Permute phenotype data.
#' @param perm.geno Permute genotype data.
#' @param ... Additional keyword arguments passed to \code{scantwoF}.
#'     
#' @return A \code{scantwoperm} list containing the results of \code{scantwoF}
#' for a single permutation. 
#' 
#' @template ref-broman-2003
#' @template ref-kwak-2014
#' 
#' @export
#' @family scan utility functions
#' @rdname nodePermScantwoF
nodePermScantwoF <- function(perm.id, cross, pheno.col=NULL, perm.pheno=TRUE, 
    perm.geno=FALSE, ...) {

    stopifnot( isSinglePositiveWholeNumber(perm.id) )
    stopifnot( allKwargs(...) )
        
    kwargs <- list(...)
    
    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)
        
    # Get list of known funqtl::scantwoF arguments.
    known.args <- union(const$scan.args[['funqtl::scantwoF']], 
        const$scan.args[['qtl::scantwo']])
    
    # Set vector of funqtl::scantwoF arguments that would cause problems here.
    unsupported.args <- c('batchsize', 'n.cluster', 'n.perm', 'perm.strata', 
        'perm.Xsp', 'pheno.cols', 'usec', 'verbose')
        
    unknown <- names(kwargs)[ ! names(kwargs) %in% known.args ]
    if ( length(unknown) > 0 ) {
        stop("unknown funqtl::scantwoF arguments passed to nodePermScantwoF - '", toString(unknown), "'")
    }
    
    unsupported <- names(kwargs)[ names(kwargs) %in% unsupported.args ]
    if ( length(unsupported) > 0 ) {
        stop("unsupported funqtl::scantwoF arguments passed to nodePermScantwoF - '", toString(unsupported), "'")
    }
    
    # Generate permutation indices for cross object.
    perm.indices <- permIndices(cross)
    
    # Permute cross data.
    cross <- permCross(cross, perm.indices=perm.indices, perm.pheno=perm.pheno, 
        perm.geno=perm.geno)
    
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
    }
    
    # Set scanone arguments. 
    args <- list(cross=cross, pheno.col=pheno.col)

    # Run permutation scanone.
    scanone.result <- do.call(qtl::scanone, c(args, kwargs))
    
    # Get LOD column indices for scanone result.
    lodcol.indices = getLodColIndices(scanone.result)
    
    # Set scanone SLOD values from the mean LOD value across phenotypes at each locus.
    scanone.slod <- rowMeans(scanone.result[, lodcol.indices, drop=FALSE])
    
    # Set scanone MLOD values from the maximum LOD value across phenotypes at each locus.
    scanone.mlod <- apply(scanone.result[, lodcol.indices, drop=FALSE], 1, max)
    
    # Run permutation scantwo.
    scantwo.result <- do.call(qtl::scantwo, c(args, kwargs))
    
    # Get scantwo SLOD values.
    scantwo.mean <- scantwo.result
    scantwo.mean$lod <- apply(scantwo.result$lod, 1:2, mean)
    mean.summary <- summary(scantwo.mean)
    
    # Get scantwo MLOD values.
    scantwo.max <- scantwo.result
    scantwo.max$lod <- apply(scantwo.result$lod, 1:2, max)
    max.summary <- summary(scantwo.max)
    
    # Create scantwoF permutation result.
    perm.result <- list(
        one = c( 
            Slods = max(scanone.slod),
            Mlods = max(scanone.mlod)
        ),
        fullvadd = c( 
            SlodsH = max(mean.summary$lod.int), 
            MlodsH = max(max.summary$lod.int)
        ),
        fv1 = c( 
            SlodsL = max(mean.summary$lod.fv1),
            MlodsL = max(max.summary$lod.fv1)
        )
    )
    
    # Convert permutation result to a matrix for each model.
    for ( perm.model in names(perm.result) ) {
        perm.result[[perm.model]] <- matrix( perm.result[[perm.model]], nrow=1, 
            byrow=TRUE, dimnames=list(perm.id, names(perm.result[[perm.model]])) )
    }
    
    # Set class of permutation result.
    class(perm.result) <- c('scantwoperm', 'list')
    
    return(perm.result)
}

# End of batchScantwoF.R #######################################################