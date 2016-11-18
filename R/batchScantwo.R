# Start of batchScantwo.R ######################################################

# batchPermScantwo -------------------------------------------------------------
#' Run \code{qtl::scantwo} on a batch of permuted \code{cross} objects.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Phenotype columns for which QTL analysis should be run. 
#' If no phenotypes are specified, all are used.
#' @template param-n.cluster
#' @param iseed Seed for random number generator.
#' @param n.perm Number of permutations.
#' @param perm.pheno Permute phenotype data.
#' @param perm.geno Permute genotype data.
#' @param ... Additional keyword arguments passed to \code{scantwo}.
#'     
#' @return A \code{scantwoperm} list containing the results of the QTL scan for 
#' all permutations.
#' 
#' @template author-thomas-walsh
#' @template author-yue-hu
#' @template ref-broman-2003
#' 
#' @export
#' @family scan utility functions
#' @rdname batchPermScantwo
batchPermScantwo <- function(cross, pheno.col=NULL, n.cluster=1, iseed=NULL, 
    n.perm=1000, perm.pheno=TRUE, perm.geno=FALSE, ...) {
    
    stopifnot( isSinglePositiveWholeNumber(n.perm) )
    stopifnot( allKwargs(...) )
    
    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)
    
    # Set batch scan arguments.
    args <- list(x=1:n.perm, scanfunction=nodePermScantwo, cross=cross, 
        pheno.col=pheno.col, n.cluster=n.cluster, iseed=iseed, 
        perm.pheno=perm.pheno, perm.geno=perm.geno)
    
    # Run permutation batch scan.
    scantwo.perms <- do.call(batchScan, c(args, list(...)))
    
    # Get model names from first permutation. 
    perm.models <- names(scantwo.perms[[1]])
    
    # Combine permutation results by model.
    combined.result <- list()
    for ( perm.model in perm.models ) {
        model.perms <- lapply(scantwo.perms, function(x) x[[perm.model]])
        combined.result[[perm.model]] <- do.call(rbind, model.perms)
        rownames(combined.result[[perm.model]]) <- NULL
    }
    
    # Set class of combined result.
    class(combined.result) <- c('scantwoperm', 'list')
    
    # Set attributes of combined result from those of first permutation.
    for ( a in const$scan.attributes[['scantwoperm']] ) {
        attr(combined.result, a) <- attr(scantwo.perms[[1]], a)
    }
    
    return(combined.result)
}

# batchPhenoScantwo ------------------------------------------------------------
#' Run \code{qtl::scantwo} on a batch of phenotypes.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Phenotype columns for which QTL analysis should be run. 
#' If no phenotypes are specified, all are used.
#' @template param-n.cluster
#' @param iseed Seed for random number generator.
#' @param ... Additional keyword arguments passed to \code{scantwo}.
#'     
#' @return A \code{scantwo} object containing the results of the QTL scan for 
#' the given phenotypes.
#' 
#' @template ref-broman-2003
#' 
#' @export
#' @family scan utility functions
#' @rdname batchPhenoScantwo
batchPhenoScantwo <- function(cross, pheno.col=NULL, n.cluster=1, iseed=NULL, 
    ...) {
    
    stopifnot( allKwargs(...) )
    
    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)

    # Set batch scan arguments.
    args <- list(x=pheno.col, scanfunction=nodePhenoScantwo, cross=cross, 
        n.cluster=n.cluster, iseed=iseed)

    # Run per-phenotype batch scan.
    scantwo.results <- do.call(batchScan, c(args, list(...)))

    # If multiple phenotypes, convert LOD matrix to array..
    if ( length(scantwo.results) > 1 ) {
        lod.list <- lapply(scantwo.results, function(scantwo.result) 
            scantwo.result$lod)
        lod.data <- array(unlist(lod.list), dim=c( dim(lod.list[[1]]), 
            length(lod.list) ) )
    } else { # ..otherwise take LOD matrix directly.
        lod.data <- scantwo.results[[1]]$lod
    }
    
    # Create combined scantwo result.
    combined.result <- list(lod=lod.data, map=scantwo.results[[1]]$map, 
        scanoneX=scantwo.results[[1]]$scanoneX)
    attr(combined.result, 'method') <- attr(scantwo.results[[1]], 'method')
    attr(combined.result, 'type') <- attr(scantwo.results[[1]], 'type')
    attr(combined.result, 'fullmap') <- attr(scantwo.results[[1]], 'fullmap')
    attr(combined.result, 'phenotypes') <- sapply(scantwo.results, 
        function(scantwo.result) attr(scantwo.result, 'phenotypes') )
    class(combined.result) <- 'scantwo'
    
    return(combined.result)
}

# nodePermScantwo --------------------------------------------------------------
#' Run \code{qtl::scantwo} on a single permuted \code{cross} object.
#' 
#' @param perm.id Permutation index.
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Phenotype columns for which QTL analysis should be run. 
#' If no phenotypes are specified, all are used.
#' @param perm.pheno Permute phenotype data.
#' @param perm.geno Permute genotype data.
#' @param ... Additional keyword arguments passed to \code{scantwo}.
#'     
#' @return A \code{scantwoperm} matrix containing the result of the QTL scan for
#' a single permutation.
#' 
#' @template ref-broman-2003
#' 
#' @export
#' @family scan utility functions
#' @rdname nodePermScantwo
nodePermScantwo <- function(perm.id, cross, pheno.col=NULL, perm.pheno=TRUE, 
    perm.geno=FALSE, ...) {

    stopifnot( isSinglePositiveWholeNumber(perm.id) )
    stopifnot( allKwargs(...) )
    
    kwargs <- list(...)
    
    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)
    
    # Get list of known qtl::scantwo arguments.
    known.args <- const$scan.args[['qtl::scantwo']]
    
    # Set vector of qtl::scantwo arguments that would cause problems here.
    unsupported.args <- c('batchsize', 'n.cluster', 'n.perm', 'perm.strata', 
        'perm.Xsp', 'verbose')
    
    unknown <- names(kwargs)[ ! names(kwargs) %in% known.args ]
    if ( length(unknown) > 0 ) {
        stop("unknown qtl::scantwo arguments passed to nodePermScantwo - '", toString(unknown), "'")
    }
    
    unsupported <- names(kwargs)[ names(kwargs) %in% unsupported.args ]
    if ( length(unsupported) > 0 ) {
        stop("unsupported qtl::scantwo arguments passed to nodePermScantwo - '", toString(unsupported), "'")
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
    
    # Set scan arguments.
    args <- list(cross=cross, pheno.col=pheno.col)
    
    # Run permutation scan.
    scantwo.result <- do.call(qtl::scantwo, c(args, kwargs))
    
    # Get phenotype names from permutation result.
    perm.pheno <- attr(scantwo.result, 'phenotypes')
    
    # Get scantwo result in permutation form.
    perm.result <- as.list( qtl::subrousummaryscantwo(scantwo.result, 
        for.perm=TRUE) )
    
    # Create scantwo permutation result.
    for ( perm.model in names(perm.result) ) {
        perm.result[[perm.model]] <- matrix( perm.result[[perm.model]], nrow=1, 
            byrow=TRUE, dimnames=list(perm.id, perm.pheno) )
    }
    class(perm.result) <- c('scantwoperm', 'list')
    for ( a in const$scan.attributes[['scantwoperm']] ) {
        attr(perm.result, a) <- attr(scantwo.result, a)
    }
    
    return(perm.result)
}

# nodePhenoScantwo -------------------------------------------------------------
#' Run \code{qtl::scantwo} on a single phenotype of a \code{cross} object.
#'  
#' @param pheno.col Phenotype column for which QTL analysis should be run.
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param ... Additional keyword arguments passed to \code{scantwo}.
#'     
#' @return A \code{scantwo} object containing the result of the QTL scan for 
#' a single phenotype. 
#' 
#' @template ref-broman-2003
#' 
#' @export
#' @family scan utility functions
#' @rdname nodePhenoScantwo
nodePhenoScantwo <- function(pheno.col, cross, ...) {
    
    stopifnot( allKwargs(...) )
    
    kwargs <- list(...)
    
    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)
    
    if ( length(pheno.col) != 1 ) {
        stop("nodePhenoScantwo cannot process multiple phenotypes")
    }
    
    # Get list of known qtl::scantwo arguments.
    known.args <- const$scan.args[['qtl::scantwo']]
    
    # Set vector of qtl::scantwo arguments that would cause problems here.
    unsupported.args <- c('batchsize', 'n.cluster', 'n.perm', 'perm.strata', 
        'perm.Xsp', 'verbose')
    
    unknown <- names(kwargs)[ ! names(kwargs) %in% known.args ]
    if ( length(unknown) > 0 ) {
        stop("unknown scantwo arguments passed to nodePhenoScantwo - '", toString(unknown), "'")
    }
    
    unsupported <- names(kwargs)[ names(kwargs) %in% unsupported.args ]
    if ( length(unsupported) > 0 ) {
        stop("unsupported scantwo arguments passed to nodePhenoScantwo - '", toString(unsupported), "'")
    }    
    
    # Set scan arguments.
    args <- list(cross=cross, pheno.col=pheno.col)
    
    # Run single-phenotype scan.
    scantwo.result <- do.call(qtl::scantwo, c(args, kwargs))
    
    return(scantwo.result)
}

# End of batchScantwo.R ########################################################