# Start of batchScanone.R ######################################################

# batchPermScanone -------------------------------------------------------------
#' Run \code{qtl::scanone} on a batch of permuted \code{cross} objects.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Phenotype columns for which QTL analysis should be run.
#' If no phenotypes are specified, all are used.
#' @template param-n.cluster
#' @param iseed Seed for random number generator.
#' @param n.perm Number of permutations.
#' @param perm.pheno Permute phenotype data.
#' @param perm.geno Permute genotype data.
#' @param perm.type Type of permutation data (see below).
#' @param ... Additional keyword arguments passed to \code{scanone}.
#' 
#' @return If \code{perm.type} is set to \code{'max'}, a regular
#' \code{scanoneperm} object is returned, containing the maximum
#' LOD values from the given permutations. If \code{perm.type} is
#' set to \code{'bins'}, a \code{scanonebins} array is returned.
#' Each row of this array corresponds to a permutation, each array
#' column corresponds to a bin spanning an interval of LOD values,
#' and each slice corresponds to a LOD column. Each element contains
#' the number of loci in the given bin interval for that LOD column
#' in that permutation.
#' 
#' @template ref-broman-2003
#' 
#' @export
#' @family scan utility functions
#' @importFrom abind abind
#' @rdname batchPermScanone
batchPermScanone <- function(cross, pheno.col=NULL, n.cluster=1, iseed=NULL, 
    n.perm=1000, perm.pheno=TRUE, perm.geno=FALSE, perm.type=c('max', 'bins'), ...) {
    
    stopifnot( isSinglePositiveWholeNumber(n.perm) )
    stopifnot( allKwargs(...) )
    
    perm.type <- match.arg(perm.type)
    
    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)
    
    # Set batch scan arguments.
    args <- list(x=1:n.perm, scanfunction=nodePermScanone, cross=cross, 
        pheno.col=pheno.col, n.cluster=n.cluster, iseed=iseed, 
        perm.pheno=perm.pheno, perm.geno=perm.geno, perm.type=perm.type)
    
    # Run permutation batch scan.
    scanone.perms <- do.call(batchScan, c(args, list(...)))
    
    # Combine permutation results.
    if ( perm.type == 'max' ) {
        
        combined.result <- do.call(rbind, scanone.perms)
        class(combined.result) <- c('scanoneperm', 'matrix')
        
    } else if ( perm.type == 'bins' ) {
        
        num.bins <- max( sapply(scanone.perms, function(x) dim(x)[2]) )
        scanone.perms <- lapply(scanone.perms, padBins, num.bins)
        combined.result <- abind::abind(scanone.perms, along=1)
        class(combined.result) <- c('scanonebins', 'array')
    }
    
    return(combined.result)
}

# batchPhenoScanone ------------------------------------------------------------
#' Run \code{qtl::scanone} on a batch of phenotypes.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Phenotype columns for which QTL analysis should be run.
#' If no phenotypes are specified, all are used.
#' @template param-n.cluster
#' @param iseed Seed for random number generator.
#' @param ... Additional keyword arguments passed to \code{scanone}.
#'     
#' @return A \code{scanone} object containing the results of \code{scanone} 
#' for the given phenotypes.
#' 
#' @template ref-broman-2003
#' 
#' @export
#' @family scan utility functions
#' @importFrom qtl phenames
#' @rdname batchPhenoScanone
batchPhenoScanone <- function(cross, pheno.col=NULL, n.cluster=1, iseed=NULL, ...) {
    
    stopifnot( allKwargs(...) )

    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)

    # Set batch scan arguments.
    args <- list(x=pheno.col, scanfunction=nodePhenoScanone, cross=cross, 
        n.cluster=n.cluster, iseed=iseed)
    
    # Run per-phenotype batch scan.
    scanone.results <- do.call(batchScan, c(args, list(...)))
    
    # If multiple LOD columns, assign phenotype names to LOD columns..
    if ( length(scanone.results) > 1 ) {
        
        pheno.names <- qtl::phenames(cross)[pheno.col]
        for ( i in seq_along(scanone.results) ) {
            lodcol.index = getDatColIndices(scanone.results[[i]])
            stopifnot( length(lodcol.index) == 1 )
            colnames(scanone.results[[i]])[lodcol.index] <- pheno.names[i]
        }
        
        # Combine per-phenotype scan results.
        combined.result <- do.call(cbind, scanone.results)
        
        # Set attributes of combined result from those of first result.
        for ( a in const$scan.attributes[['scanone']] ) {
            attr(combined.result, a) <- attr(scanone.results[[1]], a)
        }
        
    } else { # ..otherwise assign default LOD column name.
        
        combined.result <- scanone.results[[1]]
        lodcol.index = getDatColIndices(combined.result)
        stopifnot( length(lodcol.index) == 1 )
        colnames(combined.result)[lodcol.index] <- 'lod'
    }
    
    return(combined.result)
}

# nodePermScanone --------------------------------------------------------------
#' Run \code{qtl::scanone} on a single permuted \code{cross} object.
#' 
#' @param perm.id Permutation index.
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Phenotype columns for which QTL analysis should be run. 
#' If no phenotypes are specified, all are used.
#' @param perm.pheno Permute phenotype data.
#' @param perm.geno Permute genotype data.
#' @param perm.type Type of permutation data to return (see below).
#' @param ... Additional keyword arguments passed to \code{scanone}.
#' 
#' @return If \code{perm.type} is set to \code{'max'}, a regular
#' \code{scanoneperm} object is returned, containing the maximum
#' LOD value from the given permutation. If \code{perm.type} is
#' set to \code{'bins'}, a \code{scanonebins} array is returned.
#' The single row of this array corresponds to the permutation,
#' each column corresponds to a bin spanning an interval of LOD
#' values, and each slice corresponds to a LOD column. Each array
#' element contains the number of loci in the given bin interval
#' for that LOD column in that permutation.
#' 
#' @template ref-broman-2003
#' 
#' @export
#' @family scan utility functions
#' @rdname nodePermScanone
nodePermScanone <- function(perm.id, cross, pheno.col=NULL, perm.pheno=TRUE, 
    perm.geno=FALSE, perm.type=c('max', 'bins'), ...) {
    
    stopifnot( isSinglePositiveWholeNumber(perm.id) )
    stopifnot( allKwargs(...) )
    
    perm.type <- match.arg(perm.type)

    kwargs <- list(...)
    
    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)
    
    # Get list of known qtl::scanone arguments.
    known.args <- const$scan.args[['qtl::scanone']]

    # Set vector of qtl::scanone arguments that would cause problems here.
    unsupported.args <- c('batchsize', 'n.cluster', 'n.perm', 'perm.strata', 
        'perm.Xsp', 'verbose')
    
    unknown <- names(kwargs)[ ! names(kwargs) %in% known.args ]
    if ( length(unknown) > 0 ) {
        stop("unknown qtl::scanone arguments passed to nodePermScanone - '", toString(unknown), "'")
    }

    unsupported <- names(kwargs)[ names(kwargs) %in% unsupported.args ]
    if ( length(unsupported) > 0 ) {
        stop("unsupported qtl::scanone arguments passed to nodePermScanone - '", toString(unsupported), "'")
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
    
    if ( perm.type == 'max' ) {
        
        # Get max LOD value for each LOD column.
        perm.result <- apply(scanone.result[, lodcol.indices, drop=FALSE], 2, max)
        
        # Create scanone permutation result.
        perm.result <- matrix(perm.result, nrow=1, byrow=TRUE,
            dimnames=list(perm.id, names(perm.result)))
        
        # Set class of permutation result.
        class(perm.result) <- c('scanoneperm', 'matrix')
        
    } else if ( perm.type == 'bins' ) {
        
        # Bin LOD values of scanone permutation result.
        perm.result <- binLODValues(scanone.result)
        
        # Set rowname to permutation ID.
        dimnames(perm.result)[[1]] <- perm.id
    }
    
    return(perm.result)
}

# nodePhenoScanone -------------------------------------------------------------
#' Run \code{qtl::scanone} on a single phenotype of a \code{cross} object.
#'  
#' @param pheno.col Phenotype column for which QTL analysis should be run.
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param ... Additional keyword arguments passed to \code{scanone}.
#'     
#' @return A \code{scanone} object containing the result of \code{scanone}
#' for a single phenotype. 
#' 
#' @template ref-broman-2003
#' 
#' @export
#' @family scan utility functions
#' @rdname nodePhenoScanone
nodePhenoScanone <- function(pheno.col, cross, ...) {
    
    stopifnot( allKwargs(...) )
    
    kwargs <- list(...)
    
    # Get phenotype column indices.
    pheno.col <- getPhenoColIndices(cross, pheno.col)
    
    if ( length(pheno.col) != 1 ) {
        stop("nodePhenoScanone cannot process multiple phenotypes")
    } 
    
    # Get list of known qtl::scanone arguments.
    known.args <- const$scan.args[['qtl::scanone']]
    
    # Set vector of qtl::scanone arguments that would cause problems here.
    unsupported.args <- c('batchsize', 'n.cluster', 'n.perm', 'perm.strata', 
        'perm.Xsp', 'verbose')
    
    unknown <- names(kwargs)[ ! names(kwargs) %in% known.args ]
    if ( length(unknown) > 0 ) {
        stop("unknown scanone arguments passed to nodePhenoScanone - '", toString(unknown), "'")
    }
    
    unsupported <- names(kwargs)[ names(kwargs) %in% unsupported.args ]
    if ( length(unsupported) > 0 ) {
        stop("unsupported scanone arguments passed to nodePhenoScanone - '", toString(unsupported), "'")
    }
    
    # Set scan arguments.
    args <- list(cross=cross, pheno.col=pheno.col)
    
    # Run single-phenotype scan.
    scanone.result <- do.call(qtl::scanone, c(args, kwargs))
    
    return(scanone.result)
}

# End of batchScanone.R ########################################################