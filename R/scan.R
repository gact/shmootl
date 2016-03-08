# Start of scan.R ##############################################################

# Scan Utilities ---------------------------------------------------------------
#' QTL scan utilities.
#' 
#' @description The most basic scan utility is \code{\link{batchScan}}, which 
#' performs a batch QTL analysis on an \pkg{R/qtl} \code{cross} with a 
#' user-specified scan function. In general, batch scan functions should 
#' call \code{batchScan}. Each batch scan function must pass a node scan 
#' function to \code{batchScan}: the node scan function should operate on 
#' a single element of the batch input vector, which might contain phenotypes 
#' or permutation indices. 
#' 
#' For example, \code{\link{batchPermScanone}} can be used to perform a batch of
#' permutation scans with \code{qtl::scanone}. This function generates a vector
#' of permutation indices from the specified number of permutations, and passes
#' node function \code{\link{nodePermScanone}} to \code{\link{batchScan}}, which 
#' applies the node function to the vector of permutation indices. The result is
#' a \code{scanoneperm} matrix with one row for each permutation.
#' 
#' As another example, \code{\link{batchPhenoScanone}} can be used to perform a 
#' multi-phenotype batch of \code{qtl::scanone} tasks. This function generates a 
#' vector of phenotype column indices, each of which is individually passed to
#' \code{\link{nodePermScanone}} by \code{\link{batchScan}}. Each individual 
#' phenotype is analysed in a separate batch task, and the output is a 
#' \code{scanone} object containing results for the given phenotypes.
#' 
#' @docType package
#' @name Scan Utilities
NULL

# batchScan --------------------------------------------------------------------
#' Perform a batch of QTL analysis.
#' 
#' The specified \code{scanfunction} is applied to the \code{cross} for every 
#' element in vector \code{x}. This function uses the \pkg{parallel} package to
#' set up a cluster, initialiase a random-number generator stream, and run the
#' batch in parallel. 
#'    
#' @param x An atomic vector to which the scan function is applied elementwise.
#' @param scanfunction QTL scan function. The first argument of this function 
#' must be the element of \code{x} that is passed by \code{clusterApply}, and 
#' its second argument should be the \code{cross} object.
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @template param-n.cluster
#' @param iseed Seed for random number generator stream.
#' @param ... Additional keyword arguments passed to the scan function.
#'     
#' @return List of QTL scan results.
#' 
#' @export
#' @family scan utilities
#' @rdname batchScan
batchScan <- function(x, scanfunction, cross, n.cluster=1, iseed=NULL, ...) {
    
    stopifnot( is.vector(x) )
    stopifnot( is.function(scanfunction) )
    stopifnot( 'cross' %in% class(cross) )
    stopifnot( isSinglePositiveWholeNumber(n.cluster) )
    stopifnot( is.null(iseed) || is.integer(iseed) )
    stopifnot( allKwargs(...) )
    
    # Reduce cluster request for small jobs.
    n.cluster <- min(length(x), n.cluster)
    
    # Request nodes.
    nodes <- requestNodes(n.cluster)
    
    # If running on local host, run fork cluster locally..
    if ( all(nodes == 'localhost') ) 
    {
        cl <- parallel::makeForkCluster( length(nodes) )
    } 
    # ..otherwise run in parallel cluster.
    else
    {
        cl <- parallel::makeCluster(nodes)
    }
    
    # Ensure cluster will stop on exit. 
    on.exit( parallel::stopCluster(cl) )
    
    # Set pseudo-random number stream with given seed.
    parallel::clusterSetRNGStream(cl, iseed)
    
    # Run batch scan.
    batch.results <- parallel::clusterApply(cl, x, scanfunction, cross, ...)
        
    return(batch.results)
}

# End of scan.R ################################################################