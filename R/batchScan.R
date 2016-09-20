# Start of scan.R ##############################################################

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
#' @family scan utility functions
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
    
    # If running on local host, run cluster locally..
    if ( all(nodes == 'localhost') ) 
    {
        # If running on a Unix-alike platform, run a fork cluster..
        if ( .Platform$OS.type == 'unix' )
        {
            cl <- parallel::makeForkCluster( length(nodes) )
        }
        # ..otherwise run a parallel socket cluster.
        else
        {
            cl <- parallel::makePSOCKcluster( length(nodes) )
        }
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