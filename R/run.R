# Start of run.R ###############################################################

# TODO: pipeline 'scantwo'
# TODO: pipeline 'scanonef'
# TODO: pipeline 'scantwof'

# run --------------------------------------------------------------------------
#' Run shmootl pipeline from the command line.
#'   
#' @description This is the general \pkg{shmootl} pipeline-running function,
#' which provides a basic interface for standard QTL analysis pipelines to be
#' run from a command line. It cannot be called interactively within R.
#' 
#' @template section-pipelines
#'   
#' @export
#' @rdname run
run <- function() {
    
    if ( interactive() ) {
        stop("shmootl::run cannot be called interactively")
    }
    
    if ( 'package:shmootl' %in% search() ) {
        
        # Get command-line arguments.
        args <- commandArgs(trailingOnly=TRUE)
        
        # Take pipeline name from first argument, 
        # then take arguments from remainder.
        pipeline <- args[1]
        args <- args[-1]
        
        # If shmootl is available and valid pipeline is specified, run.
        if ( ! is.na(pipeline) && pipeline %in% getPkgPipelineNames() ) {
            
            # Prep argument parser.
            ap <- prepPipelineArgparser(pipeline)
            
            # Parse pipeline arguments.
            args <- argparser::parse_args(ap, args)
            
            # Process pipeline arguments.
            args <- procPipelineArgs(ap, args)
            
            # Get pipeline function name.
            pipe.func.name <- getPipelineFunctionName(pipeline)
            
            # Run pipeline.
            do.call(pipe.func.name, args)
            
            return( invisible() )
        } 
    }
 
    # Output usage info to standard error.
    message( getPkgPipelineUsage() )

    return( invisible() )
}

# End of run.R #################################################################