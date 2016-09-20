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
    
    # Pipeline Notes -----------------------------------------------------------
    # Pipeline and argument names should be simple and informative, while being 
    # as consistent as possible with existing pipelines. For more well-formatted
    # argparser help output, try to have 6-12 characters in positional argument 
    # names, 6-9 characters in optional argument names, and up to 7 characters 
    # in flag names.
    
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
        if ( ! is.na(pipeline) && pipeline %in% getPipelines() ) {
            
            # Prep argument parser.
            ap <- prepPipelineArgparser(pipeline)
            
            # Parse pipeline arguments.
            args <- argparser::parse_args(ap, args)
            
            # Process pipeline arguments.
            args <- procPipelineArgs(ap, args)
            
            # Get pipeline function.
            pipe.func <- getPipelineFunction(pipeline)
            
            # Run pipeline.
            do.call(pipe.func, args)
            
            return( invisible() )
        } 
    }
 
    # Output usage info to standard error.
    message( getPipelineUsage() )

    return( invisible() )
}

# End of run.R #################################################################