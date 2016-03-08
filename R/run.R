# Start of run.R ###############################################################

# TODO: pipelines 'estimap' and 'makegeno'

# run --------------------------------------------------------------------------
#' Run shmootl pipeline from the command line.
#'   
#' @description This is the general \pkg{shmootl} pipeline-running function, 
#' which provides a basic interface for standard QTL analysis pipelines to be 
#' run from a command line. 
#' 
#' @section Pipeline Functions:
#' 
#' Pipelines are fully-documented \pkg{shmootl} package functions that follow 
#' specific conventions regarding function name and formal arguments. They are
#' intended to be run from the command line, but it should also be possible to
#' run them from within the R environment. 
#' 
#' Pipeline function names are of the form \code{run_<pipeline>} (where 
#' \code{<pipeline>} is replaced the pipeline name). Information about a given 
#' pipeline is taken automatically from that pipeline's function definition 
#' and documentation. 
#'  
#' Pipeline function arguments are handled by \pkg{argparser}, which groups 
#' arguments into three main types:
#' 
#' \itemize{
#' \item{\strong{positional arguments}:}  {arguments identified by their position 
#' (e.g. \code{'input.csv'}).}
#' \item{\strong{optional arguments}:}  {arguments in keyword-value pairs, where
#' the argument value follows the keyword (e.g. \code{'--alpha 0.01'}).}
#' \item{\strong{flags}:}  {keyword-only arguments, where the presence of the keyword 
#' toggles a specific option (e.g. \code{'--help'}). }
#' }
#' 
#' Pipeline arguments are assigned to one of these three groups depending on
#' their default values in the pipeline function definition. Arguments without 
#' a default value are taken to be positional arguments. Those with a default 
#' value are taken to be an optional argument, unless the default value is 
#' \code{FALSE}, in which case the argument is assumed to be \code{logical} 
#' (i.e. either \code{TRUE} or \code{FALSE}) and is treated as a flag.
#' 
#' Each argument's data type is set from its default value, if possible. An 
#' argument will be taken to be of type \code{character} if it doesn't have a
#' default value from which its data type can be guessed. Pipeline functions
#' should be written to take such arguments as \code{character} variables.
#' 
#' @examples
#' \dontrun{
#' # To list available pipelines.
#' Rscript -e 'library(shmootl)' -e 'run()'
#' 
#' # To run the 'scanone' pipeline, which makes use of the qtl::scanone function.
#' Rscript -e 'library(shmootl)' -e 'run()' scanone --n.cluster 8 input.csv output.hdf5
#' }
#'   
#' @export
#' @rdname run
run <- function() {
    
    # Pipeline Notes -----------------------------------------------------------
    # Pipeline and argument names should be simple and informative, while being 
    # as consistent as possible with existing pipelines. For more understandable 
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
            
            # Create argument parser.
            ap <- getPipelineArgparser(pipeline)
            
            # Parse arguments.
            args <- argparser::parse_args(ap, args)
            
            # Remove argparser special arguments.
            args <- args[ ! names(args) %in% const$special.params ]
            
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