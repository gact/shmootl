# Start of pipeline.R ##########################################################

# getPipelineArgparser ---------------------------------------------------------
#' Get argument parser for \pkg{shmootl} pipeline.
#' 
#' @param pipeline Name of a \pkg{shmootl} pipeline.
#' 
#' @return Argument parser for the specified \pkg{shmootl} pipeline.
#' 
#' @keywords internal
#' @rdname getPipelineArgparser
getPipelineArgparser <- function(pipeline) {
    
    # Get pipeline info.
    pipe.info <- getPipelineInfo(pipeline)
    
    # Compose argument parser description by combining 
    # the function title, description, and details.
    description <- paste0(pipe.info$title, '\n')
    for ( content in c(pipe.info$description, pipe.info$details) ) {
        if ( ! is.null(content) ) {
            description <- paste0(description, content)
        }
    }
    
    # Create argument parser.
    ap <- argparser::arg_parser(name=pipe.info$command, description=description)
    
    # Add any parameters to the argument parser.
    if ( ! is.null(pipe.info$params) ) {
        
        for ( i in 1:length(pipe.info$params$arg) ) {
            
            ap <- argparser::add_argument(ap, 
                arg     = pipe.info$params$arg[i],
                short   = pipe.info$params$short[i],
                flag    = (pipe.info$params$group[i] == 'flag'),
                default = pipe.info$params$default[[i]],
                help    = pipe.info$params$help[i]
            )
        }
    }
    
    return(ap)
}

# getPipelineFunction ----------------------------------------------------------
#' Get \pkg{shmootl} pipeline function.
#' 
#' @param pipeline Name of a \pkg{shmootl} pipeline.
#' 
#' @return Function for the given \pkg{shmootl} pipeline.
#' 
#' @keywords internal
#' @rdname getPipelineFunction
getPipelineFunction <- function(pipeline) {
    stopifnot( isSingleString(pipeline) )
    stopifnot( pipeline %in% getPipelines() )
    stopifnot( 'package:shmootl' %in% search() )
    return( get( paste0('run_', pipeline) ) )
}

# getPipelineInfo --------------------------------------------------------------
#' Get info for \pkg{shmootl} pipeline.
#' 
#' This function interrogates the function definition and documentation for the 
#' specified pipeline and returns the resulting information as a list. The 
#' elements of this list are as follows:
#' 
#' \itemize{
#' \item{\code{'command'}} {Shell command used to call pipeline.}
#' \item{\code{'title'}} {Pipeline function title.}
#' \item{\code{'description'}} {Pipeline function description. This is not 
#' returned if it duplicates the function title.}
#' \item{\code{'details'}} {Pipeline function details.}
#' \item{\code{'params'}} {A \code{data.frame} containing pipeline function 
#' parameter info. If the function has no parameters, this takes a \code{NULL} 
#' value.}
#' }
#'  
#' Elements \code{'command'} and \code{'title'} are guaranteed to contain their
#' respective values. The remaining elements are only returned if corresponding 
#' documentation is available; otherwise, they are assigned a \code{NULL} value. 
#' 
#' @param pipeline Name of \pkg{shmootl} pipeline.
#' 
#' @return A \code{list} containing info for the given \pkg{shmootl} pipeline.
#' 
#' @keywords internal
#' @rdname getPipelineInfo
getPipelineInfo <- function(pipeline) {
    
    stopifnot( isSingleString(pipeline) )
    
    if ( ! pipeline %in% getPipelines() ) {
        stop("shmootl pipeline not found - '", pipeline, "'")
    }
    
    # Set the shell command used to run this pipeline.
    command <- paste("Rscript -e 'shmootl::run()'", pipeline)
    
    # Compose pipeline function name.
    pipe.func <- paste0('run_', pipeline)
    
    # Get formal argument info for pipeline function.
    f <- formals(pipe.func)
    
    # Get defined parameters.
    def.params <- names(f)
    
    # Get number of defined parameters.
    n.params <- length(def.params)
    
    # If pipeline function has any defined 
    # parameters, assemble argument info.
    if ( ! is.null(def.params) ) {
        
        # Init parameter info.
        params <- list(
            arg = character(n.params),
            short = character(n.params),
            group = factor(character(n.params), levels=const$param.groups),
            default = vector('list', n.params),
            help = character(n.params)
        )
        
        for ( i in 1:length(params) ) {
            names(params[[i]]) <- def.params
        }
  
        # Get info for each defined parameter.
        for ( p in def.params ) {
            
            # If parameter has a default value, 
            # process as optional/flag argument..
            if ( ! 'name' %in% class(f[[p]]) ) {
                
                # Get default parameter value.
                default <- f[[p]]
                
                # If default value is a 'call', evaluate that 
                # call and take first element of result.
                if ( 'call' %in% class(default) ) {
                    default <- eval(default)[[1]]
                } 
                
                # Verify default value is not NULL.
                # NB: argparser converts NULL to NA anyway. To
                # prevent confusion, use NA for unset values.
                if( is.null(default) ) {
                    stop("NULL default values in pipeline function '", 
                        pipe.func, "'")
                }
                
                # If default value is FALSE, treat as flag..
                if( isFALSE(default) ) {
                    group <- 'flag'
                    default <- NA
                    short <- NA # NB: argparser sets automatically
                } else { # ..otherwise treat as optional argument.
                    group <- 'optional'
                    short <- '-h' # NB: blocks short form
                }
                
                # If argument is itself short form, prepend a 
                # single hyphen, block argparser short form..
                if ( nchar(p) == 1 ) {
                    arg <- paste0('-', p)
                    short <- '-h' # NB: blocks short form
                } else { # ..otherwise prepend two hyphens.
                    arg <- paste0('--', p)
                }
                
            } else { # ..otherwise process as positional argument.
                arg <- p
                group <- 'positional'
                default <- NA
                short <- '-h' # NB: blocks short form
            }
            
            # Set info for this parameter.
            params$arg[p] <- arg
            params$short[p] <- short
            params$group[p] <- group
            params$default[[p]] <- default
        }
        
    } else {
        
        params <- NULL
    }
    
    # Set error message for pipeline docs processing.
    err.msg <- paste0("failed to read pipeline docs for - '", pipeline, "'")
    
    # Compose pipeline Rd file name.
    rd.file <- paste0(pipe.func, '.Rd')
    
    # Load shmootl Rd database.
    rd.db <- tools::Rd_db('shmootl')
    
    # Get Rd content for pipeline function.
    rd.content <- rd.db[[rd.file]]
    
    # Get top-level Rd tags.
    rd.tags <- sapply(rd.content, attr, 'Rd_tag')
    
    # Get title index.
    title.index <- which( rd.tags == '\\title' )
    if ( length(title.index) != 1 ) {
        stop(err.msg)
    }
    
    # Get pipeline function title.
    title <- as.character(rd.content[[title.index]])
    
    # Get description index.
    desc.index <- which( rd.tags == '\\description' )
    
    # Get pipeline function description.
    if ( length(desc.index) == 1 ) {
        description <- paste0(lapply(rd.content[[desc.index]], 
            as.character), collapse='')
    } else if ( length(desc.index) == 0 ) {
        description <- NULL
    } else {
        stop(err.msg)
    }
    
    # Remove description if essentially identical to title.
    if ( stripWhite(description) == title ) {
        description <- NULL
    }
    
    # Get details index.
    detl.index <- which( rd.tags == '\\details' )
    
    # Get pipeline function details.
    if ( length(detl.index) == 1 ) {
        details <- paste0(lapply(rd.content[[detl.index]], 
            as.character), collapse='')
    } else if ( length(detl.index) == 0 ) {
        details <- NULL
    } else {
        stop(err.msg)
    }
    
    # Get arguments index.
    args.index <- which( rd.tags == '\\arguments' )
    
    # Get pipeline function parameter info.
    if ( length(args.index) == 1 ) {
        
        # Get argument docs Rd tags.
        arg.tags <- sapply(rd.content[[args.index]], attr, 'Rd_tag')
        
        # Get indices of tags for individual argument items.
        item.indices <- which( arg.tags == '\\item' )
        
        # Get documented parameters.
        doc.params <- sapply( item.indices, function(i) as.character(
            rd.content[[args.index]][[i]][[1]] ) )
        
        undefined <- setdiff(doc.params, def.params)
        if ( length(undefined) > 0 ) {
            stop("parameters documented but not defined - '", 
                 toString(undefined), "'")
        }
        
        undocumented <- setdiff(def.params, doc.params)
        if ( length(undocumented) > 0 ) {
            stop("parameters defined but not documented - '", 
                 toString(undocumented), "'")
        }
        
        # Get parameter help docs.
        for ( i in 1:length(item.indices) ) {
            item.index <- item.indices[i]
            p <- doc.params[i]
            help.list <- lapply(rd.content[[args.index]][[item.index]][[2]], as.character)
            params$help[p] <- paste0(help.list, collapse='')
        }
        
    } else if ( length(args.index) > 1 ) {
        stop(err.msg)
    }
    
    return( list(command=command, title=title, description=description, 
        details=details, params=params) )
}

# getPipelines -----------------------------------------------------------------
#' Get \pkg{shmootl} pipeline names.
#' 
#' @return Character vector of \pkg{shmootl} pipeline names.
#' 
#' @export
#' @rdname getPipelines
getPipelines <- function() {
    
    # Load shmootl Rd database.
    db <- tools::Rd_db('shmootl')
    
    # Get names of doc files matching pipeline docs pattern.
    rd_files <- names(db)[ grepl(const$pattern$pipe.docs, names(db)) ]
    
    # Extract pipeline names from matching doc file names.
    m <- regexec(const$pattern$pipe.docs, rd_files)
    matches <- regmatches(rd_files, m)
    pipelines <- sapply(matches, getElement, 3)
    
    return(pipelines)
}

# getPipelineUsage -------------------------------------------------------------
#' Get default usage string for \pkg{shmootl} pipelines.
#' 
#' @return Character vector containing the \pkg{shmootl} pipeline default usage 
#' string.
#' 
#' @keywords internal
#' @rdname getPipelineUsage
getPipelineUsage <- function() {
    
    usage <- paste0(
    "usage: Rscript -e 'library(shmootl)' -e 'run()' <pipeline> [-h] ... \n\n",
    "Run the given shmootl pipeline.\n\n",
    "pipelines:\n")

    pipelines <- getPipelines()
    
    pipeline.widths <- nchar(pipelines)
    
    field.width <- max(pipeline.widths) + 6
    
    padding.widths <- field.width - pipeline.widths
    
    paddings <- sapply( padding.widths, function(w) 
        paste( rep(' ', w), collapse='') )
    
    titles <- sapply(pipelines, function(p) getPipelineInfo(p)$title)
    
    pipeline.listing <- paste0("  ", pipelines, paddings, 
        titles, "\n", collapse='')
    
    usage <- paste0(usage, pipeline.listing)
    
    return(usage)
}

# End of pipeline.R ############################################################