# Start of pipeline.R ##########################################################

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
    command <- paste("Rscript -e 'library(shmootl)' -e 'run()'", pipeline)
    
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
            group = factor(character(n.params), levels=const$param.groups),
            default = vector('list', n.params),
            type = character(n.params),
            help = character(n.params)
        )
        
        for ( i in getIndices(params) ) {
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
                    default <- eval(default)
                }
                
                # If default is of length one, treat
                # as flag or simple optional argument..
                if ( length(default) == 1 ) {
                    
                    # If default value is FALSE, treat as flag..
                    if( isFALSE(default) ) {
                        
                        group <- 'flag'
                        ptype <- 'logical'
                        
                    } else { # ..otherwise treat as optional argument.
                        
                        group <- 'optional'
                        
                        if ( is.logical(default) ) {
                            ptype <- 'logical'
                        } else if ( is.integer(default) ) {
                            ptype <- 'integer'
                        } else if ( is.numeric(default) ) {
                            ptype <- 'numeric'
                        } else if ( is.character(default) ) {
                            ptype <- 'character'
                        } else {
                            stop("parameter default value is of unknown type - '", default, "'")
                        }
                    }
                
                # ..otherwise if argument has multiple values,
                # treat as avector of possible choices..
                } else if ( length(default) > 1 ) {
                    
                    group <- 'choice'
                    
                    if ( is.logical(default) ) {
                        ptype <- 'logical'
                    } else if ( is.integer(default) ) {
                        ptype <- 'integer'
                    } else if ( is.numeric(default) ) {
                        ptype <- 'numeric'
                    } else if ( is.character(default) ) {
                        ptype <- 'character'
                    } else {
                        stop("parameter choice vector is of unknown type - '", default, "'")
                    }
                    
                # ..otherwise treat as a compound argument to be
                # loaded from file or directly from the command-line.
                } else {
                    
                    group <- 'compound'
                    
                    if ( is.logical(default) ) {
                        ptype <- 'logical'
                    } else if ( is.integer(default) ) {
                        ptype <- 'integer'
                    } else if ( is.numeric(default) ) {
                        ptype <- 'numeric'
                    } else if ( is.character(default) ) {
                        ptype <- 'character'
                    } else if ( is.mapping(default) ) {
                        ptype <- 'mapping'
                    }  else {
                        stop("parameter default vector is of unknown type - '", default, "'")
                    }
                }
                
                # Verify default value is not NULL.
                # NB: argparser converts NULL to NA anyway. To
                # prevent confusion, use NA for unset values.
                if( is.null(default) ) {
                    stop("parameter '", p, "' has default value NULL in pipeline function '",
                        pipe.func, "'")
                }
                
                # If argument is itself short form, prepend one hyphen..
                if ( nchar(p) == 1 ) {
                    
                    # Check that argument is a flag. NB: argparser does not
                    # display short-form arguments correctly if not a flag.
                    if ( group != 'flag' ) {
                        stop("cannot use short-form for ", group, " parameter")
                    }
                    
                    arg <- paste0('-', p)
                    
                } else { # ..otherwise prepend two hyphens.
                    
                    arg <- paste0('--', p)
                }
                
            } else { # ..otherwise process as positional argument.
                
                arg <- p
                group <- 'positional'
                default <- NA_character_
                ptype <- 'character'
                
            }
            
            # Set info for this parameter.
            params$arg[p] <- arg
            params$group[p] <- group
            params$default[[p]] <- default
            params$type[p] <- ptype
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
    title <- paste0(lapply(rd.content[[title.index]],
        as.character), collapse='')
    
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
        
        # Set parameter help doc string.
        for ( i in getIndices(item.indices) ) {
            
            # Get parameter item index.
            item.index <- item.indices[i]
            
            # Get parameter name.
            p <- doc.params[i]
            
            # Get concatenated list of help strings for this parameter.
            help.list <- lapply(rd.content[[args.index]][[item.index]][[2]],
                as.character)
            params$help[p] <- paste0(help.list, collapse='')
            
            # If parameter is multiple-choice, add list of choices to help..
            if ( params$group[p] == 'choice' ) {
                
                choice.info <- paste0('{', paste0(params$default[[p]], collapse=','), '}')
                params$help[p] <- paste(params$help[p], choice.info)
                
            } else { # ..otherwise if a single valid default value, add to help.
                
                default <- params$default[[p]]
                if ( length(default) == 1 && ! is.na(default) ) {
                    default.info <- paste0('[default: ', default, ']')
                    params$help[p] <- paste(params$help[p], default.info)
                }
            }
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

# prepPipelineArgparser --------------------------------------------------------
#' Prep argument parser for \pkg{shmootl} pipeline.
#' 
#' @param pipeline Name of a \pkg{shmootl} pipeline.
#' 
#' @return Argument parser for the specified \pkg{shmootl} pipeline.
#' 
#' @keywords internal
#' @rdname prepPipelineArgparser
prepPipelineArgparser <- function(pipeline) {
    
    # Get pipeline info.
    pinfo <- getPipelineInfo(pipeline)
    
    # Compose argument parser description by combining
    # the function title, description, and details.
    description <- paste0(pinfo$title, '\n')
    for ( content in c(pinfo$description, pinfo$details) ) {
        if ( ! is.null(content) ) {
            description <- paste0(description, content)
        }
    }
    
    # Create argument parser.
    ap <- argparser::arg_parser(name=pinfo$command, description=description)
    
    # Add any parameters to the argument parser.
    if ( ! is.null(pinfo$params) ) {
        
        for ( i in getIndices(pinfo$params$arg) ) {
            
            ap <- argparser::add_argument(ap,
                arg   = pinfo$params$arg[i],
                short = '-h', # NB: blocks short form
                flag  = (pinfo$params$group[i] == 'flag'),
                help  = pinfo$params$help[i]
            )
            
            # If compound parameter, add additional argument to load from file.
            if ( pinfo$params$group[i] == 'compound' ) {
                
                file.param <- paste0(pinfo$params$arg[i], '.file')
                stopif( file.param %in% pinfo$params$arg )
                
                ap <- argparser::add_argument(ap,
                    arg   = file.param,
                    short = '-h', # NB: blocks short form
                    flag  = FALSE,
                    help  = pinfo$params$help[i]
                )
            }
        }
    }
    
    # Add pipeline info to argument parser.
    stopif( 'pinfo' %in% names( attributes(ap) ) )
    attr(ap, 'pinfo') <- pinfo
    
    return(ap)
}

# procPipelineArgs -------------------------------------------------------------
#' Process arguments parsed by \pkg{argparser}.
#' 
#' @param ap An \pkg{argparser} parser object.
#' @param args Arguments to be processed, having been parsed by \pkg{argparser}.
#' 
#' @return Processed pipeline arguments.
#' 
#' @include const.R
#' @keywords internal
#' @rdname procPipelineArgs
procPipelineArgs <- function(ap, args) {
    
    # Get pipeline info from argument parser.
    stopifnot( 'pinfo' %in% names( attributes(ap) ) )
    pinfo <- attr(ap, 'pinfo')
    
    # Remove argparser special arguments.
    args <- args[ ! names(args) %in% const$special.params ]
    
    # Process argument(s) corresponding to each parameter.
    for ( p in names(pinfo$params$arg) ) {
        
        # Get parameter info.
        default <- pinfo$params$default[[p]]
        group <- pinfo$params$group[p]
        ptype <- pinfo$params$type[p]
        arg <- args[[p]]
        
        # If parameter is compound, get name and
        # value of corresponding file argument.
        if ( group == 'compound' ) {
            fp <- paste0(p, '.file')
            farg <- args[[fp]]
        }
        
        # If argument has been specified..
        if ( ! is.na(arg) ) {
            
            # ..and a parameter type has been specified..
            if ( ! is.na(ptype) ) {
                
                # ..set argument according to its group.
                if ( group %in% c('flag', 'optional', 'positional') ) {
                    
                    if ( ptype == 'integer' ) {
                        arg <- as.numeric(arg)
                        stopifnot( isWholeNumber(arg) )
                        arg <- as.integer(arg)
                    } else if ( ptype == 'numeric' ) {
                        arg <- as.numeric(arg)
                    } else if ( ptype == 'logical' ) {
                        arg <- as.logical(arg)
                    } else if ( ptype != 'character' ) {
                        stop("unknown parameter type - '", ptype, "'")
                    }
                    
                } else if ( group == 'choice' ) {
                    
                    if ( ptype == 'integer' ) {
                        arg <- as.numeric(arg)
                        stopifnot( isWholeNumber(arg) )
                        arg <- as.integer(arg)
                    } else if ( ptype == 'numeric' ) {
                        arg <- as.numeric(arg)
                    } else if ( ptype == 'logical' ) {
                        arg <- as.logical(arg)
                    } else if ( ptype != 'character' ) {
                        stop("unknown parameter type - '", ptype, "'")
                    }
                    
                    if ( ! arg %in% default ) {
                        stop("unsupported choice for parameter '", p, "' - '", arg, "'")
                    }
                    
                } else if ( group == 'compound' ) {
                    
                    if ( ! is.na(farg) ) {
                        stop("incompatible arguments - '", toString( c(p, fp) ), "'")
                    }
                    
                    if ( ptype == 'mapping' ) {
                        arg <- loadMapping(line=arg)
                    } else {
                        arg <- loadVector(line=arg, type=ptype)
                    }
                }
            }
            
        } else { # ..otherwise if regular argument not specified, set as appropriate.
            
            # NB: argparser will catch missing positional arguments.
            
            if ( group %in% c('flag', 'optional') ) {
                
                arg <- default
            
            } else if ( group == 'choice' ) {
                
                arg <- default[[1]]
            
            } else if ( group == 'compound' ) {
                
                # If file argument specified, set value from file..
                if ( ! is.na(farg) ) {
                    
                    if ( ptype == 'mapping' ) {
                        arg <- loadMapping(file=farg)
                    } else {
                        arg <- loadVector(file=farg, type=ptype)
                    }
                    
                    args[[fp]] <- NULL
                    
                } else { # ..otherwise set default value.
                    
                    arg <- default
                }
            }
        }
        
        args[[p]] <- arg
    }
    
    return(args)
}

# End of pipeline.R ############################################################