# Start of pipeline.R ##########################################################

# getAnalysisTitle -------------------------------------------------------------
#' Get analysis title for the given \pkg{shmootl} analysis pipeline.
#' 
#' @param pipeline Character vector of \pkg{shmootl} analysis pipeline names.
#' 
#' @return Character vector in which each element contains the title of the
#' \pkg{shmootl} analysis pipeline in the corresponding element of the input
#' vector.
#' 
#' @keywords internal
#' @rdname getAnalysisTitle
getAnalysisTitle <- function(pipeline) {
    
    stopifnot( is.character(pipeline) )
    stopifnot( length(pipeline) > 0 )
    stopifnot( all( pipeline %in% const$supported.analyses) )
    
    indices <- match(pipeline, const$supported.analyses)
    results <- names(const$supported.analyses)[indices]
    
    return(results)
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
    return( get( getPipelineFunctionName(pipeline) ) )
}

# getPipelineFunctionName ------------------------------------------------------
#' Get function name for the given \pkg{shmootl} pipeline.
#' 
#' @param pipeline Character vector of \pkg{shmootl} pipeline names.
#' 
#' @return Character vector in which each element contains the name
#' of the \pkg{shmootl} pipeline function called by the pipeline
#' named in the corresponding element of the input vector.
#' 
#' @keywords internal
#' @rdname getPipelineFunctionName
getPipelineFunctionName <- function(pipeline) {
    stopifnot( is.character(pipeline) )
    stopifnot( length(pipeline) > 0 )
    stopifnot( all( pipeline %in% getPkgPipelineNames() ) )
    return( paste0('run_', pipeline) )
}

# getPipelineName --------------------------------------------------------------
#' Get \pkg{shmootl} pipeline name for the given function.
#' 
#' @param pipe.func.name Character vector of
#' \pkg{shmootl} pipeline function names.
#' 
#' @return Character vector in which each element contains the name
#' of the \pkg{shmootl} pipeline implemented by the function named
#' in the corresponding element of the input vector.
#' 
#' @keywords internal
#' @rdname getPipelineName
getPipelineName <- function(pipe.func.name) {
    
    stopifnot( is.character(pipe.func.name) )
    stopifnot( length(pipe.func.name) > 0 )
    stopifnot( all( pipe.func.name %in% getPkgPipelineFunctionNames() ) )
    
    # NB: valid package pipeline function names are guaranteed to match.
    m <- regexec(const$pattern$pipe.func, pipe.func.name)
    matches <- regmatches(pipe.func.name, m)
    pipeline <- sapply(matches, getElement, 2)
    
    return(pipeline)
}

# getPkgAnalysisPipelineNames --------------------------------------------------
#' Get \pkg{shmootl} analysis pipeline names.
#' 
#' @return Character vector of \pkg{shmootl} analysis pipeline names.
#' 
#' @keywords internal
#' @rdname getPkgAnalysisPipelineNames
getPkgAnalysisPipelineNames <- function() {
    return( getPkgPipelineNames(group='analysis') )
}

# getPkgPipelineDocs -----------------------------------------------------------
#' Get \pkg{shmootl} pipeline docs.
#' 
#' Package pipeline docs are taken from the loaded \pkg{shmootl} package docs,
#' if available; otherwise they are taken from installed \pkg{shmootl} manual
#' pages.
#' 
#' @param pipelines Pipelines for which pipeline docs should be retrieved.
#' If no pipeline is specified, all pipeline docs are returned.
#' 
#' @return Named list of \pkg{shmootl} pipeline docs.
#' 
#' @keywords internal
#' @rdname getPkgPipelineDocs
getPkgPipelineDocs <- function(pipelines=NULL) {
    
    pkg.docs <- NULL
    
    # Try to get shmootl package docs from loaded shmootl package.
    # EG: devtools::load_all('shmootl')
    env <- environment()
    while ( ! identical( env, emptyenv() ) ) {
        
        if ( environmentName(env) == 'package:shmootl' ) {
            
            pkg.path <- attr(env, 'path')
            
            if ( ! is.null(pkg.path) ) {
                
                pkg.dirs <- list.dirs(pkg.path, full.names=FALSE, recursive=FALSE)
                
                if ( 'man' %in% pkg.dirs ) {
                    
                    pkg.docs <- tryCatch({
                        result <- tools::Rd_db(dir=pkg.path)
                    }, error=function(e) {
                        result <- NULL
                    })
                    
                    break
                }
            }
        }
        
        env <- parent.env(env)
    }
    
    # If could not find loaded shmootl package, try to get
    # shmootl package docs from installed shmootl package.
    # EG: library(shmootl)
    if ( length(pkg.docs) == 0 ) { # NULL or character(0)
        pkg.docs <- tryCatch({
            result <-  tools::Rd_db(package='shmootl')
        }, error=function(e) {
            stop("failed to access shmootl package docs")
        })
    }
    
    # Keep only pipeline function documentation.
    pipe.docs <- pkg.docs[ grepl(const$pattern$pipe.docs, names(pkg.docs)) ]
    
    # Rename pipeline docs from pattern-matched documentation names.
    m <- regexec(const$pattern$pipe.docs, names(pipe.docs))
    matches <- regmatches(names(pipe.docs), m)
    names(pipe.docs) <- sapply(matches, getElement, 2)
    
    # If pipelines specified, keep only docs for specified pipelines.
    if ( ! is.null(pipelines) ) {
        
        stopifnot( is.character(pipelines) )
        stopifnot( length(pipelines) > 0 )
        stopif( anyDuplicated(pipelines) )
        
        unknown <- pipelines[ ! pipelines %in% names(pipe.docs) ]
        if ( length(unknown) > 0 ) {
            stop("unknown shmootl pipelines - '", toString(unknown), "'")
        }
        
        pipe.docs <- pipe.docs[pipelines]
    }
    
    return(pipe.docs)
}

# getPkgPipelineFunctionNames --------------------------------------------------
#' Get \pkg{shmootl} pipeline function names.
#' 
#' @param group Pipeline group for which pipeline function names should be
#' retrieved. If no pipeline group is specified, all pipeline function names
#' are returned.
#' 
#' @return Character vector of \pkg{shmootl} pipeline function names.
#' 
#' @export
#' @rdname getPkgPipelineFunctionNames
getPkgPipelineFunctionNames <- function(group=NULL) {
    pipelines <- getPkgPipelineNames(group=group)
    return( paste0('run_', pipelines) )
}

# getPkgPipelineInfo -----------------------------------------------------------
#' Get info for \pkg{shmootl} pipelines.
#' 
#' This function interrogates the function definition and documentation for the 
#' specified pipelines and returns the resulting information as a \code{list}.
#' Each element of this \code{list} is another \code{list} containing
#' information for a specific pipeline, the elements of which are as follows:
#' 
#' \itemize{
#' \item{\code{'command'}} {Shell command used to call pipeline.}
#' \item{\code{'title'}} {Pipeline function title.}
#' \item{\code{'group'}} {Category to which the pipeline belongs.}
#' \item{\code{'description'}} {Pipeline function description. This is not 
#' returned if it duplicates the function title.}
#' \item{\code{'details'}} {Pipeline function details.}
#' \item{\code{'params'}} {A \code{data.frame} containing pipeline function 
#' parameter info. If the function has no parameters, this takes a \code{NULL} 
#' value.}
#' }
#'  
#' Elements \code{'command'}, \code{'title'}, and \code{'group'} are guaranteed
#' to contain their respective values. The remaining elements are only returned
#' if corresponding documentation is available; otherwise, they are assigned a
#' \code{NULL} value.
#' 
#' @param pipelines Names of one or more \pkg{shmootl} pipelines for which info
#' should be retrieved. If none are specified, pipeline info is returned for all
#' pipelines.
#' 
#' @return A \code{list} of \code{list} objects, each containing info for a
#' specific \pkg{shmootl} pipeline.
#' 
#' @keywords internal
#' @rdname getPkgPipelineInfo
getPkgPipelineInfo <- function(pipelines=NULL) {
    
    # TODO: display pipeline references
    
    # Get package docs for given pipelines.
    pipe.docs <- getPkgPipelineDocs(pipelines)
    
    # Replace pipelines argument with explicit pipeline names.
    pipelines <- names(pipe.docs)
    
    # Init pipeline info.
    pinfo <- structure(vector('list', length(pipelines)), names=pipelines)
    
    # Set pipeline info for each given pipeline.
    for ( pipeline in pipelines ) {
        
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
            
            for ( i in seq_along(params) ) {
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
                            
                            param.group <- 'flag'
                            param.type <- 'logical'
                            
                        } else { # ..otherwise treat as optional argument.
                            
                            param.group <- 'optional'
                            
                            if ( is.logical(default) ) {
                                param.type <- 'logical'
                            } else if ( is.integer(default) ) {
                                param.type <- 'integer'
                            } else if ( is.numeric(default) ) {
                                param.type <- 'numeric'
                            } else if ( is.character(default) ) {
                                param.type <- 'character'
                            } else {
                                stop("parameter default value is of unknown type - '",
                                    toString(default), "'")
                            }
                        }
                        
                    # ..otherwise if argument has multiple values,
                    # treat as avector of possible choices..
                    } else if ( length(default) > 1 ) {
                        
                        param.group <- 'choice'
                        
                        if ( is.logical(default) ) {
                            param.type <- 'logical'
                        } else if ( is.integer(default) ) {
                            param.type <- 'integer'
                        } else if ( is.numeric(default) ) {
                            param.type <- 'numeric'
                        } else if ( is.character(default) ) {
                            param.type <- 'character'
                        } else {
                            stop("parameter choice vector is of unknown type - '",
                                toString(default), "'")
                        }
                        
                    # ..otherwise treat as a compound argument to be
                    # loaded from file or directly from the command-line.
                    } else {
                        
                        param.group <- 'compound'
                        
                        if ( is.logical(default) ) {
                            param.type <- 'logical'
                        } else if ( is.integer(default) ) {
                            param.type <- 'integer'
                        } else if ( is.numeric(default) ) {
                            param.type <- 'numeric'
                        } else if ( is.character(default) ) {
                            param.type <- 'character'
                        } else if ( is.mapping(default) ) {
                            param.type <- 'mapping'
                        }  else {
                            stop("parameter default vector is of unknown type - '",
                                toString(default), "'")
                        }
                    }
                    
                    # Verify default value is not NULL.
                    # NB: argparser converts NULL to NA anyway. To
                    # prevent confusion, use NA for unset values.
                    if( is.null(default) ) {
                        stop("parameter '", p, "' has default value NULL in",
                            " pipeline function '", pipe.func, "'")
                    }
                    
                    # If argument is itself short form, prepend one hyphen..
                    if ( nchar(p) == 1 ) {
                        
                        # Check that argument is a flag. NB: argparser does not
                        # display short-form arguments correctly if not a flag.
                        if ( param.group != 'flag' ) {
                            stop("cannot use short-form for ", param.group,
                                " parameter")
                        }
                        
                        arg <- paste0('-', p)
                        
                    } else { # ..otherwise prepend two hyphens.
                        
                        arg <- paste0('--', p)
                    }
                    
                } else { # ..otherwise process as positional argument.
                    
                    arg <- p
                    param.group <- 'positional'
                    default <- NA_character_
                    param.type <- 'character'
                    
                }
                
                # Set info for this parameter.
                params$arg[p] <- arg
                params$group[p] <- param.group
                params$default[[p]] <- default
                params$type[p] <- param.type
            }
            
        } else {
            
            params <- NULL
        }
        
        # Set error message for pipeline docs processing.
        err.msg <- paste0("failed to process docs for pipeline - '", pipeline, "'")
        
        # Get Rd manual for pipeline function.
        pipe.manual <- pipe.docs[[pipeline]]
        
        # Get top-level Rd tags.
        rd.tags <- sapply(pipe.manual, attr, 'Rd_tag')
        
        # Get title index.
        title.index <- which( rd.tags == '\\title' )
        if ( length(title.index) != 1 ) {
            stop(err.msg)
        }
        
        # Get pipeline function title.
        title <- paste0(lapply(pipe.manual[[title.index]],
            as.character), collapse='')
        
        # Get description index.
        desc.index <- which( rd.tags == '\\description' )
        
        # Get pipeline function description.
        if ( length(desc.index) == 1 ) {
            description <- paste0(lapply(pipe.manual[[desc.index]],
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
            details <- paste0(lapply(pipe.manual[[detl.index]],
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
            arg.tags <- sapply(pipe.manual[[args.index]], attr, 'Rd_tag')
            
            # Get indices of tags for individual argument items.
            item.indices <- which( arg.tags == '\\item' )
            
            # Get documented parameters.
            doc.params <- sapply( item.indices, function(i) as.character(
                pipe.manual[[args.index]][[i]][[1]] ) )
            
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
            for ( i in seq_along(item.indices) ) {
                
                # Get parameter item index.
                item.index <- item.indices[i]
                
                # Get parameter name.
                p <- doc.params[i]
                
                # Get concatenated list of help strings for this parameter.
                help.list <- lapply(pipe.manual[[args.index]][[item.index]][[2]],
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
        
        # Get concept tags.
        concept.index <- which( rd.tags == '\\concept' )
        if ( length(concept.index) == 1 ) {
            
            # Get concept lines from pipeline Rd manual.
            concept.lines <- lapply( pipe.manual[[concept.index]], function(line)
                gsub('(?:^[[:space:]]+)|(?:[[:space:]]+$)', '', as.character(line) ) )
            
            # Filter empty concept lines.
            concept.lines <- concept.lines[ sapply(concept.lines, nzchar) ]
            
            # Get concept tags from concept lines.
            concepts <- unlist( lapply(concept.lines, strsplit, '[[:space:]]+') )
            
        } else if ( length(concept.index) == 0 ) {
            concepts <- character()
        } else {
            stop(err.msg)
        }
        
        # Get any pipeline groups in concept tags.
        # NB: assumes one capture group in pattern.
        pattern.match <- regexpr(const$pattern$pipe.group, concepts, perl=TRUE)
        indices <- which(pattern.match != -1)
        capture.first <- as.integer( attr(pattern.match, 'capture.start') )[indices]
        capture.length <- as.integer( attr(pattern.match, 'capture.length') )[indices]
        capture.last <- capture.first + capture.length - 1
        pipe.groups <- substring(concepts[indices], capture.first, capture.last)
        
        # Assign pipeline to group.
        if ( length(pipe.groups) == 1 ) {
            pipe.group <- pipe.groups[1]
        } else if ( length(pipe.groups) == 0 ) {
            pipe.group <- const$default$pipeline.group
        } else {
            stop("pipeline '", pipeline, "' assigned to multiple groups - '",
                toString(pipe.groups), "'")
        }
        
        pinfo[[pipeline]] <- list(command=command, title=title, group=pipe.group,
            description=description, details=details, params=params)
    }
    
    return(pinfo)
}

# getPkgPipelineNames ----------------------------------------------------------
#' Get \pkg{shmootl} pipeline names.
#' 
#' Package pipeline names are taken from the loaded \pkg{shmootl} namespace, if
#' available; otherwise they are taken from installed \pkg{shmootl} package.
#' 
#' @param group Pipeline group for which pipeline names should be retrieved.
#' If no pipeline group is specified, all pipeline names are returned.
#' 
#' @return Character vector of \pkg{shmootl} pipeline names.
#' 
#' @export
#' @rdname getPkgPipelineNames
getPkgPipelineNames <- function(group=NULL) {
    
    pipe.func.names <- NULL
    
    # Try to get pipeline function names from loaded shmootl package.
    # NB: devtools::load_all('shmootl')
    env <- environment()
    while ( ! identical( env, emptyenv() ) ) {
        
        if ( environmentName(env) == 'package:shmootl' ) {
            pipe.func.names <- ls(env, pattern=const$pattern$pipe.func)
            break
        }
        
        env <- parent.env(env)
    }
    
    # If could not find loaded shmootl package, try to get
    # pipeline function names from installed shmootl package.
    # NB: library(shmootl)
    if ( is.null(pipe.func.names) ) {
        pipe.func.names <- pkg.docs <- tryCatch({
            result <- ls('package:shmootl', pattern=const$pattern$pipe.func)
        }, error=function(e) {
            stop("failed to access shmootl package namespace")
        })
    }
    
    # Extract pipeline names from matching function names.
    m <- regexec(const$pattern$pipe.func, pipe.func.names)
    matches <- regmatches(pipe.func.names, m)
    pipelines <- sapply(matches, getElement, 2)
    
    if ( ! is.null(group) ) {
        
        stopifnot( isSingleString(group) )
        stopifnot( group %in% const$pipe.groups )
        
        pinfo <- getPkgPipelineInfo(pipelines)
        mask <- sapply(pipelines, function(p) pinfo[[p]]$group == group)
        pipelines <- pipelines[mask]
    }
    
    return(pipelines)
}

# getPkgPipelineUsage ----------------------------------------------------------
#' Get default usage string for \pkg{shmootl} pipelines.
#' 
#' @return Character vector containing the \pkg{shmootl} pipeline default usage 
#' string.
#' 
#' @keywords internal
#' @rdname getPkgPipelineUsage
getPkgPipelineUsage <- function() {
    
    usage <- paste0(
        "usage: Rscript -e 'library(shmootl)' -e 'run()' <pipeline> [-h] ... \n\n",
        "Run the given shmootl pipeline.\n")
    
    # Get pipeline info for each pipeline.
    pinfo <- getPkgPipelineInfo()
    
    # Get pipeline names from pipeline info.
    pipelines <- names(pinfo)
    
    # Get pipeline padding strings.
    pipeline.widths <- nchar(pipelines)
    field.width <- max(pipeline.widths) + 6
    padding.widths <- field.width - pipeline.widths
    paddings <- sapply( padding.widths, function(w) 
        paste( rep(' ', w), collapse='') )
    
    # Get title of each pipeline.
    pipeline.titles <- unlist( lapply(pipelines, function(p) pinfo[[p]]$title) )
    
    # Get group name of each pipeline.
    pipe.groups <- unlist( lapply(pipelines, function(p) pinfo[[p]]$group) )
    
    # Get non-redundant vector of pipeline group names.
    groups <- const$pipe.groups[ const$pipe.groups %in% pipe.groups ]
    
    # Create grouped pipeline listings.
    pipeline.listings <- structure(character( length(groups) ), names=groups)
    for ( group in groups ) {
        indices <- which( pipe.groups == group )
        group.head <- paste0("\n", group, ":\n")
        group.body <- paste0("  ", pipelines[indices], paddings[indices],
            pipeline.titles[indices], "\n", collapse='')
        pipeline.listings[group] <- paste0(group.head, group.body)
    }
    
    return( paste0( c(usage, pipeline.listings) ) )
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
    
    stopifnot( isSingleString(pipeline) )
    
    # Get pipeline info.
    pinfo <- getPkgPipelineInfo(pipeline)
    
    # Compose argument parser description by combining
    # the function title, description, and details.
    description <- paste0(pinfo[[pipeline]]$title, '\n')
    for ( content in c(pinfo[[pipeline]]$description, pinfo[[pipeline]]$details) ) {
        if ( ! is.null(content) ) {
            description <- paste0(description, content)
        }
    }
    
    # Create argument parser.
    ap <- argparser::arg_parser(name=pinfo[[pipeline]]$command,
        description=description)
    
    # Add any parameters to the argument parser.
    if ( ! is.null(pinfo[[pipeline]]$params) ) {
        
        for ( i in seq_along(pinfo[[pipeline]]$params$arg) ) {
            
            ap <- argparser::add_argument(ap,
                arg   = pinfo[[pipeline]]$params$arg[i],
                short = '-h', # NB: blocks short form
                flag  = (pinfo[[pipeline]]$params$group[i] == 'flag'),
                help  = pinfo[[pipeline]]$params$help[i]
            )
            
            # If compound parameter, add additional argument to load from file.
            if ( pinfo[[pipeline]]$params$group[i] == 'compound' ) {
                
                file.param <- paste0(pinfo[[pipeline]]$params$arg[i], '.file')
                stopif( file.param %in% pinfo[[pipeline]]$params$arg )
                
                ap <- argparser::add_argument(ap,
                    arg   = file.param,
                    short = '-h', # NB: blocks short form
                    flag  = FALSE,
                    help  = pinfo[[pipeline]]$params$help[i]
                )
            }
        }
    }
    
    # Add pipeline info to argument parser.
    stopif( 'pinfo' %in% names( attributes(ap) ) )
    attr(ap, 'pinfo') <- pinfo[[pipeline]]
    
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
        param.group <- pinfo$params$group[p]
        param.type <- pinfo$params$type[p]
        arg <- args[[p]]
        
        # If parameter is compound, get name and
        # value of corresponding file argument.
        if ( param.group == 'compound' ) {
            fp <- paste0(p, '.file')
            farg <- args[[fp]]
        }
        
        # If argument has been specified..
        if ( ! is.na(arg) ) {
            
            # ..and a parameter type has been specified..
            if ( ! is.na(param.type) ) {
                
                # ..set argument according to its group.
                if ( param.group %in% c('flag', 'optional', 'positional') ) {
                    
                    if ( param.type == 'integer' ) {
                        arg <- as.numeric(arg)
                        stopifnot( isWholeNumber(arg) )
                        arg <- as.integer(arg)
                    } else if ( param.type == 'numeric' ) {
                        arg <- as.numeric(arg)
                    } else if ( param.type == 'logical' ) {
                        arg <- as.logical(arg)
                    } else if ( param.type != 'character' ) {
                        stop("unknown parameter type - '", param.type, "'")
                    }
                    
                } else if ( param.group == 'choice' ) {
                    
                    if ( param.type == 'integer' ) {
                        arg <- as.numeric(arg)
                        stopifnot( isWholeNumber(arg) )
                        arg <- as.integer(arg)
                    } else if ( param.type == 'numeric' ) {
                        arg <- as.numeric(arg)
                    } else if ( param.type == 'logical' ) {
                        arg <- as.logical(arg)
                    } else if ( param.type != 'character' ) {
                        stop("unknown parameter type - '", param.type, "'")
                    }
                    
                    if ( ! arg %in% default ) {
                        stop("unsupported choice for parameter '", p, "' - '", arg, "'")
                    }
                    
                } else if ( param.group == 'compound' ) {
                    
                    if ( ! is.na(farg) ) {
                        stop("incompatible arguments - '", toString( c(p, fp) ), "'")
                    }
                    
                    if ( param.type == 'mapping' ) {
                        arg <- loadMapping(line=arg)
                    } else {
                        arg <- loadVector(line=arg, type=param.type)
                    }
                }
            }
            
        } else { # ..otherwise if regular argument not specified, set as appropriate.
            
            # NB: argparser will catch missing positional arguments.
            
            if ( param.group %in% c('flag', 'optional') ) {
                
                arg <- default
            
            } else if ( param.group == 'choice' ) {
                
                arg <- default[[1]]
            
            } else if ( param.group == 'compound' ) {
                
                # If file argument specified, set value from file..
                if ( ! is.na(farg) ) {
                    
                    if ( param.type == 'mapping' ) {
                        arg <- loadMapping(file=farg)
                    } else {
                        arg <- loadVector(file=farg, type=param.type)
                    }
                    
                } else { # ..otherwise set default value.
                    
                    arg <- default
                }
            }
        }
        
        # If parameter is compound, remove file
        # argument, as no longer needed.
        if ( param.group == 'compound' ) {
            args[[fp]] <- NULL
        }
        
        args[[p]] <- arg
    }
    
    return(args)
}

# resolveAnalysisTitle ---------------------------------------------------------
#' Resolve analysis title.
#' 
#' @param analyses Vector of \pkg{shmootl} analysis pipeline names or titles.
#' If no analysis pipelines are specified, all supported analysis titles are
#' returned.
#' 
#' @return If the \code{analyses} parameter is specified, this is a character
#' vector in which each element contains the title of the \pkg{shmootl} analysis
#' pipeline referenced in the corresponding element of the input vector. If the
#' \code{analyses} parameter is \code{NULL}, this is a character vector of all
#' supported analysis titles.
#' 
#' @keywords internal
#' @rdname resolveAnalysisTitle
resolveAnalysisTitle <- function(analyses=NULL) {
    
    if ( ! is.null(analyses) ) {
        
        stopifnot( is.character(analyses) )
        stopifnot( length(analyses) > 0 )
        
        # Get supported analysis names and titles.
        supported.names <- unname(const$supported.analyses)
        supported.titles <- names(const$supported.analyses)
        
        # Get indices of supported analyses matched by each specified analysis.
        index.list <- lapply( analyses, function(analysis)
            which( supported.names == analysis | supported.titles == analysis ) )
        
        # Get number of supported analyses matched by each specified analysis.
        match.counts <- lengths(index.list)
        
        unresolved <- analyses[ match.counts == 0 ]
        if ( length(unresolved) > 0 ) {
            stop("analyses could not be resolved - '", toString(unresolved), "'")
        }
        
        ambiguous <- analyses[ match.counts > 1 ]
        if ( length(ambiguous) > 0 ) {
            stop("analyses could not be unambiguously resolved - '",
                toString(ambiguous), "'")
        }
        
        res <- names(const$supported.analyses)[ unlist(index.list) ]
        
    } else {
        
        res <- names(const$supported.analyses)
    }
    
    return(res)
}

# End of pipeline.R ############################################################