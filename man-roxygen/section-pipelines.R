#' @section Package Pipelines:
#' 
#' Within \pkg{\link{shmootl}}, package pipelines are fully-documented functions that
#' follow specific conventions regarding function name and formal arguments.
#' Pipeline function names are of the form \code{'run_<pipeline>'} (where
#' \code{'<pipeline>'} is replaced by the pipeline name).
#' 
#' Package pipelines can be run from a command line with \code{Rscript} as follows:
#' 
#' \code{Rscript -e 'library(shmootl)' -e 'run()' <pipeline> [-h] [<args>]}
#' 
#' ...where \code{\link{run}} is the general \pkg{shmootl} pipeline-running
#' function, \code{<pipeline>} is the name of the pipeline to run, and
#' \code{<args>} represents any arguments to be passed to the given pipeline.
#' To see the available options for a pipeline, input the help flag (\code{-h})
#' after the name of the pipeline. To list all available pipelines, input the
#' help flag (\code{-h}) without a pipeline name. (See usage examples below.)
#' 
#' @examples
#' \dontrun{
#' # To list available pipelines from the command line.
#' Rscript -e 'library(shmootl)' -e 'run()' -h
#' 
#' # To run the 'scanone' pipeline from the command line.
#' Rscript -e 'library(shmootl)' -e 'run()' scanone \
#'     --infile input.csv --h5file output.hdf5
#' 
#' # To run the 'scanone' pipeline from within R.
#' run_scanone(infile='input.csv', h5file='output.hdf5')
#' }
