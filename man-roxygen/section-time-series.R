#' @section Time-Series Phenotypes:
#' 
#' A set of phenotypes can be designated as a time-series by naming each
#' phenotype with the time point at which phenotype observations were made
#' (e.g. \code{'0.0'}, \code{'1.0'}, \code{'2.0'}). Time points can be in
#' any unit, but must be non-negative, monotonically increasing, and have
#' a consistent time step. If some time points are missing, the resulting
#' gap in time must be a multiple of the time step.
