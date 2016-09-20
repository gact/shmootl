#' @section Tetradic Samples:
#' 
#' Sample IDs can be used to indicate tetrad membership, even in a \code{cross}
#' object with some missing samples. In a tetradic dataset, sample IDs with a
#' numeric suffix (e.g. \code{'FS101'}) are taken as segregant numbers and used
#' to infer the tetrad to which each sample ID belongs, assuming that tetrads
#' are labelled sequentially, with four samples per tetrad. Sample IDs can also
#' have an alphanumeric suffix (e.g. \code{'FS01A'}), where the numeric part is
#' a tetrad number and the final letter (i.e. \code{'A'}, \code{'B'},
#' \code{'C'}, or \code{'D'}) identifies the individual tetrad member.
