# Start of qtlfeatures.R #######################################################

# `[.qtlfeatures` --------------------------------------------------------------
#' @export
#' @keywords internal
`[.qtlfeatures` <- function(x, i) {
    others <- otherattributes(x)
    x <- unclass(x)
    x <- x[i]
    class(x) <- c('qtlfeatures', 'list')
    otherattributes(x) <- others
    return(x)
}

# End of qtlfeatures.R #########################################################