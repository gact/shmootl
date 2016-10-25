# Start of scanonebins.R #######################################################

# TODO: plot.scanonebins

# `[.scanonebins` --------------------------------------------------------------
#' @export
#' @keywords internal
`[.scanonebins` <- function(x, i, j, k, drop=NULL) {
    stopifnot( is.null(drop) )
    others <- otherattributes(x)
    x <- unclass(x)
    x <- x[i, j, k, drop=FALSE]
    class(x) <- c('scanonebins', 'array')
    otherattributes(x) <- others
    return(x)
}

# makeLODBinLabels -------------------------------------------------------------
#' Make LOD bin labels from bin starts.
#'
#' @param bin.starts Numeric vector of bin starts.
#' 
#' @return Character vector of LOD bin labels corresponding to the input bins.
#' 
#' @keywords internal
#' @rdname makeLODBinLabels
makeLODBinLabels <- function(bin.starts) {
    
    if ( length(bin.starts) > 0 ) {
        
        stopifnot( is.numeric(bin.starts) )
        
        if ( bin.starts[1] != const$lod.bin$min.start ) {
            stop("invalid LOD bins start ", bin.starts[1])
        }
        
        if ( length(bin.starts) > 1 ) {
            
            bin.sizes <- diff(bin.starts)
            
            neg.bin.sizes <- bin.sizes[ bin.sizes < 0 ]
            if ( length(neg.bin.sizes) > 0 ) {
                stop("negative LOD bin sizes - '", toString(neg.bin.sizes), "'")
            }
            
            bin.delta <- abs(bin.sizes - const$lod.bin$size)
            err.bin.sizes <- bin.sizes[ bin.delta > .Machine$double.eps^0.5 ]
            if ( length(err.bin.sizes) > 0 ) {
                stop("invalid LOD bin sizes - '", toString(err.bin.sizes), "'")
            }
        }
        
        bin.ends <- bin.starts + const$lod.bin$size
        
        # Get number of digits needed to display bin size number.
        fractional.part <- const$lod.bin$size %% 1.0
        if (fractional.part) {
            digits <- abs( floor( log10(fractional.part) ) )
        } else {
            digits <- 0
        }
        
        # Format bin start and end strings.
        bin.starts.str <- sprintf(paste0('%.', digits, 'f'), bin.starts)
        bin.ends.str <- sprintf(paste0('%.', digits, 'f'), bin.ends)
        
        bin.labels <- paste0('LOD[', bin.starts.str, ',', bin.ends.str, ')')
        
    } else {
        
        bin.labels <- character()
    }
    
    return(bin.labels)
}

# parseLODBinLabels ------------------------------------------------------------
#' Parse bin starts from LOD bin labels.
#'
#' @param bin.labels Character vector of LOD bin labels to a set of LOD bins.
#' 
#' @return Numeric vector of bin starts.
#' 
#' @keywords internal
#' @rdname parseLODBinLabels
parseLODBinLabels <- function(bin.labels) {
    
    if ( length(bin.labels) > 0 ) {
        
        stopifnot( is.character(bin.labels) )
        
        lod.bin.regex <- '^LOD\\[([0-9]+(?:[.][0-9]+)?),([0-9]+(?:[.][0-9]+)?)\\)$'
        
        m <- regexec(lod.bin.regex, bin.labels)
        regmatch.list <- regmatches(bin.labels, m)
        
        invalid.labels <- bin.labels[ lengths(regmatch.list) == 0 ]
        if ( length(invalid.labels) > 0 ) {
            stop("invalid LOD bin labels - '", toString(invalid.labels), "'")
        }
        
        bin.starts <- as.numeric( sapply(regmatch.list, function(x) x[2]) )
        bin.ends <- as.numeric( sapply(regmatch.list, function(x) x[3]) )
        
        if ( bin.starts[1] != const$lod.bin$min.start ) {
            stop("invalid LOD bins start ", bin.starts[1])
        }
        
        bin.sizes <- bin.ends - bin.starts
        
        neg.bin.sizes <- bin.sizes[ bin.sizes < 0 ]
        if ( length(neg.bin.sizes) > 0 ) {
            stop("negative LOD bin sizes - '", toString(neg.bin.sizes), "'")
        }
        
        bin.delta <- abs(bin.sizes - const$lod.bin$size)
        err.bin.sizes <- bin.sizes[ bin.delta > .Machine$double.eps^0.5 ]
        if ( length(err.bin.sizes) > 0 ) {
            stop("invalid LOD bin sizes - '", toString(err.bin.sizes), "'")
        }
    
    } else {
        
        bin.starts <- numeric()
    }
    
    return(bin.starts)
}

# mergeLODBinLabels ------------------------------------------------------------
#' Merge multiple vectors of LOD bin labels.
#'
#' @param ... One or more character vectors containing LOD bin labels.
#' 
#' @return Merged bin label vector.
#' 
#' @keywords internal
#' @rdname mergeLODBinLabels
mergeLODBinLabels <- function(...) {
    
    bin.label.list <- list(...)
    stopifnot( is.list(bin.label.list) )
    stopifnot( length(bin.label.list) > 0 )
    
    # Remove empty bin label vectors.
    bin.label.list <- bin.label.list[ ! unlist( lapply(bin.label.list, is.null) ) ]
    
    if ( length(bin.label.list) > 0 ) {
    
        stopifnot( all( sapply(bin.label.list, is.character) ) )
        
        max.idx <- which.max( lengths(bin.label.list) )
        
        merged.labels <- bin.label.list[[max.idx]]
        
        parseLODBinLabels(merged.labels) # NB: implicitly validates merged labels
        
        for ( bin.labels in bin.label.list ) {
            
            l <- length(bin.labels)
            
            if ( l > 0 && any( bin.labels != merged.labels[1:l] ) ) {
                stop("LOD bin label mismatch")
            }
        }
    
    } else {
        
        merged.labels <- character()
    }
    
    return(merged.labels)
}

# padBins ----------------------------------------------------------------------
#' Pad \code{scanonebins} object to given number of bins.
#'
#' @param x A \code{scanonebins} object.
#' @param num.bins Number of bins to which \code{scanonebins} object
#' should be padded.
#' 
#' @return Input \code{scanonebins} object, padded with zeroes to the
#' specified size.
#' 
#' @importFrom abind abind
#' @keywords internal
#' @rdname padBins
padBins <- function(x, num.bins) {
    
    stopifnot( 'scanonebins' %in% class(x) )
    
    others <- otherattributes(x)
    
    bin.starts <- seq(from=const$lod.bin$min.start,
        by=const$lod.bin$size, length.out=num.bins)
    
    bin.labels <- makeLODBinLabels(bin.starts)
    
    x.labels <- dimnames(x)[[2]]
    
    x.dim <- dim(x)
    
    n <- x.dim[2]
    
    bin.labels <- do.call(mergeLODBinLabels, list(x.labels, bin.labels))
    
    if ( n < num.bins ) {
        
        padding.indices <- (n + 1):num.bins
        padding.labels <- bin.labels[padding.indices]
        
        padding.dimnames <- dimnames(x)
        padding.dimnames[[2]] <- padding.labels
        
        padding.dim <- x.dim
        padding.dim[2] <- length(padding.labels)
        
        padding <- array(0, dim=padding.dim,
            dimnames=padding.dimnames)
        
        x <- abind::abind(x, padding, along=2)
        
        class(x) <- c('scanonebins', 'array')
    }
    
    otherattributes(x) <- others
    
    return(x)
}

# print.summary.scanonebins ----------------------------------------------------
#' Print \code{summary.scanonebins} object.
#'
#' @param x A \code{summary.scanonebins} object.
#' @param ... Unused arguments.
#' 
#' @export print.summary.scanonebins
#' @keywords internal
#' @method print summary.scanonebins
#' @rdname print.summary.scanonebins
print.summary.scanonebins <- function(x, ...) {
    num.perms <- attr(x, 'n.perm')
    cat("FDR LOD thresholds (", num.perms, " permutations)\n", sep='')
    print( matrix(x, dimnames=dimnames(x)) )
}

# summary.scanonebins ----------------------------------------------------------
#' Estimate FDR LOD thresholds from a \code{scanonebins} object.
#' 
#' This function is based on the \pkg{R/qtl} function \code{mqmscanfdr} (see
#' reference below), but uses a precalculated \code{scanonebins} object
#' instead of running permutations directly.
#' 
#' Thresholds are set from the bin edges of the \code{scanonebins} object, and
#' empirical false-discovery rates (FDRs) are calculated at these thresholds.
#' For each false-discovery rate specified in the \code{fdr} parameter, the
#' threshold with the next lowest FDR estimate is returned where possible, or
#' next highest FDR estimate otherwise. Summary results take \code{NA} values
#' in cases where no sensible FDR can be returned (e.g. fewer loci above a
#' threshold in the main scanone result than in the corresponding permutation
#' results).
#' 
#' @param object A \code{scanonebins} object.
#' @param scanone.result A \code{scanone} object created from the same data that
#' was used to create the \code{scanonebins} object.
#' @template param-lodcolumns
#' @param fdr False discovery rates for which an approximate LOD thresholds
#' should be estimated.
#' @param ... Unused arguments.
#' 
#' @return A \code{summary.scanonebins} matrix containing LOD thresholds that
#' correspond approximately to the specified false discovery rates.
#' 
#' @references Arends D, Prins P, Jansen RC, Broman KW (2010) R/qtl:
#' high-throughput multiple QTL mapping. Bioinformatics 26(23):2990-2.
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/20966004}{PubMed})
#' 
#' @export
#' @keywords internal
#' @method summary scanonebins
#' @rdname summary.scanonebins
summary.scanonebins <- function(object, scanone.result, lodcolumns=NULL, fdr=0.01, ...) {
    
    stopifnot( 'scanone' %in% class(scanone.result) )
    stopifnot( all( ! is.na(fdr) & is.finite(fdr) & fdr > 0 & fdr < 1 ) )
    stopifnot( emptyArgs(...) )
    
    # Get sorted unique FDR values. 
    fdr <- unique( sort(fdr) )
    
    # Check that LOD column names match between scanone and scanonebins objects.
    scan.lodcols <- getLodColIndices(scanone.result)
    scan.lodcol.names <- colnames(scanone.result)[scan.lodcols]
    perm.lodcol.names <- dimnames(object)[[3]]
    if ( length(perm.lodcol.names) != length(scan.lodcol.names) ||
        any( perm.lodcol.names != scan.lodcol.names ) ) {
        stop("LOD column mismatch")
    }
    
    # Get specified LOD column indices.
    lodcol.indices <- getLodColIndices(scanone.result,
        lodcolumns=lodcolumns, strict=TRUE)
    stopifnot( length(lodcol.indices) > 0 )
    
    # If no error occurred above, there is a single LOD column.
    if ( is.null(lodcolumns) ) {
        lodcolumns <- scan.lodcols - 2 # 'chr' and 'pos' columns
    }
    
    # Subset scanone permutations for LOD columns of interest.
    object <- object[,, lodcolumns]
    
    # Get number of permutations.
    num.perms <- dim(object)[1]
    
    # Bin LOD values of scanone result.
    scan.bins <- binLODValues(scanone.result)
    
    # Get threshold values from merged bin labels.
    scan.bin.labels <- dimnames(scan.bins)[[2]]
    perm.bin.labels <- dimnames(object)[[2]]

    bin.labels <- do.call(mergeLODBinLabels, list(scan.bin.labels, perm.bin.labels))
    thresholds <- parseLODBinLabels(bin.labels)
    
    # Get number of bins in scanone result (i.e. last scan bin index).
    num.scan.bins <- last.scan.index <- dim(scan.bins)[2]
    
    # Get number of bins in scanone permutations.
    num.perm.bins <- dim(object)[2]
    
    # Set number of bins from greater of scanone and permutations.
    num.bins <- max(num.scan.bins, num.perm.bins)
    
    # Pad scanone bins if needed.
    if ( num.scan.bins < num.bins ) {
        scan.bins <- padBins(scan.bins, num.bins)
    }
    
    # Pad permutation bins if needed.
    if ( num.perm.bins < num.bins ) {
        object <- padBins(object, num.bins)
    }
    
    # Get counts of above-threshold loci in scanone and permutations.
    perm.counts <- scan.counts <- rep_len(0, num.bins)
    for ( i in seq_along(thresholds) ) {
        perm.counts[i] <- sum(object[, i:num.bins, ])
        scan.counts[i] <- sum(scan.bins[, i:num.bins, ])
    }
    
    # Get mean permutation counts.
    mean.perm.counts <- perm.counts / num.perms
    
    # Get FDRs at each fixed LOD threshold.
    fdrs <- mean.perm.counts / scan.counts
    
    # Get indices of positive finite FDR estimates.
    # NB: this leaves FDR estimates greater than one, but we filter those later.
    util.indices <- which( ! is.na(fdrs) & is.finite(fdrs) & fdrs > 0 )
    
    # If there are useable FDR values, for each requested FDR, find the
    # closest-matching FDR estimate and its corresponding LOD threshold..
    if ( length(util.indices) > 0 ) {
        
        # Get useable FDR estimates.
        util.fdrs <- fdrs[util.indices]
        
        # Warn if FDR estimates are noisy.
        # NB: at low scan counts, the distribution of significant scanone
        # results is relatively noisy compared to the distribution of
        # significant permutations. This can mean that as the threshold
        # decreases, the proportional increase in mean permutation counts can
        # be greater than the proportional increase in scan counts, which can
        # actually result in a higher threshold having an artefactually higher
        # FDR estimate. To circumvent this issue, for each specified FDR, where
        # there are multiple closest-matching FDR estimates, we take the one
        # associated with the highest LOD threshold. This may result in a
        # conservative overestimate of LOD threshold, but that is often
        # preferable to an underestimate.
        if ( is.unsorted( rev(util.fdrs) ) ) {
            warning("FDR estimates are noisy for scanone with ", length(lodcolumns),
                " phenotypes - more accurate LOD thresholds may be obtained with",
                " more phenotypes, or by setting stringency with 'alpha' instead")
        }
        
        # Init closest-match FDR values, meaning those FDR estimates that
        # are closest available to the requested false-discovery rates.
        prox.fdrs <- rep(NA_real_, length(fdr))
        
        # Search for closest-match FDR values.
        for ( i in seq_along(fdr) ) {
            
            prox.fdr <- NA_real_
            
            for ( util.fdr in rev(util.fdrs) ) {
                
                if ( util.fdr <= fdr[i] ) {
                    
                    # If FDR estimate is less than or equal to requested FDR,
                    # check FDR estimate for next highest LOD threshold..
                    prox.fdrs[i] <- util.fdr
                    
                } else {
                    
                    # ..otherwise the previously-obtained FDR estimate was the
                    # closest match, so stop searching. If no valid FDR estimate
                    # was obtained that is less than or equal to the requested
                    # FDR, and this estimate is valid, accept this FDR estimate
                    # as being the closest-match FDR.
                    # NB: here we filter FDR estimates greater than one.
                    if ( is.na(prox.fdr) && util.fdr < 1.0 ) {
                        prox.fdrs[i] <- util.fdr
                    }
                    break
                }
            }
        }
        
        # Warn if LOD threshold saturation at high FDR and low LOD threshold.
        overestimates <- fdr[ ! is.na(prox.fdrs) & fdr > max(prox.fdrs) ]
        if ( length(overestimates) > 0 ) {
            warning("LOD threshold for FDR value(s) (", toString(overestimates),
                ") may be an overestimate - please choose lower FDR value(s)")
        }
        
        # Warn if LOD threshold saturation at low FDR and high LOD threshold.
        underestimates <- fdr[ ! is.na(prox.fdrs) & fdr < min(prox.fdrs) ]
        if ( length(underestimates) > 0 ) {
            if ( num.scan.bins <= num.perm.bins ) {
                warning("LOD threshold for FDR value(s) (", toString(underestimates),
                    ") is maximum possible with this scanone LOD profile",
                    " - may be an underestimate")
            } else {
                warning("LOD thresholds for FDR value(s) (", toString(underestimates),
                    ") may be an underestimate - more permutations may help")
            }
        }
        
        # Warn if the closest-match FDR differs from the
        # requested FDR by an order of magnitude or more.
        oom.indices <- which( abs( log10(fdr) - log10(prox.fdrs) ) >= 1.0 )
        for ( i in oom.indices ) {
            warning("requested FDR (", fdr[i], ") is very different from its",
                " closest-matching FDR (", prox.fdrs[i], ")")
        }
        
        # Get non-redundant set of valid closest-match FDR values.
        prox.fdrs <- unique( prox.fdrs[ ! is.na(prox.fdrs) ] )
        
        # If any closest-match FDR values found, get corresponding thresholds..
        if ( length(prox.fdrs) > 0 ) {
            
            # Get rightmost index for each distinct closest-match FDR value.
            prox.rindices <- unique(match(prox.fdrs, rev(fdrs)))
            prox.indices <- length(fdrs) - prox.rindices + 1
            
            # Get thresholds and FDR estimates for closest-match FDR values.
            thresholds <- thresholds[prox.indices]
            fdrs <- fdrs[prox.indices]
            
        } else { # ..otherwise clear LOD thresholds and FDR estimates.
            
            thresholds <- numeric()
            fdrs <- numeric()
        }
        
    } else { # ..otherwise clear LOD thresholds and FDR estimates.
        
        thresholds <- numeric()
        fdrs <- numeric()
    }
    
    # If no FDR values found, set requested thresholds to NA values.
    if ( length(fdrs) == 0 ) {
        warning("LOD thresholds not found for FDR values - '", toString(fdr), "'")
        thresholds <- rep_len(NA_real_, length(fdr))
        fdrs <- fdr
    }
    
    # Set significance levels to ensure no two FDR values clash.
    if ( length(fdrs) > 1 ) {
        fdr.diffs <- abs( diff(fdrs) )
        fdr.delta <- pmin(c(fdr.diffs, Inf), c(Inf, fdr.diffs), na.rm=TRUE)
        digits <- max(3, abs(floor(log10(fdr.delta))) )
    } else {
        digits <- 3
    }
    
    # Get percent FDR values, rounding to the given number of significant digits.
    pct.fdrs <- paste0(signif(fdrs, digits=digits) * 100, '%')
    
    result <- matrix(thresholds, nrow=length(fdrs), ncol=1,
        dimnames=list(pct.fdrs, 'lod'))
    
    attr(result, 'n.perm') <- num.perms
    
    class(result) <- 'summary.scanonebins'
    
    return(result)
}

# End of scanonebins.R #########################################################