# Start of cross.R #############################################################

# TODO: improve identification of tetrads.

# Cross Utilities --------------------------------------------------------------
#' Utilities for handling \pkg{R/qtl} \code{cross} objects.
#' 
#' @description Cross utilities can query or modify an existing \pkg{R/qtl} 
#' \code{cross} object.
#'   
#' @details An \pkg{R/qtl} \code{cross} object contains phenotype data for a set 
#' of individual samples, along with genotype data for the same samples for 
#' multiple marker loci across one or more sequences.
#' 
#' These utilities include functions to pull sequences
#' (\code{\link{pull.chr}}), individual sample IDs (\code{\link{pull.ind}}), or
#' alleles (\code{\link{pull.alleles}}). They also include functions to get
#' column indices within an \pkg{R/qtl} \code{cross$pheno} \code{data.frame}, 
#' for the ID column (\code{\link{getIdColIndex}}) and phenotype columns 
#' \code{\link{getPhenoColIndices}}.
#' 
#' There are functions to infer tetrad indices 
#' (\code{\link{inferTetradIndices}}) or strain indices 
#' (\code{\link{inferStrainIndices}}) for a \code{cross} object, based on 
#' sample IDs (if available) and genotype data.
#' 
#' A \code{cross} object can be modified by permuting either its phenotype or 
#' genotype data (\code{\link{permCross}}). When permuting phenotypes, it may
#' also be necessary to permute any associated data (e.g. covariates) in the 
#' same way. In such cases, \code{\link{permIndices}} can be used to generate
#' a set of permutation indices, which can in turn be passed to 
#' \code{\link{permCross}} or used to permute associated data.
#' 
#' There are several functions for handling a \code{cross} with a time-series 
#' phenotype, such that phenotypes have been observed at regular time points, 
#' and successive phenotype columns are named for the respective observation 
#' times (e.g. '0.0' for an observation at time zero). These allow for the time 
#' step to be inferred automatically (\code{\link{inferTimeStep}}), for the 
#' phenotype data to be padded in gaps of the time series 
#' (\code{\link{padTimeSeries}}), and for the phenotype values to 
#' be interpolated in time series gaps (\code{\link{interpTimeSeries}}).
#' 
#' @docType package
#' @name Cross Utilities
NULL

# getIdColIndex ----------------------------------------------------------------
#' Get sample ID column index.
#' 
#' Get the index of the sample ID column in the \code{cross} phenotype 
#' \code{data.frame}.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#'     
#' @return Sample ID column index.
#' 
#' @export
#' @family cross utilities
#' @rdname getIdColIndex
getIdColIndex <- function(cross) {
    
    stopifnot( 'cross' %in% class(cross) )
    
    id.col <- which( tolower( names(cross$pheno) ) == 'id' )
    
    if ( length(id.col) == 0 ) {
        id.col <- NULL
    } else if ( length(id.col) > 1 ) {
        stop("multiple ID columns found")
    }
    
    return(id.col)
}

# getFlankingPhenoColIndices ---------------------------------------------------
#' Get flanking phenotype column indices.
#' 
#' Get the indices of the phenotype columns in the \code{cross} phenotype data 
#' frame that are flanking the specified index, excluding special columns such 
#' as the sample ID column.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param i Column index for which the flanking phenotype column indices 
#' should be returned.
#'     
#' @return Flanking phenotype column indices. If the specified index refers to
#' the first or last column, an NA value is returned for the lower or upper 
#' flanking index, respectively.
#' 
#' @keywords internal
#' @rdname getFlankingPhenoColIndices
getFlankingPhenoColIndices <- function(cross, i) {
    
    stopifnot( isSinglePositiveNumber(i) )
    
    # Get indices of phenotype columns in cross.
    pheno.col <- getPhenoColIndices(cross)
    
    if ( i < 1 || i > ncol(cross$pheno) ) {
        stop("cross phenotype column index out of bounds - '", i, "'")
    }
    
    # Get lower flanking phenotype column index, or NA value if none available.
    lower.diff <- i - pheno.col
    if ( any(lower.diff > 0) ) {
        lower.index <- which( lower.diff == min(lower.diff[ lower.diff > 0 ]) )
    } else {
        lower.index <- NA
    }
    
    # Get upper flanking phenotype column index, or NA value if none available.
    upper.diff <- pheno.col - i
    if ( any(upper.diff > 0) ) {
        upper.index <- which( upper.diff == min(upper.diff[ upper.diff > 0 ]) )
    } else {
        upper.index <- NA
    }
    
    return( c(lower.index, upper.index) )
}

# getPhenoColIndices -----------------------------------------------------------
#' Get phenotype column indices.
#' 
#' Get the indices of the phenotype columns in the \code{cross} phenotype data 
#' frame, excluding special columns such as the sample ID column.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Vector indicating the phenotype indices to return. This can 
#' be an integer vector of phenotype column indices with respect to the columns 
#' of the \code{cross} phenotype \code{data.frame}, or a character vector that 
#' contains phenotype IDs or their syntactically valid names. If none are 
#' specified, all phenotype column indices are returned.
#'     
#' @return Phenotype column indices.
#' 
#' @export
#' @family cross utilities
#' @rdname getPhenoColIndices
getPhenoColIndices <- function(cross, pheno.col=NULL) {

    stopifnot( 'cross' %in% class(cross) )
    
    reserved.indices <- which( tolower( names(cross$pheno) ) %in% 
        const$reserved.phenotypes )
    
    if ( ! is.null(pheno.col) ) {
        
        stopifnot( length(pheno.col) > 0 )
        stopif( anyNA(pheno.col) )
        
        if ( is.character(pheno.col) ) {
            
            indices <- vector('integer', length(pheno.col))
            pheno.names <- names(cross$pheno)
            
            for ( i in getIndices(pheno.col) ) {
                
                p <- pheno.col[i]
                
                pheno.index <- which( pheno.names == make.names(p) )
                
                if ( length(pheno.index) == 0 ) {
                    stop("index not found for pheno.col - '", p, "'")
                }
                
                if ( length(pheno.index) > 1 ) {
                    stop("multiple indices found for pheno.col - '", p, "'")
                }
                
                if ( pheno.index %in% reserved.indices ) {
                    stop("pheno.col heading matches reserved phenotype - '", p, "'")
                }
                
                indices[i] <- pheno.index
            }
            
        } else if ( is.numeric(pheno.col) && all( isWholeNumber(pheno.col) ) ) {
            
            exrange <- pheno.col[ pheno.col < 1 | pheno.col > ncol(cross$pheno) ]
            if ( length(exrange) > 0 ) {
                stop("pheno.col indices out of range - '", toString(exrange), "'")
            }
            
            reserved <- pheno.col[ pheno.col %in% reserved.indices ]
            if ( length(reserved) > 0 ) {
                stop("pheno.col indices point to reserved phenotype headings - '", toString(reserved), "'")
            }
            
            indices <- pheno.col 
            
        } else {
            
            stop("pheno.col must contain phenotype indices or phenotype names")
        }
        
    } else {
        
        indices <- getColIndices(cross$pheno)
        indices <- indices[ ! indices %in% reserved.indices ]
    }
    
    return(indices)
}

# inferLodColNames -------------------------------------------------------------
#' Infer \code{scanone} result LOD column names.
#' 
#' Given a \code{cross} and set of phenotypes, get the expected the LOD column 
#' names of the corresponding scanone result.
#'
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param pheno.col Vector indicating the phenotypes for which LOD column names 
#' should be inferred.. This can be an integer vector of phenotype column 
#' indices with respect to the columns of the \code{cross} phenotype 
#' \code{data.frame}, or a character vector that contains phenotype IDs or their 
#' syntactically valid names. If none are specified, LOD names are inferred for 
#' all \code{cross} phenotypes.
#'  
#' @return Vector of the LOD column headings that would be expected in a 
#' \code{scanone} result for the given \code{cross} and phenotypes.
#' 
#' @keywords internal
#' @rdname inferLodColNames
inferLodColNames <- function(cross, pheno.col=NULL) {
    
    pheno.col <- getPhenoColIndices(cross, pheno.col)
    
    if ( length(pheno.col) == 1 ) {
        lodcol.names <- 'lod'
    } else if ( length(pheno.col) > 1 ) {
        lodcol.names <- colnames(cross$pheno)[pheno.col]
    } else {
        stop("cannot infer LOD column names - no phenotypes specified")
    }
    
    return(lodcol.names)
}

# inferStrainIndices -----------------------------------------------------------
#' Infer sample strain indices.
#'
#' Infers strain indices of a set of samples from the given object. A character 
#' vector of sample IDs can be passed as an argument, in which case the strain 
#' indices will be inferred from each set of identical consecutive sample IDs. 
#' If a \code{cross} object is passed as an argument, strain indices will be 
#' inferred from sample IDs, if present, or otherwise by comparison of sample 
#' genotypes.
#'
#' @param x An \pkg{R/qtl} \code{cross} object, or a character vector of sample
#' IDs.
#'  
#' @return Vector of sample strain indices inferred from the input object.
#' 
#' @export
#' @family cross utilities
#' @rdname inferStrainIndices
inferStrainIndices <- function (x) {
    UseMethod('inferStrainIndices', x)
}

# inferStrainIndices.character -------------------------------------------------
#' @export
#' @rdname inferStrainIndices
inferStrainIndices.character <- function(x) {
    
    runs <- rle(x)
    
    if ( anyDuplicated(runs$values) ) {
        stop("non-consecutive replicate sample IDs")
    }
    
    runs$values <- 1:length(runs$values)
    
    return( inverse.rle(runs) )
}

# inferStrainIndices.cross -----------------------------------------------------
#' @export
#' @rdname inferStrainIndices
inferStrainIndices.cross <- function(x) {
    
    # For each sample, combine genotypes into a single string.
    geno <- apply(qtl::pull.geno(x), 1, paste0, collapse='')
    
    # Get cross sample IDs.
    cross.info <- attr(x, 'info')
    if ( ! is.null(cross.info) ) {
        compareCrossInfo(x, cross.info)
        sample.ids <- getSamples(cross.info)
    } else {
        sample.ids <- pull.ind(x)
    }
    
    stopifnot( length(sample.ids) > 0 )
    
    # If cross has sample IDs, infer strain indices from these, then 
    # check genotypes are identical in the inferred strain replicates..
    if ( ! is.null(sample.ids) ) {
        
        strain.indices <- inferStrainIndices(sample.ids)
        
        for ( s in unique(strain.indices) ) {
            strain.geno <- geno[ which( strain.indices == s ) ]
            if ( any(strain.geno != strain.geno[1]) ) {
                stop("genotype mismatch in replicate samples")
            }
        } 
        
    } else { # ..otherwise infer strain indices directly from genotypes.
        runs <- rle(geno)
        runs$values <- 1:length(runs$values)
        strain.indices <- inverse.rle(runs)
    }
    
    return(strain.indices)
}

# inferTetradIndices -----------------------------------------------------------
#' Infer sample tetrad indices.
#'
#' @description Infers tetrad indices of a set of samples from the given object. 
#' A character vector of sample IDs can be passed as the first argument, in 
#' which case the tetrad indices will be inferred from the sample IDs themselves.
#' If sample IDs have a numeric suffix, these are taken as segregant numbers 
#' and used to infer which tetrad each sample ID might belong to. 
#'   
#' If a \code{cross} is the first argument, tetrad indices will be inferred 
#' from sample IDs, if present, or otherwise by comparison of sample genotypes.
#'   
#' @param x An \pkg{R/qtl} \code{cross} object, or a character vector of sample IDs.
#'  
#' @return Vector of sample tetrad indices inferred from the input object. 
#' Returns \code{NULL} if tetrad pattern could not be found.
#' 
#' @export
#' @family cross utilities
#' @rdname inferTetradIndices
inferTetradIndices <- function(x) {
    UseMethod('inferTetradIndices', x)
}

# inferTetradIndices.character -------------------------------------------------
#' @export
#' @rdname inferTetradIndices
inferTetradIndices.character <- function(x) {

    # Tetrad indices corresponding to the input sample IDs.
    sample.tindices <- NULL
    
    # Get 'exemplar' sample IDs.
    exemplar.ids <- unique(x)
    
    # Get indices of exemplars with respect to the vector of sample IDs.
    exemplar.sindices <- sapply(exemplar.ids, match, x)
    
    # Get indices of replicates with respect to the vector of sample IDs.
    terminal.sindices <- c(sapply(exemplar.sindices[2:length(exemplar.sindices)], 
        function(e) e - 1), length(x) )
    replicate.sindices <- mapply(function(i, j) i:j, exemplar.sindices, 
        terminal.sindices, SIMPLIFY=FALSE)
    
    # Indices of exemplars with respect to the 'notional' exemplar sample IDs, 
    # in the ideal situation where every sample is present in every tetrad. 
    # Notional exemplar indices are taken from the numeric/alphanumeric suffix 
    # of each sample ID. This is used to place exemplar samples in tetrads,
    # and can work even if there are missing samples.
    exemplar.nindices <- NULL
    
    # If exemplars have numeric suffixes, try to assign notional indices from these..
    if ( all( grepl(const$pattern$tetrad['numeric'], exemplar.ids) ) ) {
        
        # Get sample numbers from numeric sample ID suffixes.
        sample.numbers <- as.integer( sub(const$pattern$tetrad['numeric'], '\\1', exemplar.ids) )
        
        # Test if sample numbers are ordered.
        samples.ordered <- ! is.unsorted(sample.numbers)
        
        # If sample numbers are sorted, get notional indices from the sample 
        # numbers, as offset by the value of the first sample number.
        if (samples.ordered) {
            exemplar.nindices <- sample.numbers - sample.numbers[1] + 1
        }
    
    # ..otherwise if exemplars have alphanumeric suffixes, try to assign notional indices from those.    
    } else if ( all( grepl(const$pattern$tetrad['alphanumeric'], exemplar.ids) ) ) {
        
        # Get tetrad numbers from numeric part of sample ID suffixes.
        tetrad.numbers <- as.integer( sub(const$pattern$tetrad['alphanumeric'], '\\1', exemplar.ids) )
        
        # Test if tetrad numbers are ordered.
        tetrads.ordered <- ! is.unsorted(tetrad.numbers)
        
        # Get tetrad sample labels, assign a label index to each (e.g. 1 for 'A', 2 for 'B').
        sample.labels <- toupper( sub(const$pattern$tetrad['alphanumeric'], '\\2', exemplar.ids) )
        exemplar.lindices <- sapply(sample.labels, function(sample.label) 
            which( const$tetrad.sample.labels == sample.label ), USE.NAMES=FALSE)
        
        # Test if sample labels are ordered within tetrads.
        labels.ordered <- all( sapply(unique(tetrad.numbers), function(i) 
            ! is.unsorted(exemplar.lindices[ which( tetrad.numbers == i ) ]) ) )
        
        # If tetrad and sample numbers are sorted, get notional indices from the tetrad 
        # and sample numbers, as offset by the value of the first tetrad number.
        if ( tetrads.ordered && labels.ordered ) {
            exemplar.tindices <- tetrad.numbers - tetrad.numbers[1] + 1
            exemplar.nindices <- 4 * (exemplar.tindices - 1) + exemplar.lindices
        }
    }
    
    # If notional exemplar indices were identified, use these to assign tetrad indices.
    if ( ! is.null(exemplar.nindices) ) {
        
        num.tetrads <- ceiling( max(exemplar.nindices) / 4 )
        
        expected.xindices <- lapply(1:num.tetrads, function(i) 
            { o <- 4*(i-1); seq(o + 1, o + 4) } )
        
        tetrad.xindices <- lapply(1:num.tetrads, function(i) 
            exemplar.nindices[ exemplar.nindices %in% expected.xindices[[i]] ])
        
        sample.tindices <- integer( length=length(x) )
        for ( i in 1:num.tetrads ) {
            tetrad.sindices <- unlist( replicate.sindices[ tetrad.xindices[[i]] ] )
            sample.tindices[tetrad.sindices] <- i
        }
    }    
    
    return(sample.tindices)
}

# inferTetradIndices.cross -----------------------------------------------------
#' @export
#' @rdname inferTetradIndices
inferTetradIndices.cross <- function(x) {
    
    # Tetrad indices corresponding to the input sample IDs.
    sample.tindices <- NULL
    
    # Get cross sample IDs.
    cross.info <- attr(x, 'info')
    if ( ! is.null(cross.info) ) {
        compareCrossInfo(x, cross.info)
        sample.ids <- getSamples(cross.info)
    } else {
        sample.ids <- pull.ind(x)
    }
    
    stopifnot( length(sample.ids) > 0 )
    
    # Set minimum fraction of tetradic genotypes required to infer tetrad indices.
    threshold <- 0.95
    
    # Infer strain indices, as replicates will be included in each tetrad.
    strain.indices <- inferStrainIndices(x)
    
    # Get indices of exemplars with respect to tthe set of samples.
    exemplar.sindices <- sapply(unique(strain.indices), match, strain.indices)
    
    # Get indices of replicates with respect to the set of samples.
    terminal.sindices <- c(sapply(exemplar.sindices[2:length(exemplar.sindices)], 
        function(e) e - 1), length(strain.indices) )
    replicate.sindices <- mapply(function(i, j) i:j, exemplar.sindices, 
        terminal.sindices, SIMPLIFY=FALSE)
    
    # Tetrad indices with respect to the set of exemplar 
    # sample indices. NB: no need to look at full sample 
    # data, as replicate samples have identical genotypes.
    tetrad.xindices <- NULL
    
    # If cross has sample IDs, try to infer tetrad indices from these.
    if ( ! is.null(sample.ids) ) {
        id.tindices <- inferTetradIndices(sample.ids)
        num.tetrads <- max(id.tindices)
        exemplar.tindices <- id.tindices[exemplar.sindices]
        tetrad.xindices <- lapply(1:num.tetrads, function(i) 
            which( exemplar.tindices == i ) )
    } 
    
    # If we failed to infer tetrad exemplar indices from sample IDs, 
    # fall back to assuming that samples form complete tetrads.
    if ( is.null(tetrad.xindices) ) {
        num.tetrads <- length(exemplar.sindices) / 4
        if ( isWholeNumber(num.tetrads) ) {
            tetrad.xindices <- lapply(1:num.tetrads, function(i) 
            { o <- 4*(i-1); seq(o + 1, o + 4) } )
        }
    }
    
    # If tetrad indices identified from sample IDs or assumed from number of 
    # exemplars, check this against genotype data and get tetrad indices 
    # corresponding to samples. 
    if ( ! is.null(tetrad.xindices) ) {
        
        # Calculate total number of possible tetrad genotypes.
        total.tetrad.genotypes <- sum(lengths(tetrad.xindices) > 0) * qtl::totmar(x)
        
        total.counted <- total.tetradic <- 0
        
        geno <- qtl::pull.geno(x)
        
        markers <- colnames(geno)
        
        # For each candidate tetrad, count tetradic genotypes.
        for ( i in 1:num.tetrads) {
            
            tetrad.exemplar.sindices <- exemplar.sindices[ tetrad.xindices[[i]] ]
            
            # If we have genotypes for a complete tetrad, check if appears to be tetradic.
            if ( length(tetrad.exemplar.sindices) == 4 ) {
                
                tetrad.geno <- geno[tetrad.exemplar.sindices, ]
                
                tetrad.tables <- lapply(markers, function(m) table(tetrad.geno[, m]))
                
                geno.counted <- sapply(tetrad.tables, sum) == 4
                
                geno.homozygous <- ! geno.counted | sapply(tetrad.tables, function(x) all(x == 4))
                
                # Create run-length encoding of homozygous genotype mask.., 
                runs <- rle(geno.homozygous)
                # ..set internal runs of homozygosity as FALSE..
                runs$values[ -c(1, length(runs$values)) ] <- FALSE
                # ..and then recreate the original homozygous genotype mask.
                geno.homozygous <- inverse.rle(runs)
                
                geno.counted <- geno.counted & ! geno.homozygous
                
                geno.tetradic <- geno.counted & sapply(tetrad.tables, function(x) all(x == 2))
                
                total.tetradic <- total.tetradic + sum(geno.tetradic)
                total.counted <- total.counted + sum(geno.counted)
            }
        }
        
        # If enough genotypes checked and enough of those 
        # were tetradic, get tetrad indices of samples.
        if ( total.counted > 0.5 * total.tetrad.genotypes ) {
            
            if ( total.tetradic / total.counted >= threshold ) {
                
                sample.tindices <- integer( length=qtl::nind(x) )
                
                for ( i in 1:num.tetrads ) {
                    tetrad.sindices <- unlist( replicate.sindices[ tetrad.xindices[[i]] ] )
                    sample.tindices[tetrad.sindices] <- i
                }
            } 
        }  
    }
    
    return(sample.tindices)
}  

# inferTimeStep ----------------------------------------------------------------
#' Infer time step of time-series phenotypes.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param allow.gaps Allow gaps in time series, provided that any gaps are a 
#' multiple of the inferred time step. 
#' @param tol Tolerance for time step equality.
#'     
#' @return Inferred time step. Returns \code{NULL} if time step could not be inferred.
#' 
#' @export
#' @family cross utilities
#' @rdname inferTimeStep
inferTimeStep <- function(cross, allow.gaps=TRUE, tol=1e-5) {
    
    stopifnot( isBOOL(allow.gaps) )
    
    # Get time-series phenotypes.
    cross.info <- attr(cross, 'info')
    if ( ! is.null(cross.info) ) {
        compareCrossInfo(cross, cross.info)
        phenotypes <- getPhenotypes(cross.info)
        times <- as.numeric(phenotypes)
    } else {
        pheno.col <- getPhenoColIndices(cross)
        phenotypes <- colnames(cross$pheno)[pheno.col]
        times <- makeNumbers(phenotypes)
    }
    
    stopifnot( length(phenotypes) > 0 )
    
    # One time point cannot form a series.
    if ( length(times) == 1 ) {
        return(NULL)
    }
    
    # Time points must be valid non-negative numbers.
    if ( ! all ( isNonNegativeNumber(times) ) ) {
        return(NULL)
    }
    
    # Times must be strictly monotonically increasing.
    if ( is.unsorted(times, strictly=TRUE) ) {
        return(NULL)
    }
    
    # Get vector of time steps.
    time.steps <- diff(times)
    
    # Infer size of time step.
    time.step <- inferStepSize(time.steps, tol=tol)
    
    if ( is.null(time.step) ) {
        return(NULL)
    }
    
    # Time steps must be consistent.
    if ( any(time.step - time.steps > tol) ) {
        return(NULL)
    }
    
    # Identify gaps where difference between two times deviates from time step.
    gap.indices <- which( time.steps - time.step > tol )
    
    # If there are any gaps, verify that gaps are allowed and 
    # that all gaps are a multiple of the apparent time step.
    if ( length(gap.indices) > 0 ) {  
        
        if ( ! allow.gaps ) {
            return(NULL)
        }
        
        step.counts <- time.steps[gap.indices] / time.step - 1
        
        if ( any( ! isWholeNumber(step.counts, tol=tol) ) ) {
            return(NULL)
        }
    }
    
    return(time.step)
}

# interpTimeSeries -------------------------------------------------------------
#' Interpolate gaps in a time-series.
#' 
#' Interpolate values for gaps in time-series phenotypes of a \code{cross} 
#' object.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param tol Tolerance for time step equality.
#'     
#' @return The input \code{cross} object is returned with gaps in time-series 
#' phenotypes filled with values interpolated from the gap endpoints.
#'  
#' @export
#' @family cross utilities
#' @rdname interpTimeSeries
interpTimeSeries <- function(cross, tol=1e-5) {
    
    cross <- padTimeSeries(cross, tol=tol)
    
    # Get indices of phenotype columns in cross.
    pheno.col <- getPhenoColIndices(cross)
    
    # Get time-series phenotypes.
    cross.info <- attr(cross, 'info')
    if ( ! is.null(cross.info) ) {
        compareCrossInfo(cross, cross.info)
        phenotypes <- getPhenotypes(cross.info)
        times <- as.numeric(phenotypes)
    } else {
        phenotypes <- colnames(cross$pheno)[pheno.col]
        times <- makeNumbers(phenotypes)
    }
    
    numeric.phenotypes <- sapply(pheno.col, 
        function (i) is.numeric(cross$pheno[[i]]) )
    
    if ( ! all(numeric.phenotypes) ) {
        stop("cannot interpolate non-numeric phenotypes")
    }
    
    if ( anyNA(cross$pheno[, pheno.col[1] ]) ) {
        stop("cannot interpolate data at start of time series")
    }
    
    if ( anyNA(cross$pheno[, pheno.col[ length(pheno.col) ] ]) ) {
        stop("cannot interpolate data at end of time series")
    }
    
    for ( i in 1:nrow(cross$pheno) ) {
        
        # Get time-series values for sample (including blanks).
        time.series <- as.numeric( cross$pheno[i, pheno.col] )
        
        # Get runs of blank values, then get list of vectors, such that each  
        # vector contains the indices of a consecutive set of blank columns.
        blank.runs <- rle( is.na(time.series) )
        J <- cumsum(blank.runs$lengths)
        I <- c( 1, sapply(J[1:(length(J)-1)], function(j) j + 1) )
        b <- blank.runs$values == TRUE
        blank.index.list <- mapply(function(i, j) i:j, I[b], J[b], SIMPLIFY=FALSE)
        
        # We know from the code immediately above that consecutive runs of blank 
        # values are grouped together, and we have checked above that none of 
        # these run off either end of the time series. Therefore, we can safely 
        # take the flanking columns as our interpolation endpoints (below).
        
        for ( blank.indices in blank.index.list ) {
            
            # Take values from flanking columns as interpolation endpoints.
            x <- c(blank.indices[1] - 1, blank.indices[ length(blank.indices) ] + 1)
            
            # Get time series values for this blank region.
            y <- time.series[x]
            
            # Interpolate phenotype values for this blank region.
            interpolation <- approx(x, y, blank.indices, 'linear')
            
            # Set blank values from interpolation results.
            time.series[blank.indices] <- interpolation$y
        }
        
        # Set time-series values for sample (with blanks interpolated).
        cross$pheno[i, pheno.col] <- time.series
    }
    
    return(cross)
}

# makeCross --------------------------------------------------------------------
#' Make an \pkg{R/qtl} \code{cross} object.
#' 
#' @param geno A \code{geno} object.
#' @param pheno A \code{pheno} object.
#' 
#' @return An \pkg{R/qtl} \code{cross} object.
#' 
#' @keywords internal
#' @rdname makeCross
makeCross <- function(geno, pheno) {
    
    # TODO: handle replicate phenotype measurements
    
    stopifnot( 'geno' %in% class(geno) )
    stopifnot( 'pheno' %in% class(pheno) )
    
    # Get alleles of geno.
    alleles <- attr(geno, 'alleles')
    attr(geno, 'alleles') <- NULL
    
    # Get CrossInfo for geno data.
    geno.info <- attr(geno, 'info')
    compareCrossInfo(geno, geno.info)
    attr(geno, 'info') <- NULL

    # Get CrossInfo for pheno data.
    pheno.info <- attr(pheno, 'info')
    compareCrossInfo(pheno, pheno.info)
    attr(pheno, 'info') <- NULL
    
    if ( ! hasSampleIDs(geno.info) ) {
        stop("cannot make cross - no genotype sample IDs")
    }
    
    if ( ! hasSampleIDs(pheno.info) ) {
        stop("cannot make cross - no phenotype sample IDs")
    }
    
    if ( any(getSamples(geno.info) != getSamples(pheno.info)) ) {
        stop("cannot make cross - sample mismatch")
    }
    
    if ( is.na(geno.info@crosstype) ) {
        stop("cannot make cross - unknown cross type")
    }
    
    # Create CrossInfo object from genotype and phenotype info.
    cross.info <- methods::new('CrossInfo', seq=geno.info@seq,
        pheno     = pheno.info@pheno,
        markers   = geno.info@markers,
        samples   = geno.info@samples,
        alleles   = geno.info@alleles,
        genotypes = geno.info@genotypes,
        crosstype = geno.info@crosstype
    )
    
    # Prepare genotype and phenotype objects.
    class(pheno) <- 'data.frame'
    class(geno) <- 'list'
    
    # Create cross list.
    cross <- list(geno=geno, pheno=pheno)
    
    # Set cross alleles.
    attr(cross, 'alleles') <- alleles
    
    # Set cross info.
    attr(cross, 'info') <- cross.info
    
    # Set cross class.
    class(cross) <- c(cross.info@crosstype, 'cross')
    
    return(cross)
}

# padTimeSeries ----------------------------------------------------------------
#' Pad gaps in a time-series.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param tol Tolerance for time step equality.
#'     
#' @return The input \code{cross} object is returned with all gaps in 
#' time-series phenotypes padded with NA values.
#'  
#' @export
#' @family cross utilities
#' @rdname padTimeSeries
padTimeSeries <- function(cross, tol=1e-5) {
    
    # Get indices of phenotype columns in cross.
    pheno.col <- getPhenoColIndices(cross)
    
    # Get index of ID column in cross.
    id.col <- getIdColIndex(cross)
    
    # Get ID column name.
    id.col.name <- colnames(cross$pheno)[id.col]
    
    # If ID column found, get flanking phenotypes. The ID column  
    # is reinserted between these in the padded phenotype data.
    if ( ! is.null(id.col) ) {
        flanking.indices <- getFlankingPhenoColIndices(cross, id.col)
        flanking.phenames <- colnames(cross$pheno)[flanking.indices]
    }
    
    # Get phenotype data.
    pheno <- qtl::pull.pheno(cross, pheno.col)
    
    # Get phenotype data class.
    pheno.class <- unique( sapply(pheno, class) )
    
    if ( length(pheno.class) != 1 ) {
        stop("inconsistent phenotype data types")
    }
    
    # Get missing value for this phenotype data.
    missing.value <- getMissingValueFromClassS3(pheno.class)
    
    # Get phenotypes and sample IDs.
    cross.info <- attr(cross, 'info')
    if ( ! is.null(cross.info) ) {
        compareCrossInfo(cross, cross.info)
        phenotypes <- getPhenotypes(cross.info)
        times <- as.numeric(phenotypes)
        sample.ids <- getSamples(cross.info)
    } else {
        phenotypes <- colnames(cross$pheno)[pheno.col]
        times <- makeNumbers(phenotypes)
        sample.ids <- pull.ind(cross)
    }
    
    stopifnot( length(sample.ids) > 0 )
    
    # Time step is the mode of time differences.
    time.step <- inferTimeStep(cross, tol=tol)
    
    # Get differences between times.
    time.diff <- diff(times)
    
    # Identify gaps where difference between two times deviates from time step.
    gap.indices <- which( time.diff - time.step > tol )
    
    # Get number of such gaps.
    num.gaps <- length(gap.indices)
    
    # If there are gaps, pad them with NA values.
    if ( num.gaps > 0 ) {
        
        # Calculate step counts for all gaps.
        step.counts <- as.integer( round(time.diff[gap.indices] / time.step - 1) )
        
        # Set missing time points for each gap.
        gap.steps <- vector('list', num.gaps)
        for ( i in 1:num.gaps ) {
            g <- gap.indices[i]
            first.step <- times[g] + time.step
            last.step <- times[g + 1] - time.step
            gap.steps[[i]] <- seq(first.step, last.step, by=time.step)
        }
        
        # Get expected number of time points after gaps are padded.
        num.times <- length(times) + sum(step.counts)
        
        # Pad phenotype data and time points.
        padded.pheno <- data.frame( matrix(nrow=nrow(pheno), ncol=num.times), 
            stringsAsFactors=FALSE)
        padded.times <- rep(NaN, length=num.times)
        
        padded.pheno[, 1] <- pheno[, 1]
        padded.times[1] <- times[1]
        
        t <- 2
        
        for ( d in 1:length(time.diff) ) {
            
            if ( d %in% gap.indices ) {
                
                g <- which( gap.indices == d )
                gap.times <- gap.steps[[g]]
                
                for ( gap.time in gap.times ) {
                    
                    padded.pheno[, t] <- rep(missing.value, nrow(pheno))
                    padded.times[t] <- gap.time
                    t <- t + 1L
                }
            } 
            
            i <- d + 1L
            padded.pheno[, t] <- pheno[, i]
            padded.times[t] <- times[i]
            t <- t + 1L
        }
        
        # Set syntactically valid column names for padded phenotypes.
        colnames(padded.pheno) <- make.names(padded.times)
        
        # If phenotype has ID column, reinsert ID column 
        # as close as possible to its original position.
        if ( ! is.null(id.col) ) {
            
            if ( ! is.na(flanking.phenames[1]) ) {
                i <- which( colnames(padded.pheno) == flanking.phenames[1] ) + 1
            } else {
                i <- which( colnames(padded.pheno) == flanking.phenames[2] ) - 1
            }
            
            padded.pheno <- insertColumn(padded.pheno, i, col.name=id.col.name, data=sample.ids)
        }
        
        # Replace padded phenotype data.
        cross$pheno <- padded.pheno
        
        # If cross has CrossInfo, update phenotypes to reflect padded data.
        if ( ! is.null(cross.info) ) {
            cross.info <- setPhenotypes(cross.info, padded.times)
            attr(cross, 'info') <- cross.info
        }
    } 
    
    return(cross)
}

# permIndices ------------------------------------------------------------------
#' Generate permutation indices for a \code{cross} object.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#'     
#' @return Vector of sample indices that can be used to permute a 
#' \code{cross} object or associated data. 
#' 
#' @export
#' @rdname permIndices
permIndices <- function(cross) {
  
    stopifnot( 'cross' %in% class(cross) )
    
    num.samples <- nrow(cross$pheno)
    sample.indices <- 1:num.samples
    
    cross.info <- attr(cross, 'info')
    
    if ( ! is.null(cross.info) ) {
        compareCrossInfo(cross, cross.info)
        replicate.indices <- getStrainIndices(cross.info)
        tetrad.indices <- getTetradIndices(cross.info)
    } else {
        replicate.indices <- inferStrainIndices(cross)
        tetrad.indices <- inferTetradIndices(cross)
    }
    
    if ( ! is.null(tetrad.indices) ) {
        block.indices <- tetrad.indices
    } else {
        block.indices <- rep.int(1, num.samples)
    }
    
    blocks <- lapply( 1:max(block.indices), 
        function(s) which( block.indices == s ) )
    
    perm.indices <- rep(0, num.samples)
    
    for ( block in blocks ) {
        
        rep.freq <- table(replicate.indices[block])
        
        rep.freq <- sort(rep.freq, decreasing=TRUE)
        
        if( length(rep.freq) == 1 || rep.freq[1] > sum(rep.freq[2:length(rep.freq)]) ) {
            stop("cannot generate permutation indices - imbalanced replicates")
        }
        
        perm.block <- sample(sample.indices[block])

        while ( any(replicate.indices[perm.block] == replicate.indices[block]) ) {
            perm.block <- sample(sample.indices[block])
        }
        
        perm.indices[block] <- perm.block
    }
    
    return(perm.indices)
}
    
# permCross -----------------------------------------------------------------
#' Permute \code{cross} phenotype or genotype data.
#' 
#' @details Given an input \code{cross} object, this function permutes either
#' the phenotype or genotype data (but not both), and returns the permuted 
#' \code{cross} object. Note that when permuting phenotypes, any covariates
#' must also be permuted in the same way.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param perm.indices Permutation indices.
#' @param perm.pheno Permute phenotype data.
#' @param perm.geno Permute genotype data.
#'     
#' @return Permuted \code{cross} object.
#' 
#' @export
#' @family cross utilities
#' @rdname permCross
permCross <- function(cross, perm.indices=NULL, perm.pheno=TRUE, 
    perm.geno=FALSE) {

    stopifnot( 'cross' %in% class(cross) )
    stopifnot( isBOOL(perm.pheno) )
    stopifnot( isBOOL(perm.geno) )

    if ( ! xor(perm.pheno, perm.geno) ) {
        stop("permCross must permute either cross phenotypes or cross genotypes")
    }
    
    num.samples <- nrow(cross$pheno)
    
    if ( ! is.null(perm.indices) ) {
        stopifnot( length(perm.indices) == num.samples )
        stopifnot( all( perm.indices %in% 1:num.samples ) )
        stopif( anyDuplicated(perm.indices) )
    } else {
        perm.indices <- permIndices(cross)
    }
    
    if (perm.pheno) {
        
        pheno.col <- getPhenoColIndices(cross)
        
        cross$pheno[, pheno.col] <- cross$pheno[perm.indices, pheno.col]
    }
    
    if (perm.geno) {
        
        geno.seqs <- names(cross$geno)
        
        for ( geno.seq in geno.seqs ) {
            cross$geno[[geno.seq]]$data <- cross$geno[[geno.seq]]$data[perm.indices, ]
        }
        
        cross <- refreshCross(cross)
    }
    
    return(cross)
}

# pull.alleles -----------------------------------------------------------------
#' Pull alleles from \pkg{R/qtl} \code{cross}.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#'     
#' @return Vector of alleles in the \code{cross} object.
#' 
#' @export
#' @family cross utilities
#' @rdname pull.alleles
pull.alleles <- function(cross) {
    stopifnot( 'cross' %in% class(cross) )
    return( attr(cross, 'alleles') )
}

# pull.chr ---------------------------------------------------------------------
#' Pull sequences from \pkg{R/qtl} \code{cross}.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#'     
#' @return Vector of sequences in the \code{cross} object.
#' 
#' @export
#' @family cross utilities
#' @rdname pull.chr
pull.chr <- function(cross) {
    stopifnot( 'cross' %in% class(cross) )
    return( names(cross$geno) )
}

# pull.crosstype ---------------------------------------------------------------
#' Pull cross type from \pkg{R/qtl} \code{cross}.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#'     
#' @return Cross type of the \code{cross} object.
#' 
#' @export
#' @family cross utilities
#' @rdname pull.crosstype
pull.crosstype <- function(cross) {
    stopifnot( 'cross' %in% class(cross)  )
    return( class(cross)[1] )
}

# pull.ind ---------------------------------------------------------------------
#' Pull individual IDs from \pkg{R/qtl} \code{cross}.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#'     
#' @return Vector of individual IDs in the \code{cross} object.
#' 
#' @export
#' @family cross utilities
#' @rdname pull.ind
pull.ind <- function(cross) {
    
    stopifnot( 'cross' %in% class(cross) )
    
    id.col <- getIdColIndex(cross)
    
    if ( ! is.null(id.col) ) {
        ind <- as.character(cross$pheno[[id.col]])
    } else {
        ind <- NULL
    }
    
    return(ind)
}

# refreshCross -----------------------------------------------------------------
#' Refresh derived data attributes of \pkg{R/qtl} \code{cross}.
#' 
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param derived.attributes Names of derived data attributes to refresh. If 
#' none are specified, all derived data attributes in the \code{cross} are 
#' refreshed. 
#'  
#' @return Input \code{cross} object with derived data attributes refreshed.
#' 
#' @keywords internal
#' @rdname refreshCross
refreshCross <- function(cross, derived.attributes=NULL) {
    
    # This function is modelled on R/qtl function clean.cross; if that function
    # changes, this should be changed to reflect it. Also, this function relies 
    # on internal flag 'onlylod' for determining which function to use.
    
    stopifnot( 'cross' %in% class(cross) )
    
    known.derived <- c('argmax', 'draws', 'errorlod', 'prob', 'rf')
    
    if ( ! is.null(derived.attributes) ) {
        
        unknown <- derived.attributes[ ! derived.attributes %in% known.derived ]
        if ( length(unknown) > 0 ) {
            stop("cross has unknown derived data elements - '", 
                toString(unknown), "'")
        }
    
    } else {
        derived.attributes <- known.derived
    }
    
    geno.attr <- unique( as.vector( sapply( cross$geno, function (g) names(g) ) ) )
    
    if ( 'argmax' %in% derived.attributes && 'argmax' %in% geno.attr ) {
        
        arg.names <- c('step', 'off.end', 'error.prob', 'map.function', 'stepwidth')
        args <- lapply(arg.names, function (a) { unique( sapply( cross$geno, 
            function (g) attr(g$argmax, a) ) ) } )
        names(args) <- arg.names
        
        cross <- qtl::argmax.geno(cross, step=args$step, off.end=args$off.end, 
            error.prob=args$error.prob, map.function=args$map.function, 
            stepwidth=args$stepwidth)
    }
    
    if ( 'draws' %in% derived.attributes && 'draws' %in% geno.attr ) {
        
        n.draws <- unique( sapply( cross$geno, function (g) dim(g$draws)[3] ) )
        
        arg.names <- c('step', 'off.end', 'error.prob', 'map.function', 'stepwidth')
        args <- lapply(arg.names, function (a) { unique( sapply( cross$geno, 
            function (g) attr(g$draws, a) ) ) } )
        names(args) <- arg.names
        
        cross <- qtl::sim.geno(cross, n.draws=n.draws, step=args$step, 
            off.end=args$off.end, error.prob=args$error.prob, 
            map.function=args$map.function, stepwidth=args$stepwidth)
    }
    
    if ( 'errorlod' %in% derived.attributes && 'errorlod' %in% geno.attr ) {
        
        arg.names <- c('error.prob', 'map.function')
        args <- lapply(arg.names, function (a) { unique( sapply( cross$geno, 
            function (g) attr(g$errorlod, a) ) ) } )
        names(args) <- arg.names
        
        cross <- qtl::calc.errorlod(cross, error.prob=args$error.prob, 
            map.function=args$map.function)
    }
    
    if ( 'prob' %in% derived.attributes && 'prob' %in% geno.attr ) {
        
        arg.names <- c('step', 'off.end', 'error.prob', 'map.function', 'stepwidth')
        args <- lapply(arg.names, function (a) { unique( sapply( cross$geno, 
            function (g) attr(g$prob, a) ) ) } )
        names(args) <- arg.names
        
        cross <- qtl::calc.genoprob(cross, step=args$step, off.end=args$off.end, 
            error.prob=args$error.prob, map.function=args$map.function, 
            stepwidth=args$stepwidth)
    }
    
    if ( 'rf' %in% derived.attributes && 'rf' %in% names(cross) ) {
        
        lod.flag <- attr(cross$rf, 'onlylod')
        
        if ( ! is.null(lod.flag) && lod.flag ) {
            cross <- qtl::markerlrt(cross)
        } else {
            cross <- qtl::est.rf(cross)
        }
    }
    
    return(cross)
}

# End of cross.R ###############################################################