# Start of geno.R ##############################################################

# as.data.frame.geno -----------------------------------------------------------
#' Coerce \code{geno} object to \code{data.frame}.
#' 
#' @param x A \code{geno} object.
#' @param ... Unused arguments.
#' @param chr Vector of sequences for which genotype data should be included in 
#' the genotype \code{data.frame}. If none are specified, a genotype 
#' \code{data.frame} is returned for all available sequences.
#' @param digits If specified, round genetic map positions to the specified 
#' number of digits.
#' @param include.mapunit Include map unit information in map positions.
#' 
#' @return A \code{data.frame} corresponding to the input \code{geno} object.
#' 
#' @keywords internal
#' @method as.data.frame geno
#' @rdname as.data.frame.geno
as.data.frame.geno <- function(x, ..., chr=NULL, digits=NULL,
    include.mapunit=TRUE) {
    
    stopifnot( isBOOL(include.mapunit) )
    
    map.unit <- 'cM'
    
    # Get CrossInfo object.
    cross.info <- attr(x, 'info')
    
    compareCrossInfo(x, cross.info)
    
    # Get relevant CrossInfo.
    alleles <- getAlleles(cross.info)
    locus.ids <- getMarkers(cross.info)
    locus.names <- make.names(locus.ids)
    samples <- getSamples(cross.info)
    
    # Get specified sequences.
    chr <- subsetBySeq(names(x), chr)
    stopifnot( length(chr) > 0 )
    
    # Subset geno object by the specified chromosomes.
    x <- x[chr]
    
    # Pull map from geno object.
    map.table <- as.mapframe( lapply(x, function(obj)
        obj$map), map.unit=map.unit )
    
    # If digits specified, round map positions.
    if ( ! is.null(digits) ) {
        stopifnot( isSinglePositiveWholeNumber(digits) )
        map.table$pos <- round(map.table$pos, digits=digits)
    }
    
    # If including map unit, add map unit info to map table.
    if (include.mapunit) {
        map.table <- setPosColDataMapUnit(map.table, map.unit)
    }
    
    # Pull genotype matrix from geno object.
    geno.matrix <- do.call(cbind, lapply(x, function(obj) obj$data))
    
    # Replace encoded genotypes with actual genotype values.
    for ( i in getIndices(alleles) ) {
        geno.matrix[ geno.matrix == i ] <- alleles[i]
    }
    
    # Prepare map matrix.
    map.table <- insertColumn(map.table, col.index=1, 
        col.name='id', data=locus.ids)
    rownames(map.table) <- locus.names
    map.matrix <- t(map.table)
    
    # Prepare genotype matrix.
    dimnames(geno.matrix) <- list(samples, locus.names)
    
    # Bind map and genotype matrices into one data frame.
    geno.frame <- as.data.frame(rbind(map.matrix, geno.matrix), 
        stringsAsFactors=FALSE)
    
    return(geno.frame)
}

# as.geno ----------------------------------------------------------------------
#' Coerce to \code{geno} object.
#' 
#' @param from Object containing genotype data.
#' 
#' @return An \code{geno} object corresponding to the input object.
#' 
#' @keywords internal
#' @rdname as.geno
as.geno <- function(from) {
    UseMethod('as.geno', from)
}

# as.geno.data.frame -----------------------------------------------------------
#' @method as.geno data.frame
#' @rdname as.geno
as.geno.data.frame <- function(from) {
    
    crosstype <- 'bc' # TODO: support other cross types
    
    # Set putative map positions from third row.
    putative.positions <- as.character(from[3, ])
    
    # Try to get map unit from putative map positions.
    map.unit <- getMapUnitSuffix(putative.positions)
    
    if ( map.unit != 'cM' ) {
        stop("cross geno map positions must be in centiMorgans (e.g. '47 cM')")
    }
    
    # Denote map positions as detected if a map unit was identified.
    map.positions.found <- ! is.na(map.unit)
    
    # If map positions not detected, double-check putative map positions.
    if ( ! map.positions.found ) {
        
        # Get set of values in putative mapline.
        mapline.values <- unique( as.character(putative.positions) )
        
        # Get set of invalid values in mapline.
        invalid.values <- mapline.values[ ! ( isValidGenotype(mapline.values) |
            is.na(mapline.values) ) ]
        
        # Check mapline contains only valid genotypes and missing values.
        if ( length(invalid.values) > 0 ) {
            mapline.numbers <- mapline.values[ ! is.na( suppressWarnings(
                as.numeric(mapline.values) ) ) ]
            if ( length(mapline.numbers) == length(mapline.values) ) {
                stop("cross geno map positions must include centiMorgan units (e.g. '47 cM')")
            } else {
                stop("invalid genotype symbols - '", toString(invalid.values), "'")
            }
        }
    }
    
    # Get indices of initial rows.
    head.rows <- if (map.positions.found) {1:3} else {1:2}
    
    # Get offset of main body of genotype data.
    dat.offset = length(head.rows)
    
    # Get index of first and last data rows.
    first.data.row <- dat.offset + 1
    last.data.row <- nrow(from)
    
    # Trim any empty rows from the bottom.
    while ( allNA( from[last.data.row, ] ) ) {
        last.data.row <- last.data.row - 1
    }
    
    stopifnot( last.data.row >= first.data.row )
    
    # Get vector of data row indices.
    dat.rows <- first.data.row : last.data.row
    
    # Get matrix of genotype data.
    geno.matrix <- as.matrix(from[dat.rows, ])
    
    # Get set of symbols in genotype matrix.
    symbols <- unique( as.character( unlist(geno.matrix) ) )
    
    # Decompose symbols into different types.
    founder.symbols <- symbols[ isFounderGenotype(symbols) ]
    rawgeno.symbols <- symbols[ isRawGenotype(symbols) ]
    invalid.values <- symbols[ ! ( isValidGenotype(symbols) | is.na(symbols) ) ]
    
    # Check for invalid values.
    if ( length(invalid.values) > 0 ) {
        stop("invalid genotype symbols - '", toString(invalid.values), "'")
    }
    
    # Check for map positions that resemble raw genotypes.
    if ( ! map.positions.found && length(rawgeno.symbols) > 0 &&
         ! all( is.na(mapline.values) | mapline.values == '1' ) ) {
        stop("cross geno map positions must include centiMorgan units (e.g. '47 cM')")
    }
    
    # Get genotype symbols.
    if ( length(founder.symbols) > 0 && length(rawgeno.symbols) > 0 ) {
        stop("cross geno contains both raw and founder genotypes")
    } else if ( length(founder.symbols) > 0 ) {
        genotypes <- founder.symbols
    } else if ( length(rawgeno.symbols) > 0 ) {
        genotypes <- rawgeno.symbols
    } else {
        stop("cross geno has no genotype data")
    }
    
    # Get allele symbols from characters in genotype symbols.
    alleles <- unique( unlist( strsplit(genotypes, '') ) )
    
    # Convert genotype character matrix to a numeric matrix, with
    # genotype numbers assigned to corresponding genotype symbols.
    for ( i in getIndices(genotypes) ) {
        geno.matrix[ geno.matrix == genotypes[i] ] <- i
    }
    geno.matrix <- apply(geno.matrix, 2, as.numeric)
    
    # If rownames present, set vector of sample IDs from these..
    if ( hasRownames(from) ) {
        samples <- rownames(from)[dat.rows]
    } else { # ..otherwise set vector of sample indices.
        samples <- getIndices(dat.rows)
    }
    
    # Create map table from initial rows.
    map.table <- as.data.frame( t(from[head.rows, ]), stringsAsFactors=FALSE)
    colnames(map.table) <- const$maptable.colnames[head.rows]
    map.table <- setRownamesFromColumn(map.table, col.name='id')
    
    # If map positions found, coerce map table to a map object..
    if (map.positions.found) {
        geno.map <- as.map(map.table)
    } else { # ..otherwise create a dummy map, as in R/qtl.
        geno.map <- makeDummyMap(map.table$chr, rownames(map.table))
    }
    
    # Get individual locus IDs.
    locus.ids <- pullLocusIDs(geno.map)
    
    #Get sequences corresponding to individual map loci.
    locus.seqs <- pullLocusSeq(geno.map)
    
    # Get map sequences.
    geno.seqs <- unique(locus.seqs)
    
    if ( length(geno.seqs) < const$min.spm  ) {
        stop("cannot coerce data frame to cross geno - too few sequences (min=",
             const$min.spm, ")")
    }
    
    # Create CrossInfo object.
    cross.info <- methods::new('CrossInfo')
    cross.info <- setMarkers(cross.info, markers=locus.ids)
    cross.info <- setSequences(cross.info, geno.seqs)
    cross.info <- setMarkerSeqs(cross.info, sequences=locus.seqs)
    cross.info <- setAlleles(cross.info, alleles)
    cross.info <- setCrossType(cross.info, crosstype)
    cross.info <- setSamples(cross.info, samples)
    
    # Init cross genotype object.
    cross.geno <- list()
    
    for ( geno.seq in geno.seqs ) {
        
        # Get genotype data for this sequence.
        seq.dat <- geno.matrix[, locus.seqs == geno.seq ]
        
        # Get map info for this sequence.
        seq.map <- geno.map[[geno.seq]]
        class(seq.map) <- 'numeric'
        
        if ( length(seq.map) < const$min.lps ) {
            stop("cannot coerce data frame to cross geno - sequence has too few loci - '",
                geno.seq, "'")
        }
        
        # Assign geno data and map for this sequence.
        cross.geno[[geno.seq]] <- list(data=seq.dat, map=seq.map)
        class(cross.geno[[geno.seq]]) <- 'A' # NB: assumes no 'X' chromosomes.
    }
    
    attr(cross.geno, 'info') <- cross.info
    
    class(cross.geno) <- c(crosstype, 'geno', 'list')
    
    return(cross.geno)
}

# makeFounderGenoMatrix (S3) ---------------------------------------------------
#' Make founder genotype matrix from genotype data.
#' 
#' Given input genotype data for cross segregant samples and their founder 
#' strains, this function assigns a symbol to each locus according to the 
#' inferred founder allele. 
#' 
#' @param x Sample genotype data.
#' @param founder.geno Founder genotype data.
#'   
#' @return A genotype matrix, with genotypes encoded as integers and their
#' corresponding allele symbols in the attribute \code{'alleles'}.
#' 
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings QualityScaledDNAStringSet
#' @keywords internal
#' @rdname makeFounderGenoMatrix
makeFounderGenoMatrix <- function(x, founder.geno) {
    UseMethod('makeFounderGenoMatrix', x)
}

# makeFounderGenoMatrix.DNAStringSet -------------------------------------------
#' @rdname makeFounderGenoMatrix
makeFounderGenoMatrix.DNAStringSet <- function(x, founder.geno) {
    
    stopifnot( 'DNAStringSet' %in% class(founder.geno) || 
        'QualityScaledDNAStringSet' %in% class(founder.geno) )
    stopifnot( 'loci' %in% names(x@metadata) )
    stopifnot( 'mapframe' %in% class(x@metadata[['loci']]) )
    stopifnot( 'loci' %in% names(founder.geno@metadata) )
    stopifnot( 'mapframe' %in% class(founder.geno@metadata[['loci']]) )
    
    if ( ! hasRownames(x@metadata[['loci']]) ) {
        stop("cannot make founder genotype matrix - no sample locus IDs found")
    }
   
    if ( ! hasRownames(founder.geno@metadata[['loci']]) ) {
        stop("cannot make founder genotype matrix - no founder locus IDs found")
    } 
    
    # Get sample IDs.
    sample.ids <- names(x)
    
    dup.samples <- sample.ids[ duplicated(sample.ids) ]
    if ( length(dup.samples) > 0 ) {
        stop("duplicate sample IDs - '", toString(dup.samples), "'")
    }
    
    invalid.ids <- sample.ids[ ! isValidID(sample.ids) ]
    if ( length(invalid.ids) > 0 ) {
        stop("invalid sample IDs - '", toString(invalid.ids), "'")
    }
    
    # Get founder IDs.
    founder.ids <- names(founder.geno)
    
    if ( length(founder.ids) != 2 ) { 
        stop("unsupported number of founder genotypes - '", length(founder.ids), "'")
    } 
    
    dup.founders <- sample.ids[ duplicated(founder.ids) ]
    if ( length(dup.founders) > 0 ) {
        stop("duplicate founder IDs - '", toString(dup.founders), "'")
    }
    
    invalid.ids <- founder.ids[ ! isValidID(founder.ids) ]
    if ( length(invalid.ids) > 0 ) {
        stop("invalid founder IDs - '", toString(invalid.ids), "'")
    }
    
    clashing.ids <- intersect(sample.ids, founder.ids)
    if ( length(clashing.ids) > 0 ) {
        stop("sample/founder ID clash - '", toString(clashing.ids), "'")
    }   
    
    if ( length(founder.ids) > length(const$founder.allele.charset) ) {
        stop("number of founders (", length(founder.ids),
            ") exceeds number of available allele symbols (",
            length(const$founder.allele.charset), ")")
    }
    
    # Get genotype loci common to samples and founders.
    loc.list <- list(x@metadata[['loci']], founder.geno@metadata[['loci']])
    intersection <- do.call(intersectLoci, loc.list)
    
    # Keep only common loci.
    founder.geno <- subsetByLoci(founder.geno, intersection)
    x <- subsetByLoci(x, intersection)
    
    # Get common locus IDs.
    loc.ids <- rownames(x@metadata[['loci']])
    
    # Convert sample genotype data to matrix.
    x <- Biostrings::as.matrix(x)
    colnames(x) <- loc.ids
    
    # Convert founder genotype data to matrix.
    founder.geno <- Biostrings::as.matrix(founder.geno)
    colnames(founder.geno) <- loc.ids
    
    # Init founder genotype matrix.
    geno.matrix <- matrix(nrow=length(sample.ids), ncol=length(loc.ids), 
        dimnames=list(NULL, loc.ids))
    
    for ( col in getIndices(loc.ids) ) {
        
        # Init genotype numbers for this locus.
        geno.numbers <- rep(NA_integer_, length(sample.ids))
        
        # Get sample symbols and genotypes for this locus.
        sample.symbols <- x[, col]
        sample.genotypes <- unique(sample.symbols[ sample.symbols != '.' ])
        
        # Get founder symbols and genotypes for this locus.
        founder.symbols <- founder.geno[, col]
        founder.genotypes <- unique(founder.symbols[ founder.symbols != '.' ])
        
        # Assign locus genotypes from matching founder.
        if ( length(founder.genotypes) == 2 ) { # TODO: support polyallelic markers.
            if ( length(sample.genotypes) > 1 ) {
                if ( all( sample.genotypes %in% founder.genotypes ) ) {
                    for ( i in getIndices(founder.genotypes) ) {
                        geno.numbers[ sample.symbols == founder.genotypes[i] ] <- i
                    }
                }
            }
        }
        
        # Set genotype numbers for this locus.
        geno.matrix[, col] <- geno.numbers
    }
    
    # Remove null loci.
    geno.matrix <- geno.matrix[, ! apply(geno.matrix, 2, allNA) ]
    
    if ( ncol(geno.matrix) == 0 ) {
        stop("cannot make founder genotype matrix - no diallelic loci found")
    }
    
    # Set allele symbols for founders.
    attr(geno.matrix, 'alleles') <- const$founder.allele.charset[
        getIndices(founder.ids) ]
    
    return(geno.matrix)
}

# makeFounderGenoMatrix.QualityScaledDNAStringSet ------------------------------
#' @rdname makeFounderGenoMatrix
makeFounderGenoMatrix.QualityScaledDNAStringSet <- function(x, founder.geno) {
    return( makeFounderGenoMatrix.DNAStringSet(x, founder.geno) )
}

# makeFounderGenoMatrix (S4) ---------------------------------------------------
#' @rdname makeFounderGenoMatrix
setGeneric('makeFounderGenoMatrix', makeFounderGenoMatrix)

# DNAStringSet::makeFounderGenoMatrix ------------------------------------------
#' @rdname makeFounderGenoMatrix
setMethod('makeFounderGenoMatrix', signature='DNAStringSet', 
    definition=makeFounderGenoMatrix.DNAStringSet)

# QualityScaledDNAStringSet::makeFounderGenoMatrix -----------------------------
#' @rdname makeFounderGenoMatrix
setMethod('makeFounderGenoMatrix', signature='QualityScaledDNAStringSet', 
    definition=makeFounderGenoMatrix.DNAStringSet)

# makeRawGenoMatrix (S3) -------------------------------------------------------
#' Make raw genotype matrix from genotype data.
#' 
#' Given input sample genotype data, this function assigns an arbitrary symbol 
#' at each locus according to the observed raw SNP genotype. So for example, if
#' the SNVs at a given locus are 'A' and 'C', samples are assigned the genotypes 
#' '1' and '2', respectively. The mapping of SNV to genotype is performed
#' independently for each locus, so a given raw genotype does not have the same
#' meaning across loci.
#' 
#' @param x Sample genotype data.
#'   
#' @return A genotype matrix, with genotypes encoded as integers and their
#' corresponding allele symbols in the attribute \code{'alleles'}.
#' 
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings QualityScaledDNAStringSet
#' @keywords internal
#' @rdname makeRawGenoMatrix
makeRawGenoMatrix <- function(x) {
    UseMethod('makeRawGenoMatrix', x)
}

# makeRawGenoMatrix.DNAStringSet -----------------------------------------------
#' @rdname makeRawGenoMatrix
makeRawGenoMatrix.DNAStringSet <- function(x) {
    
    stopifnot( 'loci' %in% names(x@metadata) )
    stopifnot( 'mapframe' %in% class(x@metadata[['loci']]) )
    
    if ( ! hasRownames(x@metadata[['loci']]) ) {
        stop("cannot make raw genotype matrix - no locus IDs found")
    }

    # Get sample IDs.
    sample.ids <- names(x)
    
    dup.samples <- sample.ids[ duplicated(sample.ids) ]
    if ( length(dup.samples) > 0 ) {
        stop("duplicate sample IDs - '", toString(dup.samples), "'")
    }
    
    invalid.ids <- sample.ids[ ! isValidID(sample.ids) ]
    if ( length(invalid.ids) > 0 ) {
        stop("invalid sample IDs - '", toString(invalid.ids), "'")
    }
    
    # Get locus IDs.
    loc.ids <- rownames(x@metadata[['loci']])
    
    # Get number of loci, check consistent.
    num.loci <- unique( Biostrings::width(x) )
    stopifnot( length(num.loci) == 1 )
    stopifnot( length(num.loci) > 0 )
    stopifnot( length(loc.ids) == num.loci )
    
    # Convert sample genotype data to matrix.
    x <- Biostrings::as.matrix(x)
    colnames(x) <- loc.ids
    
    # Init raw genotype matrix.
    geno.matrix <- matrix(nrow=length(sample.ids), ncol=num.loci, 
        dimnames=list(NULL, loc.ids))
    
    for ( col in 1:num.loci ) {
        
        # Init genotype numbers for this locus.
        geno.numbers <- rep(NA_integer_, length(sample.ids))
        
        # Get sample symbols and alleles for this locus.
        sample.symbols <- x[, col]
        sample.genotypes <- unique( sample.symbols[ sample.symbols != '.' ] )
        
        # Skip loci with more raw genotypes than can be represented.
        if ( length(sample.genotypes) > length(const$raw.allele.charset) ) {
            next
        }
        
        # Assign locus genotypes in alphabetical order of symbol.
        if ( length(sample.genotypes) == 2 ) { # TODO: support polyallelic markers.
            for ( i in getIndices(sample.genotypes) ) {
                geno.numbers[ sample.symbols == sample.genotypes[i] ] <- i
            }
        }
        
        # Set genotype numbers for this locus.
        geno.matrix[, col] <- geno.numbers
    }
    
    # Remove null loci.
    geno.matrix <- geno.matrix[, ! apply(geno.matrix, 2, allNA) ]
    
    if ( ncol(geno.matrix) == 0 ) {
        stop("cannot make raw genotype matrix - no diallelic loci found")
    }
    
    # Set allele symbols.
    attr(geno.matrix, 'alleles') <- as.character( 1:max(geno.matrix, na.rm=TRUE) )
    
    return(geno.matrix)
}

# makeRawGenoMatrix.QualityScaledDNAStringSet ----------------------------------
#' @rdname makeRawGenoMatrix
makeRawGenoMatrix.QualityScaledDNAStringSet <- function(x) {
    return( makeRawGenoMatrix.DNAStringSet(x) )
}

# makeRawGenoMatrix (S4) -------------------------------------------------------
#' @rdname makeRawGenoMatrix
setGeneric('makeRawGenoMatrix', makeRawGenoMatrix)

# DNAStringSet::makeRawGenoMatrix ----------------------------------------------
#' @rdname makeRawGenoMatrix
setMethod('makeRawGenoMatrix', signature='DNAStringSet', 
    definition=makeRawGenoMatrix.DNAStringSet)

# QualityScaledDNAStringSet::makeRawGenoMatrix ---------------------------------
#' @rdname makeRawGenoMatrix
setMethod('makeRawGenoMatrix', signature='QualityScaledDNAStringSet', 
    definition=makeRawGenoMatrix.DNAStringSet)

# makeGeno (S3) ----------------------------------------------------------------
#' Make an \pkg{R/qtl} \code{cross} \code{geno} object.
#' 
#' @param sample.geno Sample genotype data.
#' @param founder.geno Founder genotype data.
#' 
#' @return A \code{geno} object.
#' 
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings QualityScaledDNAStringSet
#' @keywords internal
#' @rdname makeGeno
makeGeno <- function(sample.geno, founder.geno=NULL) {
    UseMethod('makeGeno', sample.geno)
}

# makeGeno.DNAStringSet --------------------------------------------------------
#' @rdname makeGeno
makeGeno.DNAStringSet <- function(sample.geno, founder.geno=NULL) {
    
    crosstype <- 'bc' # TODO: support other cross types
    
    # If founder genotypes available, get founder genotype matrix..
    if ( ! is.null(founder.geno) ) {
        geno.matrix <- makeFounderGenoMatrix(sample.geno, founder.geno)
    } else { # ..otherwise get raw genotype matrix.
        geno.matrix <- makeRawGenoMatrix(sample.geno)
    }
    
    # Subset sample genotype data to contain only loci in the genotype matrix.
    sample.geno <- subsetByLocusID(sample.geno, 
        function(loc.id) loc.id %in% colnames(geno.matrix))
    
    # Get locus info.
    locus.ids <- colnames(geno.matrix)
    
    # Make a map from sample genotype loci.
    # NB: resolves and sorts sequence IDs.
    geno.map <- makeMap(sample.geno)
    
    # Get allele symbols from genotype matrix.
    alleles <- attr(geno.matrix, 'alleles')
    attr(geno.matrix, 'alleles') <- NULL
    
    # Get sequences of map loci.
    locus.seqs <- pullLocusSeq(geno.map)
    
    # Get map sequences.
    geno.seqs <- unique(locus.seqs)
    
    if ( length(geno.seqs) < const$min.spm  ) {
        stop("cannot make cross geno - too few sequences (min=",
            const$min.spm, ")")
    }
    
    # Create CrossInfo object.
    cross.info <- methods::new('CrossInfo')
    cross.info <- setMarkers(cross.info, markers=locus.ids)
    cross.info <- setMarkerSeqs(cross.info, sequences=locus.seqs)
    cross.info <- setAlleles(cross.info, alleles)
    cross.info <- setCrossType(cross.info, crosstype)
    cross.info <- setSamples(cross.info, names(sample.geno))
    cross.info <- setSequences(cross.info, geno.seqs)
    
    # Init cross genotype object.
    cross.geno <- list()
    
    for ( geno.seq in geno.seqs ) {
   
        # Get genotype data for this sequence.
        seq.dat <- geno.matrix[, locus.seqs == geno.seq ]
        
        # Get map info for this sequence.
        seq.map <- geno.map[[geno.seq]]
        class(seq.map) <- 'numeric'
        
        if ( length(seq.map) < const$min.lps ) {
            stop("cannot make cross geno - sequence has too few loci - '",
                geno.seq, "'")
        }
        
        # Assign geno data and map for this sequence.
        cross.geno[[geno.seq]] <- list(data=seq.dat, map=seq.map)
        class(cross.geno[[geno.seq]]) <- 'A' # NB: assumes no 'X' chromosomes.
    }
    
    attr(cross.geno, 'info') <- cross.info
    
    class(cross.geno) <- c(crosstype, 'geno', 'list')
    
    return(cross.geno)
}

# makeGeno.QualityScaledDNAStringSet -------------------------------------------
#' @rdname makeGeno
makeGeno.QualityScaledDNAStringSet <- function(sample.geno, founder.geno=NULL) {
    return( makeGeno.DNAStringSet(sample.geno, founder.geno=founder.geno) )
}

# makeGeno (S4) ----------------------------------------------------------------
#' @rdname makeGeno
setGeneric('makeGeno', makeGeno)

# DNAStringSet::makeGeno -------------------------------------------------------
#' @rdname makeGeno
setMethod('makeGeno', signature='DNAStringSet', 
    definition=makeGeno.DNAStringSet)

# QualityScaledDNAStringSet::makeGeno ------------------------------------------
#' @rdname makeGeno
setMethod('makeGeno', signature='QualityScaledDNAStringSet', 
    definition=makeGeno.DNAStringSet)

# End of geno.R ################################################################