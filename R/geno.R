# Start of geno.R ##############################################################

# as.data.frame.geno -----------------------------------------------------------
#' Convert \code{geno} object to \code{data.frame}.
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
    genotypes <- getGenotypes(cross.info)
    locus.ids <- getMarkers(cross.info)
    locus.names <- make.names(locus.ids)
    
    if ( hasSampleIDs(cross.info) ) {
        samples <- getSamples(cross.info)
    } else {
        samples <- rep.int('', getNumSamples(cross.info))
    }
    
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
    geno.matrix <- decodeGenotypes(geno.matrix, genotypes)
    
    # Prepare map matrix.
    map.table <- insertColumn(map.table, col.index=1, 
        col.name='id', data=locus.ids)
    rownames(map.table) <- locus.names
    map.matrix <- t(map.table)
    rownames(map.matrix) <- NULL
    
    # Prepare genotype matrix.
    colnames(geno.matrix) <- locus.names
    
    # Bind map and genotype matrices into one data frame.
    geno.frame <- as.data.frame(rbind(map.matrix, geno.matrix), 
        stringsAsFactors=FALSE)
    
    # Insert sample ID column.
    geno.frame <- insertColumn(geno.frame, 1, col.name='id',
        data=c('id', '', '', samples))
    
    return(geno.frame)
}

# as.geno ----------------------------------------------------------------------
#' Convert to \code{geno} object.
#' 
#' @param from Object containing genotype data.
#' @param require.mapunit Require map unit information in map positions.
#' 
#' @return A \code{geno} object corresponding to the input object.
#' 
#' @keywords internal
#' @rdname as.geno
as.geno <- function(from, require.mapunit=TRUE) {
    UseMethod('as.geno', from)
}

# as.geno.data.frame -----------------------------------------------------------
#' @method as.geno data.frame
#' @rdname as.geno
as.geno.data.frame <- function(from, require.mapunit=TRUE) {
    
    stopifnot( isBOOL(require.mapunit) )
    
    crosstype <- 'haploid' # TODO: support other cross types
    
    if ( nrow(from) < 3 || ncol(from) < 2 ) {
        stop("invalid genotype data frame")
    }
    
    # Check for ID heading in first column.
    id.col <- which( tolower(from[1, ]) == 'id' )
    if ( length(id.col) == 0 || id.col[1] != 1 ) {
        stop("ID column not found in genotype data frame")
    }
    
    if ( from[2, id.col] != '' ) {
        stop("second row of ID column must be blank")
    }
    
    map.present <- ! is.na(from[3, id.col]) && from[3, id.col] == ''
    
    # Get indices of initial rows.
    head.rows <- if (map.present) {1:3} else {1:2}
    
    # Get offset of main body of genotype data.
    dat.offset = length(head.rows)
    
    # Get index of first and last data rows.
    first.data.row <- dat.offset + 1
    last.data.row <- nrow(from)
    
    stopifnot( last.data.row >= first.data.row )
    
    # Get vector of data row indices.
    dat.rows <- first.data.row : last.data.row
    
    if (map.present) {
        
        map.unit <- getMapUnitSuffix( as.character(from[3, -id.col]) )
        
        if ( ! is.na(map.unit) ) {
            
            if ( map.unit != 'cM' ) {
                stop("cross map positions must be in centiMorgans (e.g. '47 cM')")
            }
            
        } else {
            
            if (require.mapunit) {
                stop("cross map positions must include map units (e.g. '47 cM')")
            }
            
            map.unit <- 'cM'
        }
    }
    
    # Get matrix of genotype data.
    geno.matrix <- as.matrix(from[dat.rows, -id.col])
    colnames(geno.matrix) <- from[1, -id.col]
    
    # Get set of symbols in genotype matrix.
    genotypes <- sort( # NB: sort removes NA values.
        unique( as.character( unlist(geno.matrix) ) ) )
    
    # Make allele set from genotype set.
    alleles <- makeAlleleSet(genotypes) # NB: validates genotypes
    
    # Convert genotype character matrix to a numeric matrix, with
    # genotype numbers assigned to corresponding genotype symbols.
    geno.matrix <- encodeGenotypes(geno.matrix, genotypes)
    
    # Get sample IDs or indices.
    samples <- as.character(from[dat.rows, id.col])
    if ( allNA(samples) ) {
        samples <- seq_along(samples)
    } else if ( anyNA(samples) || any( samples == '' ) ) {
        stop("ID column is incomplete in genotype data frame")
    }
    
    # Create map table from initial rows.
    map.table <- as.data.frame( t(from[head.rows, -id.col]), stringsAsFactors=FALSE)
    colnames(map.table) <- const$maptable.colnames[head.rows]
    map.table <- setRownamesFromColumn(map.table, col.name='id')
    
    # If map positions found, convert map table to a map object..
    if (map.present) {
        geno.map <- as.map(map.table, map.unit=map.unit)
    } else { # ..otherwise create a placeholder map, as in R/qtl.
        geno.map <- makePlaceholderMap(map.table$chr, rownames(map.table))
    }
    
    # Get individual locus IDs.
    locus.ids <- pullLocusIDs(geno.map)
    
    #Get sequences corresponding to individual map loci.
    locus.seqs <- pullLocusSeq(geno.map)
    
    # Get map sequences.
    geno.seqs <- unique(locus.seqs)
    
    if ( length(geno.seqs) < const$min.spm  ) {
        stop("cannot convert data frame to cross geno - too few sequences (min=",
             const$min.spm, ")")
    }
    
    # Create CrossInfo object.
    cross.info <- methods::new('CrossInfo')
    cross.info <- setMarkers(cross.info, markers=locus.ids)
    cross.info <- setSequences(cross.info, geno.seqs)
    cross.info <- setMarkerSeqs(cross.info, sequences=locus.seqs)
    cross.info <- setAlleles(cross.info, alleles)
    cross.info <- setGenotypes(cross.info, genotypes)
    cross.info <- setCrosstype(cross.info, crosstype)
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
            stop("cannot convert data frame to cross geno - sequence has too few loci - '",
                geno.seq, "'")
        }
        
        # Assign geno data and map for this sequence.
        cross.geno[[geno.seq]] <- list(data=seq.dat, map=seq.map)
        class(cross.geno[[geno.seq]]) <- 'A' # NB: assumes no 'X' chromosomes.
    }
    
    attr(cross.geno, 'alleles') <- alleles
    
    attr(cross.geno, 'info') <- cross.info
    
    class(cross.geno) <- c('geno', 'list')
    
    return(cross.geno)
}

# makeFounderGenoMatrix --------------------------------------------------------
#' Make founder genotype matrix from genotype data.
#' 
#' Given input genotype data for cross segregant samples and their founder 
#' strains, this function assigns a symbol to each locus according to the 
#' inferred founder genotype.
#' 
#' @param x Raw sample genotype \code{array}.
#' @param y Raw founder genotype \code{array}.
#' @param alleles Founder allele symbols.
#'   
#' @return A founder genotype matrix, with genotypes encoded as integers and
#' their corresponding genotype symbols in the attribute \code{'genotypes'}.
#' 
#' @keywords internal
#' @rdname makeFounderGenoMatrix
makeFounderGenoMatrix <- function(x, y, alleles=NULL) {
    
    validateRawGenoArray(x)
    validateRawGenoArray(y)
    
    x.names <- dimnames(x)
    x.samples <- x.names[[1]]
    x.snps <- x.names[[2]]
    
    y.names <- dimnames(y)
    y.samples <- y.names[[1]]
    y.snps <- y.names[[2]]
    
    if ( ! is.null(alleles) ) {
        
        if ( length(alleles) != length(y.samples) ) {
            stop("founder allele length mismatch")
        }
        
        dup.alleles <- alleles[ duplicated(alleles) ]
        if ( length(dup.alleles) > 0 ) {
            stop("duplicate founder alleles - '", toString(dup.alleles), "'")
        }
        
        err.alleles <- alleles[ ! isFounderAllele(alleles) ]
        if ( length(err.alleles) > 0 ) {
            stop("invalid founder alleles - '", toString(err.alleles), "'")
        }
        
    } else {
        
        alleles <- const$founder.allele.charset[ seq_along(y.samples) ]
    }
    
    # Get genotype loci common to samples and founders.
    common.snps <- intersect(x.snps, y.snps)
    
    # Keep only common loci.
    x <- x[, common.snps, , drop=FALSE]
    y <- y[, common.snps, , drop=FALSE]
    
    # Extract raw genotype data.
    x.geno <- as.matrix(x[, , 'geno'])
    y.geno <- as.matrix(y[, , 'geno'])
    
    # TODO: use genotype error probabilities
    # x.prob <- as.matrix(x[, , 'prob'])
    # mode(x.prob) <- 'numeric'
    # y.prob <- as.matrix(y[, , 'prob'])
    # mode(y.prob) <- 'numeric'
    
    # Check no samples are marked as both founder and segregant.
    clashing.samples <- intersect(x.samples, y.samples)
    if ( length(clashing.samples) > 0 ) {
        stop("segregant/founder sample ID clash - '", toString(clashing.samples), "'")
    }   
    
    # TODO: handle more than two founder samples.
    if ( length(y.samples) != 2 ) {
        stop("unsupported number of founder genotypes - '", length(y.samples), "'")
    } 
    
    # Check upper limit of founder sample count.
    # TODO: uncomment this if not checking founder sample count elsewhere.
    # allele.count <- length(y.samples)
    # max.allele.count <- length(const$founder.allele.charset)
    # if ( allele.count > max.allele.count ) {
    #     stop("number of founders (", allele.count, ") exceeds number of "
    #         "available symbols (", max.allele.count), ")")
    # }
    
    # Verify that segregant genotypes are haploid.
    # TODO: handle other segregant ploidies.
    x.ploidy <- unique( nchar( unique( as.character(x.geno) ) ) )
    if ( x.ploidy != 1 ) {
        stop("unsupported segregant genotype ploidy - '", x.ploidy, "'")
    }
    
    # Verify that founder genotypes are haploid.
    # TODO: review issue of founder ploidy.
    y.ploidy <- unique( nchar( unique( as.character(y.geno) ) ) )
    if ( y.ploidy != 1 ) {
        stop("unsupported founder genotype ploidy - '", y.ploidy, "'")
    }
    
    # Init founder genotype matrix.
    geno.matrix <- matrix( NA_integer_, nrow=nrow(x.geno), ncol=ncol(x.geno),
        dimnames=dimnames(x.geno) )
    
    for ( j in getColIndices(geno.matrix) ) {
        
        # Get segregant symbols and genotypes for this locus.
        x.symbols <- x.geno[, j]
        x.levels <- unique(x.symbols)
        x.genotypes <- x.levels[ x.levels != const$missing.value ]
        
        # Get founder symbols and genotypes for this locus.
        y.genotypes <- y.symbols <- y.geno[, j]
        
        # Skip monomorphic/incomplete markers.
        # TODO: support polyallelic markers.
        if ( anyDuplicated(y.symbols) || const$missing.value %in% y.symbols ) {
            next
        }
        
        # Assign locus genotypes from matching founder.
        if ( length(x.genotypes) > 1 ) {
            if ( all( x.genotypes %in% y.genotypes ) ) {
                for ( i in seq_along(y.genotypes) ) {
                    geno.matrix[ x.symbols == y.genotypes[i], j ] <- i
                }
            }
        }
    }
    
    # Remove null loci.
    geno.matrix <- removeColsNA(geno.matrix)
    
    if ( ncol(geno.matrix) == 0 ) {
        stop("cannot make founder genotype matrix - no diallelic loci found")
    }
    
    # Set genotype symbols for founders.
    # TODO: handle other segregant ploidies.
    attr(geno.matrix, 'genotypes') <- alleles
    
    return(geno.matrix)
}

# makeEnumGenoMatrix -----------------------------------------------------------
#' Make enumerated genotype matrix from genotype data.
#' 
#' Given input sample genotype data, this function creates an enumerated genotype
#' matrix, in which each the genotype data at each locus are coded as numbers
#' in order of occurrence. So for example, if the first sample at a locus has
#' raw genotype \code{'A'}, and the second sample has raw genotype \code{'C'},
#' these will be assigned genotypes \code{'1'} and \code{'2'}, respectively.
#' The enumeration of genotypes is performed independently for each locus,
#' so a given enumerated genotype does not have the same meaning across loci.
#' 
#' @param x Raw sample genotype \code{array}.
#'   
#' @return An enumerated genotype matrix, with genotypes encoded as integers and
#' their corresponding genotype symbols in the attribute \code{'genotypes'}.
#' 
#' @keywords internal
#' @rdname makeEnumGenoMatrix
makeEnumGenoMatrix <- function(x) {
    
    validateRawGenoArray(x)
    
    # Extract raw genotype data as matrix.
    x.geno <- as.matrix(x[, , 'geno'])
    
    # Init enumerated genotype matrix.
    geno.matrix <- matrix(NA_integer_, nrow=nrow(x.geno), ncol=ncol(x.geno),
        dimnames=dimnames(x.geno) )
    
    for ( i in getColIndices(geno.matrix) ) {
        
        # Get sample symbols and genotypes for this locus.
        x.symbols <- x.geno[, i]
        x.levels <- unique(x.symbols)
        x.genotypes <- x.levels[ x.levels != const$missing.value ]
        
        # Skip loci with more genotypes than can be represented.
        if ( length(x.genotypes) > length(const$enum.geno.charset) ) {
            next
        }
        
        # Set genotype numbers for this locus.
        geno.matrix[, i] <- match(x.symbols, x.genotypes)
    }
    
    # Remove null loci.
    geno.matrix <- removeColsNA(geno.matrix)
    
    if ( ncol(geno.matrix) == 0 ) {
        stop("cannot make enumerated genotype matrix - no diallelic loci found")
    }
    
    # Set genotype symbols.
    attr(geno.matrix, 'genotypes') <- as.character( 1:max(geno.matrix, na.rm=TRUE) )
    
    return(geno.matrix)
}

# makeGeno ---------------------------------------------------------------------
#' Make an \pkg{R/qtl} \code{cross} \code{geno} object.
#' 
#' @param x Raw sample genotype \code{array}.
#' @param y Raw founder genotype \code{array}.
#' @param alleles Founder allele symbols.
#' 
#' @return A \code{geno} object.
#' 
#' @keywords internal
#' @rdname makeGeno
makeGeno <- function(x, y=NULL, alleles=NULL) {
    
    crosstype <- 'haploid' # TODO: support other cross types
    
    # If founder genotypes available, get founder genotype matrix..
    if ( ! is.null(y) ) {
        
        geno.matrix <- makeFounderGenoMatrix(x, y, alleles=alleles)
        
    } else { # ..otherwise get enumerated genotype matrix.
        
        if ( ! is.null(alleles) ) {
            stop("founder alleles specified without founders")
        }
        
        geno.matrix <- makeEnumGenoMatrix(x)
    }
    
    # Get genotype symbols from genotype matrix.
    genotypes <- attr(geno.matrix, 'genotypes')
    attr(geno.matrix, 'genotypes') <- NULL
    
    # Make allele set from genotype set.
    alleles <- makeAlleleSet(genotypes) # NB: validates genotypes
    
    # Get sample and SNP IDs from genotype matrix.
    sample.ids <- rownames(geno.matrix)
    snp.ids <- colnames(geno.matrix)
    dimnames(geno.matrix) <- NULL
    
    # Make a physical map from sample genotype loci.
    # NB: resolves and sorts sequence IDs.
    geno.map <- makeMapFromDefaultMarkerIDs(snp.ids)
    
    # Convert physical map to genetic map.
    geno.map <- setMapUnit(geno.map, 'cM')
    
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
    cross.info <- setMarkers(cross.info, markers=snp.ids)
    cross.info <- setMarkerSeqs(cross.info, sequences=locus.seqs)
    cross.info <- setAlleles(cross.info, alleles)
    cross.info <- setGenotypes(cross.info, genotypes)
    cross.info <- setCrosstype(cross.info, crosstype)
    cross.info <- setSamples(cross.info, sample.ids)
    cross.info <- setSequences(cross.info, geno.seqs)
    
    # Init cross genotype object.
    cross.geno <- list()
    
    for ( geno.seq in geno.seqs ) {
        
        # Get genotype data for this sequence.
        seq.dat <- geno.matrix[, locus.seqs == geno.seq, drop=FALSE]
        
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
    
    attr(cross.geno, 'alleles') <- alleles
    
    attr(cross.geno, 'info') <- cross.info
    
    class(cross.geno) <- c('geno', 'list')
    
    return(cross.geno)
}

# validateRawGenoArray ---------------------------------------------------------
#' Validate \code{array} as containing raw genotype data.
#' 
#' @param x Raw genotype \code{array}.
#' 
#' @return \code{TRUE} if object is an \code{array} containing raw
#' genotype data in expected form; otherwise, raises first error.
#' 
#' @keywords internal
#' @rdname validateRawGenoArray
validateRawGenoArray <- function(x) {
    
    stopifnot( is.array(x) )
    stopifnot( all( dim(x) > 0 ) )
    
    x.names <- dimnames(x)
    
    if ( length(x.names) != 3 ) {
        stop("raw genotype array has incorrect number of dimensions - '",
             length(x.names), "'")
    }
    
    sample.ids <- x.names[[1]]
    snp.ids <- x.names[[2]]
    slices <- x.names[[3]]
    
    if ( ! length(slices) %in% 1:2 || slices[1] != 'geno' ||
        ( length(slices) == 2 && slices[2] != 'prob' ) ) {
        stop("raw genotype array has invalid slice names")
    }
    
    if ( is.null(sample.ids) ) {
        stop("no sample IDs found in raw genotype array")
    }
    
    dup.samples <- sample.ids[ duplicated(sample.ids) ]
    if ( length(dup.samples) > 0 ) {
        stop("duplicate sample IDs - '", toString(dup.samples), "'")
    }
    
    err.samples <- sample.ids[ ! isValidID(sample.ids) ]
    if ( length(err.samples) > 0 ) {
        stop("invalid sample IDs - '", toString(err.samples), "'")
    }
    
    if ( is.null(snp.ids) ) {
        stop("no SNP IDs found in raw genotype array")
    }
    
    dup.snps <- snp.ids[ duplicated(snp.ids) ]
    if ( length(dup.snps) > 0 ) {
        stop("duplicate SNP IDs - '", toString(dup.snps), "'")
    }
    
    err.snps <- snp.ids[ ! isDefaultMarkerID(snp.ids) ]
    if ( length(err.snps) > 0 ) {
        stop("invalid SNP IDs - '", toString(err.snps), "'")
    }
    
    # Get symbols from genotype matrix.
    g.symbols <- unique( as.character(x[, , 'geno']) )
    
    # Get genotype symbols.
    genotypes <- g.symbols[ g.symbols != const$missing.value ]
    
    ploidies <- unique( nchar(genotypes) )
    if ( length(ploidies) > 1 ) {
        stop("mixed ploidy in raw genotype array")
    }
    
    # Get characters in genotype symbols.
    a.symbols <- unique( unlist( strsplit(genotypes, '') ) )
    
    # Get allele symbols.
    alleles <- a.symbols[ a.symbols != const$missing.value ]
    
    unknown <- alleles[ ! isRawAllele(alleles) ]
    if ( length(unknown) > 0 ) {
        stop("raw genotype array data contains unknown alleles - '",
            toString(unknown), "'")
    }
    
    if ( 'prob' %in% slices ) {
        prob <- as.numeric(x[, , 'prob'])
        if ( ! all( isProbability(prob) ) ) {
            stop("raw genotype array data contains invalid error probabilities")
        }
    }
}

# End of geno.R ################################################################
