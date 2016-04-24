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
#' @param missing.value Missing data value. This can be any single character
#' that is not a possible phenotype or genotype value.
#' @param include.mapunit Include map unit information in map positions.
#' 
#' @return A \code{data.frame} corresponding to the input \code{geno} object.
#' 
#' @keywords internal
#' @method as.data.frame geno
#' @rdname as.data.frame.geno
as.data.frame.geno <- function(x, ..., chr=NULL, digits=NULL, missing.value='-',
    include.mapunit=TRUE) {
    
    stopifnot( isBOOL(include.mapunit) )
    
    map.unit <- 'cM'
    
    # Get CrossInfo object.
    cross.info <- attr(x, 'info')
    
    if ( is.null(cross.info) ) {
        stop("no CrossInfo found in cross geno object")
    }
    
    # Get relevant CrossInfo.
    alleles <- getAlleles(cross.info)
    locus.ids <- getMarkers(cross.info)
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
    
    geno.matrix[ is.na(geno.matrix) ] <- missing.value
    
    # Prepare map matrix.
    map.table <- insertColumn(map.table, col.index=1, 
        col.name='id', data=locus.ids)
    rownames(map.table) <- NULL
    map.matrix <- t(map.table)
    
    # Prepare genotype matrix.
    dimnames(geno.matrix) <- list(samples, NULL)
    
    # Bind map and genotype matrices into one data frame.
    geno.frame <- as.data.frame(rbind(map.matrix, geno.matrix), 
        stringsAsFactors=FALSE)
    
    return(geno.frame)
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
    locus.names <- make.names(locus.ids)
    
    # Ensure data has syntactially valid locus IDs. The original
    # locus IDs will be stored separately in the CrossInfo object.
    rownames(sample.geno@metadata[['loci']]) <- locus.names
    colnames(geno.matrix) <- locus.names
    
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
    
    # Create CrossInfo object.
    cross.info <- methods::new('CrossInfo')
    cross.info <- setMarkers(cross.info, markers=locus.ids)
    cross.info <- setMarkerSeqs(cross.info, sequences=locus.seqs)
    cross.info <- setAlleles(cross.info, alleles)
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
        
        # Assign geno data and map for this sequence.
        cross.geno[[geno.seq]] <- list(data=seq.dat, map=seq.map)
        class(cross.geno[[geno.seq]]) <- 'A' # NB: assumes no 'X' chromosomes.
    }
    
    attr(cross.geno, 'info') <- cross.info
    
    class(cross.geno) <- c('geno', 'list')
    
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