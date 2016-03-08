# Start of geno.R ##############################################################

# TODO: function makeFounderGenoMatrix(); must use common loci.

# makeRawGenoMatrix (S3) -------------------------------------------------------
#' Make raw genotype matrix from genotype data.
#' 
#' Given input sample genotype data, this function assigns an arbitrary symbol 
#' at each locus according to the observed raw SNP genotype. So for example, if
#' the SNVs at a given locus are 'A' and 'C', samples are assigned the genotypes 
#' 'X1' and 'X2', respectively. The mapping of SNV to genotype is performed 
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
    
    if ( ! hasRownames(x) ) {
        stop("cannot make raw genotype matrix - no locus IDs found")
    }

    loc.ids <- rownames(x@metadata[['loci']])
    
    sample.ids <- names(x)
    
    dup.samples <- sample.ids[ duplicated(sample.ids) ]
    if ( length(dup.samples) > 0 ) {
        stop("duplicate sample IDs - '", toString(dup.samples), "'")
    }
    
    invalid.ids <- sample.ids[ ! isValidID(sample.ids) ]
    if ( length(invalid.ids) > 0 ) {
        stop("invalid sample IDs - '", toString(invalid.ids), "'")
    }
    
    num.samples <- length(x)
    num.loci <- unique( Biostrings::width(x) )
    stopifnot( length(num.loci) == 1 )
    stopifnot( length(loc.ids) == num.loci )
    
    orig.dat <- Biostrings::as.matrix(x)
    colnames(orig.dat) <- loc.ids
    
    geno.matrix <- matrix(nrow=num.samples, ncol=num.loci, 
        dimnames=list(sample.ids, loc.ids))
    
    num.alleles <- 0
    
    for ( col in 1:num.loci ) {
        
        geno.symbols <- orig.dat[, col]
        
        geno.numbers <- rep(NA_integer_, num.samples)
        
        locus.alleles <- sort( unique( geno.symbols[ geno.symbols != '.' ] ) )
        
        if ( length(locus.alleles) > num.alleles ) {
            num.alleles <- length(locus.alleles)
        }
        
        for ( i in 1:length(locus.alleles) ) {
            geno.numbers[ geno.symbols == locus.alleles[i] ] <- i
        }
        
        geno.matrix[, col] <- geno.numbers
    }
        
    if ( num.alleles != 2 ) { 
        stop("unsupported number of raw genotypes - '", num.alleles, "'")
    } 
        
    attr(geno.matrix, 'alleles') <- make.names(1:num.alleles)

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
#' @return An \pkg{R/qtl} \code{cross} \code{geno} object.
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
    
    geno.map <- makeMap(sample.geno)
    
    if ( ! is.null(founder.geno) ) {
        # TODO: geno.matrix <- makeFounderGenoMatrix(sample.geno, founder.geno=founder.geno)
        stop("cannot make founder genotype data (yet)")
    } else {
        geno.matrix <- makeRawGenoMatrix(sample.geno)
    }
    
    cross.geno <- list()
    
    attr(cross.geno, 'alleles') <- attr(geno.matrix, 'alleles')
    attr(geno.matrix, 'alleles') <- NULL
    
    locus.seqs <- pullLocusSeq(sample.geno@metadata[['loci']])
    
    geno.seqs <- unique(locus.seqs)
    
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

# makeGenoTable (S3) -----------------------------------------------------------
#' Make genotype table from genotype data.
#' 
#' @param x An \pkg{R/qtl} \code{cross} \code{geno} object.
#' @param chr Vector of sequences for which genotype data should be included in 
#' the genotype table. If none are specified, a genotype table is returned for 
#' all available sequences.
#' @param digits If specified, round genetic map positions to the specified 
#' number of digits.
#' @param include.mapunit Include map unit information in map positions.
#' 
#' @return A \code{data.frame} containing a genotype table.
#' 
#' @keywords internal
#' @rdname makeGenoTable
makeGenoTable <- function(x, chr=NULL, digits=NULL, include.mapunit=TRUE) {
    UseMethod('makeGenoTable', x)
}

# makeGenoTable.list -----------------------------------------------------------
#' @rdname makeGenoTable
makeGenoTable.list <- function(x, chr=NULL, digits=NULL, include.mapunit=TRUE) {
    
    # TODO: validateGeno(x)
    
    stopifnot( isBOOL(include.mapunit) )
    
    map.unit <- 'cM'
    
    if ( ! 'alleles' %in% names( attributes(x) ) ) {
        stop("no alleles found in cross geno object")
    }
    
    # Get alleles from geno object.
    alleles <- attr(x, 'alleles')

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
    for ( i in 1:length(alleles) ) {
        geno.matrix[ geno.matrix == i ] <- alleles[i]
    }
    
    # Prepare map.
    map.table <- setColumnFromRownames(map.table)
    map.matrix <- t(map.table)
    rownames(map.matrix) <- NULL
    map.matrix <- insertColumn(map.matrix, 1, data=c('id', '', ''))
    
    # Prepare genotype matrix.
    geno.matrix <- setColumnFromRownames(geno.matrix)
    colnames(geno.matrix) <- NULL
    
    # Bind map and genotype matrices into one table.
    geno.table <- as.data.frame(rbind(map.matrix, geno.matrix), 
        stringsAsFactors=FALSE)
    
    return(geno.table)
}

# End of geno.R ################################################################