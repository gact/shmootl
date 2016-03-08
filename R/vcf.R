# Start of vcf.R ###############################################################

# getSamplesVCF ----------------------------------------------------------------
#' Get sample IDs from VCF file.
#' 
#' @param file A VCF file path.
#' 
#' @return Vector of the sample IDs in the given VCF file.
#'   
#' @keywords internal
#' @rdname getSamplesVCF
getSamplesVCF <- function(file) {
    stopifnot( isSingleString(file) )
    stopifnot( all( file.exists(file) ) )
    header <- VariantAnnotation::scanVcfHeader(file)
    return( VariantAnnotation::samples(header) )
}

# hasGenoQualVCF ---------------------------------------------------------------
#' Test if VCF file contains GQ scores.
#' 
#' @param file A VCF file path.
#' 
#' @return TRUE if the specified VCF file contains genotype quality (GQ) scores;
#' FALSE otherwise.
#'   
#' @keywords internal
#' @rdname hasGenoQualVCF
hasGenoQualVCF <- function(file) {
    stopifnot( isSingleString(file) )
    stopifnot( all( file.exists(file) ) )
    header <- VariantAnnotation::scanVcfHeader(file)
    return( 'GQ' %in% rownames( VariantAnnotation::geno(header) ) )
}

# readSnpsVCF ------------------------------------------------------------------
#' Read SNP genotypes from VCF files.
#' 
#' This function reads SNP genotype data from one or more VCF files, and returns 
#' these as a sequence of raw SNP genotypes - effectively variant base calls - 
#' for each sample. If all relevant input VCF files contain genotype quality (GQ) 
#' scores for the samples of interest, these are used - along with variant 
#' quality scores - to calculate a quality score for each SNP genotype, and 
#' the result is returned as a \pkg{Biostrings} \code{QualityScaledDNAStringSet} 
#' object containing raw SNP genotypes and their corresponding Phred quality 
#' scores, encoded as an ASCII string. Otherwise, SNP genotypes are returned 
#' as a \pkg{Biostrings} \code{DNAStringSet} object containing only raw SNP 
#' genotype values. In any case, the returned object will have an attribute 
#' \code{'loci'}: a \code{mapframe} of marker loci corresponding to the 
#' genotyped variants.
#' 
#' @param ... Input VCF file paths.
#' @param samples Vector of samples for which SNP genotypes should be obtained. If
#' not specified, genotypes are returned for all available samples.
#' @param require.complete Remove variants that are not completely genotyped 
#' with respect to the given samples.
#' @param require.polymorphic Remove variants that do not have at least two 
#' different genotype calls among the given samples. 
#' 
#' @return A \pkg{Biostrings} \code{QualityScaledDNAStringSet} object if 
#' genotype quality scores are available; otherwise a \pkg{Biostrings} 
#' \code{DNAStringSet} object containing only raw SNP genotype values. 
#'   
#' @importFrom BiocGenerics lengths
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges ranges
#' @keywords internal
#' @rdname readSnpsVCF
readSnpsVCF <- function(..., samples=NULL, require.complete=FALSE, 
    require.polymorphic=FALSE) {
    
    infiles <- list(...)
    stopifnot( length(infiles) > 0 )
    stopifnot( isBOOL(require.complete) )
    stopifnot( isBOOL(require.polymorphic) )
    
    # Set NA string to the 'not available' 
    # character of `Biostrings::DNA_ALPHABET`.
    na.string <- '.'
    
    # Set reference genome.
    # TODO: check contig info against reference genome.
    genome <- genomeOpt()
    
    # Get samples in each input VCF file.
    sample.list <- lapply(infiles, getSamplesVCF)
    
    # If samples specified, filter sample list, keep only those specified. 
    if ( ! is.null(samples) ) {
        
        stopifnot( is.character(samples) )
        stopifnot( length(samples) > 0 )
        
        for ( i in 1:length(infiles) ) {
            sample.list[[i]] <- sample.list[[i]][ sample.list[[i]] %in% samples ]
        }
        
        unfound <- samples[ ! samples %in% unlist(sample.list) ]
        if ( length(unfound) > 0 ) {
            stop("samples not found - '", toString(unfound), "'")
        }
    }
    
    samples <- unlist(sample.list)

    if ( length(samples) == 0 ) {
        stop("no samples found")
    }
        
    dup.samples <- samples[ duplicated(samples) ]
    if ( length(dup.samples) > 0 ) {
        stop("duplicate samples - '", toString(dup.samples), "'")
    }
    
    invalid.ids <- samples[ ! isValidID(samples) ]
    if ( length(invalid.ids) > 0 ) {
        stop("invalid sample IDs - '", toString(invalid.ids), "'")
    }
    
    # Keep only files relevant for the given samples.
    relevant <- lengths(sample.list) > 0
    sample.list <- sample.list[relevant]
    infiles <- infiles[relevant]
    
    # Get mask indicating which relevant files contain genotype qualities.
    gq.mask <- sapply(infiles, hasGenoQualVCF)
    
    # If all relevant files contain genotype qualities, 
    # calculate SNP genotype quality scores..
    if ( all(gq.mask) ) {
        reading.qualities <- TRUE
    } else { # ..otherwise return only SNP genotypes.
        reading.qualities <- FALSE
        if ( any(gq.mask) ) {
            warning("some VCF files lack GQ scores, reading genotypes only")
        }
    }

    snps <- vector('list', length(infiles))
    
    # Read SNP genotypes from each input VCF file.
    for ( i in 1:length(infiles) ) {
        
        file.samples <- sample.list[[i]]
        infile <- infiles[[i]]
        
        # Set genotype fields to fetch.
        if (reading.qualities) {
            geno.fields <- c('GT', 'GQ')
        } else {
            geno.fields <- 'GT'
        }
        
        # Set VCF parameters.
        param <- VariantAnnotation::ScanVcfParam(fixed=c('ALT', 'QUAL'), 
            geno=geno.fields, samples=file.samples)
        
        # Read VCF.
        vcf <-  VariantAnnotation::readVcf(infile, genome=genome, param=param)
        stopifnot( length(vcf) > 0 )
        
        # Get variant data from VCF.
        var.seqs <- normSeq( as.character(                     # reference sequence
            GenomeInfoDb::seqnames(vcf) ) )
        var.pos <- BiocGenerics::start( IRanges::ranges(vcf) ) # position  
        var.ref <- as.character( VariantAnnotation::ref(vcf) ) # reference allele
        var.alt <- lapply( IRanges::CharacterList(             # alternative alleles
            VariantAnnotation::alt(vcf) ), as.character)
        var.qual <- VariantAnnotation::qual(vcf)               # variant quality
        var.geno <- VariantAnnotation::geno(vcf)               # genotype info 
        var.gt <- var.geno$GT                                  # genotype calls
        
        invalid <- unique( var.gt[ ! grepl(const$pattern$vcf.haploid.geno, var.gt) ] )
        if ( length(invalid) > 0 ) {
            stop("invalid haploid VCF genotypes - '", toString(invalid), "'")
        }
        
        # Convert genotype allele codes to integers, VCF missing values to NA.
        var.gt <- suppressWarnings( apply(var.gt, 2, as.integer) )
        
        # Set indices of SNP variants.
        indices <- which( VariantAnnotation::isSNV(vcf) )
        
        # Get variant allele values.
        allele.values <- lapply(indices, function(i) 
            c(var.ref[i], var.alt[[i]]) )   
        
        # Get variant allele codes, converting 0-offset allele code  
        # used in VCF format to the 1-offset allele code used in R.
        allele.codes <- sapply(indices, function(i) var.gt[i, ] + 1)
        
        # Get base call for each sample genotype.
        base.matrix <- sapply(1:length(indices), function(i) 
            allele.values[[i]][ allele.codes[, i] ])
        
        # Replace NA values.
        base.matrix[ is.na(base.matrix) ] <- na.string
        
        # Get base call string for each sample.
        bases <- sapply(1:length(file.samples), function(s) 
            paste0(base.matrix[s, ], collapse='') )
        
        # If genotype qualities available, create  
        # object with quality-scaled genotypes..
        if (reading.qualities) {
            
            # Set vector of variant error probabilities.
            variant.error.probs <- 10 ^ ( -0.1 * var.qual[indices] )
            
            # Set matrix of genotype error probabilities.
            genotype.error.probs <- 10 ^ ( -0.1 * var.geno$GQ[indices, ] )
            
            # Calculate a Phred quality for each sample genotype  
            # from combined variant/genotype quality scores.
            qual.matrix <- sapply(1:length(indices), function(i) 
                -10 * log10( variant.error.probs[i] + genotype.error.probs[i, ] ) )
            
            # Replace NA values with minimum quality score.
            qual.matrix[ is.na(qual.matrix) ] <- const$qual$phred$range[1]
            
            # Clamp quality scores within Phred range, add Phred offset.
            qual.matrix <- sapply(1:length(file.samples), function(s) 
                clamp(qual.matrix[s, ], const$qual$phred$range) + 
                const$qual$phred$offset )
            
            # Convert Phred quality scores to ASCII characters.
            quals <- sapply(1:length(file.samples), function(s) 
                intToUtf8( round(qual.matrix[, s]) ) )
            
            snps[[i]] <- Biostrings::QualityScaledDNAStringSet( 
                Biostrings::DNAStringSet(bases), 
                Biostrings::PhredQuality(quals) )
            names(snps[[i]]@quality) <- file.samples
            names(snps[[i]]) <- file.samples
            
        } else { # ..otherwise create object with only genotypes.
            
            snps[[i]] <- Biostrings::DNAStringSet(bases)
            names(snps[[i]]) <- file.samples
        }
        
        # Add marker metadata.
        loc <- mapframe(chr=var.seqs[indices], 
            pos=var.pos[indices], map.unit='bp')
        rownames(loc) <- makeDefaultMarkerIDs(loc)
        snps[[i]]@metadata[['loci']] <- loc
    }
    
    # If multiple relevant files, combine SNP genotypes..
    if ( length(infiles) > 1 ) {
        
        # Get common SNP genotype loci.
        loc.list <- lapply(snps, function(x) x@metadata[['loci']])
        intersection <- do.call(intersectLoci, loc.list)
        
        if ( nrow(intersection) == 0 ) {
            stop("no common loci found in VCF files - ", 
                toString(infiles), "'")
        }
        
        # Keep only common loci.
        for ( i in 1:length(snps) ) {
            snps[[i]] <- subsetByLoci(snps[[i]], intersection)
        }
        
        # Combine SNP genotypes.
        result <- snps[[1]]
        if ( length(snps) > 1 ) {
            for ( i in 2:length(snps) ) {
                result <- append(result, snps[[i]])
            }
        }
        
    } else { # ..otherwise take single set of SNP genotypes.
        
        result <- snps[[1]]
    }
    
    # Get number of SNP variants.
    num.snps <- unique( BiocGenerics::lengths(result) )
    
    stopifnot( length(num.snps) == 1 )
    stopifnot( num.snps > 0 )
    
    # If filters specified, filter SNP variants.
    if ( require.complete || require.polymorphic ) {
        
        mask <- rep(TRUE, num.snps)
        
        for ( i in 1:num.snps ) {
            
            gt <- sapply(result, function(x) as.character(x[i]))
            
            if (require.complete) {
                if ( any( gt == na.string ) ) {
                    mask[i] <- FALSE
                    next
                }
            }
            
            if (require.polymorphic) {
                if ( length( unique( gt[ gt != na.string ]) ) > 1 ) {
                    mask[i] <- FALSE
                    next
                }
            }
        }
        
        # Apply filter mask.
        result@metadata[['loci']] <- result@metadata[['loci']][mask, ]
        for ( i in 1:length(result) ) {
            result[[i]] <- result[[i]][mask] 
            if ( 'QualityScaledDNAStringSet' %in% class(result) ) {
                result@quality[[i]] <- result@quality[[i]][mask]
            }
        }
    }

    return(result)
}

# readGenoVCF ------------------------------------------------------------------
#' Read genotype data from a VCF file.
#' 
#' This function reads SNP genotype data from one or more VCF files, and returns 
#' these as an \pkg{R/qtl} \code{cross} \code{geno} object. If cross samples are 
#' not specified, all input samples are used. 
#' 
#' If no founder samples are specified, this function assigns an arbitrary symbol 
#' at each locus according to the observed raw SNP genotype. So for example, if
#' the SNVs at a given locus are 'A' and 'C', samples are assigned the genotypes 
#' 'X1' and 'X2', respectively. The mapping of SNV to genotype is performed 
#' independently for each locus, so a given raw genotype does not have the same
#' meaning across loci.
#' 
#' If founder samples are specified, this function assigns a genotype symbol 
#' according to the inferred founder for each genotype at each locus. In this
#' case, a given genotype represents the same founder strain across markers.
#' 
#' @param ... Input VCF file paths.
#' @param samples Cross sample IDs.
#' @param founders Founder sample IDs.
#' 
#' @return An \pkg{R/qtl} \code{cross} \code{geno} object.
#' 
#' @export
#' @rdname readGenoVCF
readGenoVCF <- function(..., samples, founders=NULL) {
    
    # TODO: optimise.
    
    sample.geno <- readSnpsVCF(..., samples=samples, 
        require.polymorphic=TRUE)
    
    if ( ! is.null(founders) ) {
        
        clashing.samples <- intersect(samples, founders)
        if ( length(clashing.samples) > 0 ) {
            stop("clashing sample IDs - ", toString(clashing.samples), "'")
        }
        
        founder.geno <- readSnpsVCF(..., samples=founders, 
            require.complete=TRUE, require.polymorphic=TRUE)
        
    } else {
        
        founder.geno <- NULL
    }
    
    geno <- makeGeno(sample.geno, founder.geno=founder.geno)
    
    return(geno)
}

# End of vcf.R #################################################################