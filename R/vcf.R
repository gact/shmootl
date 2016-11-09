# Start of vcf.R ###############################################################

# getSamplesVCF ----------------------------------------------------------------
#' Get sample IDs from VCF file.
#' 
#' @param file A VCF file path.
#' 
#' @return Vector of the sample IDs in the given VCF file.
#'   
#' @export
#' @family VCF functions
#' @rdname getSamplesVCF
getSamplesVCF <- function(file) {
    stopifnot( isSingleString(file) )
    stopifnot( file.exists(file) )
    header <- VariantAnnotation::scanVcfHeader(file)
    return( VariantAnnotation::samples(header) )
}

# hasGenoQualVCF ---------------------------------------------------------------
#' Test if VCF file contains GQ scores.
#' 
#' @param file A VCF file path.
#' 
#' @return \code{TRUE} if the specified VCF file contains
#' genotype quality (GQ) scores; \code{FALSE} otherwise.
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
#' Read raw SNP genotypes from VCF files.
#' 
#' This function reads SNP genotype data from one or more VCF files, and returns
#' these as a sequence of raw SNP genotypes - effectively variant base calls -
#' for each sample. If all relevant input VCF files contain genotype quality
#' (GQ) scores for the samples of interest, these are used - along with variant
#' quality scores - to calculate an error probability for each SNP genotype, and
#' the result is returned as an \code{array} with two slices - \code{'geno'} and
#' \code{'prob'} - containing raw SNP genotypes and their corresponding error
#' probabilities, respectively. Otherwise, SNP genotypes are returned as an
#' \code{array} with one slice, \code{'geno'}, which contains raw SNP genotype
#' values. In any case, the rows of each slice contain data for a given sample,
#' while the columns of each slice contain data for a given SNP locus.
#' 
#' @param ... Input VCF file paths.
#' @param samples Vector of samples for which SNP genotypes should be obtained.
#' If not specified, genotypes are returned for all available samples.
#' @param require.all Remove variants that are not completely genotyped
#' with respect to the given samples.
#' @param require.any Remove variants that do not have at least one genotype
#' call among the given samples.
#' @param require.polymorphic Remove variants that do not have at least two
#' different genotype calls among the given samples.
#' 
#' @return An \code{array} object containing raw genotype data, and if available,
#' genotype error probabilities.
#'   
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges CharacterList
#' @importFrom IRanges ranges
#' @keywords internal
#' @rdname readSnpsVCF
readSnpsVCF <- function(..., samples=NULL, require.all=FALSE, require.any=FALSE,
    require.polymorphic=FALSE) {
    
    infiles <- unlist( list(...) )
    stopifnot( length(infiles) > 0 )
    stopifnot( isBOOL(require.all) )
    stopifnot( isBOOL(require.polymorphic) )

    # Set reference genome.
    # TODO: check contig info against reference genome.
    genome <- genomeOpt()
    
    # Get samples in each input VCF file.
    sample.list <- lapply(infiles, getSamplesVCF)
    
    # If samples specified, filter sample list, keep only those specified. 
    if ( ! is.null(samples) ) {
        
        stopifnot( is.character(samples) )
        stopifnot( length(samples) > 0 )
        
        for ( i in seq_along(infiles) ) {
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
    gq.mask <- sapply(infiles, hasGenoQualVCF, USE.NAMES=FALSE)
    
    # If all relevant files contain genotype qualities, 
    # calculate SNP genotype quality scores..
    if ( all(gq.mask) ) {
        reading.qualities <- TRUE
        slices <- c('geno', 'prob')
    } else { # ..otherwise return only SNP genotypes.
        reading.qualities <- FALSE
        slices <- c('geno')
        if ( any(gq.mask) ) {
            warning("some VCF files lack GQ scores, reading genotypes only")
        }
    }
    
    # Init list of SNP data.
    snps <- vector('list', length(infiles))
    
    # Init list for mapping genotype strings to allele indices.
    geno.memo <- list()
    
    # Read SNP genotypes from each input VCF file.
    for ( i in seq_along(infiles) ) {
        
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
        var.seqs <- as.character( GenomeInfoDb::seqnames(vcf) ) # reference sequence
        var.pos <- BiocGenerics::start( IRanges::ranges(vcf) )  # position
        var.ref <- as.character( VariantAnnotation::ref(vcf) )  # reference allele
        var.alt <- lapply( IRanges::CharacterList(              # alternative alleles
            VariantAnnotation::alt(vcf) ), as.character)
        var.qual <- VariantAnnotation::qual(vcf)                # variant quality
        var.geno <- VariantAnnotation::geno(vcf)                # genotype info
        geno.matrix <- t( var.geno$GT )                         # genotype calls
        
        # Set indices of SNP variants.
        snp.indices <- which( VariantAnnotation::isSNV(vcf) )
        num.snps <- length(snp.indices)
        stopifnot( num.snps > 0 )
        
        # Filter genotype data to retain only SNP variants.
        geno.matrix <- geno.matrix[, snp.indices, drop=FALSE]
        
        # Combine REF and ALT alleles for each SNP.
        var.alleles <- lapply(snp.indices, function(j)
            c(var.ref[j], var.alt[[j]]))
        
        # Create dataframe with variant locus info.
        snp.loc <- data.frame(chr=var.seqs[snp.indices],
            pos=var.pos[snp.indices])
        
        # Sanity check for ordered VCF records.
        if ( any( sapply( unique(var.seqs), function(var.seq)
            is.unsorted(snp.loc[snp.loc$chr == var.seq, 'pos']) ) ) ) {
            stop("unordered variants in file - '", infile, "'")
        }
        
        # Get default marker IDs for SNPs in this file.
        file.snps <- makeDefaultMarkerIDs( as.mapframe(snp.loc,
            map.unit='bp') )
        
        # Check for multiple variants coinciding at same locus.
        # TODO: handle coinciding variants?
        if ( anyDuplicated(file.snps) ) {
            stop("coinciding variants in file - '", infile, "'")
        }
        
        colnames(geno.matrix) <- file.snps
        
        # Resolve raw genotypes for each variant as concatenated base calls.
        for ( j in 1:num.snps ) {
            
            # Get vector of alleles for this variant.
            variant.alleles <- var.alleles[[j]]
            
            # Get genotype strings for each sample in this variant record.
            geno.data <- unname( geno.matrix[, j] )
            
            # Get mask of genotypes that may have been called.
            # NB: this will misidentify diploid/polyploid missing values
            # as being genotyped, but such cases are handled below.
            called <- geno.data != const$vcf.missing.value
            
            # Set uncalled genotypes to missing value.
            geno.data[ ! called ] <- const$missing.value
            
            if ( any(called) ) {
                
                # Resolve previously unseen genotype strings
                # to their corresponding allele indices.
                for ( unique.call in unique(geno.data[called]) ) {
                    if ( ! unique.call %in% names(geno.memo) ) {
                        indicesC <- unlist( strsplit(unique.call, '[/|]') )
                        indices0 <- suppressWarnings( as.integer(indicesC) )
                        indices1 <- indices0 + 1
                        geno.memo[[unique.call]] <- indices1
                    }
                }
                
                # Get resolved allele indices for called genotypes.
                allele.list <- lapply(geno.data[called], function(geno.call)
                    geno.memo[[geno.call]])
                
                # Check for consistent ploidy.
                ploidy <- unique( lengths(allele.list) )
                if ( length(ploidy) > 1 ) {
                    k <- snp.indices[j]
                    stop("mixed ploidy in variant at position ", var.pos[k],
                        " of ", var.seqs[k], " in file - '", infile, "'")
                }
                
                # Create matrix of allele indices, set alleles (or missing values).
                allele.matrix <- matrix(unlist(allele.list), ncol=ploidy, byrow=TRUE)
                for ( a in seq_along(variant.alleles) ) {
                    allele.matrix[ allele.matrix == a ] <- variant.alleles[a]
                }
                allele.matrix[ is.na(allele.matrix) ] <- const$missing.value
                
                # Concatenate alleles in each called genotype.
                geno.data[called] <- apply(allele.matrix, 1, paste0, collapse='')
            }
            
            geno.matrix[, j] <- geno.data
        }
        
        # Get symbols from genotype matrix.
        g.symbols <- unique( as.character(geno.matrix) )
        
        # Get genotype symbols.
        genotypes <- g.symbols[ g.symbols != const$missing.value ]
        
        # Get characters in genotype symbols.
        a.symbols <- unique( unlist( strsplit(genotypes, '') ) )
        
        # Get allele symbols.
        alleles <- a.symbols[ a.symbols != const$missing.value ]
        
        unknown <- alleles[ ! isRawAllele(alleles) ]
        if ( length(unknown) > 0 ) {
            stop("unknown raw SNP alleles (", toString(unknown),
                ") in file - '", infile, "'")
        }
        
        # If genotype qualities available, create  
        # object with quality-scaled genotypes..
        if (reading.qualities) {
            
            # Set vector of variant error probabilities.
            variant.error.probs <- 10 ^ ( -0.1 * var.qual[snp.indices] )
            
            # Set matrix of genotype error probabilities.
            genotype.error.probs <- 10 ^ ( -0.1 * var.geno$GQ[snp.indices, , drop=FALSE] )
            
            # Calculate error probability for each sample genotype
            # from converted variant/genotype quality scores.
            qual.matrix <- sapply( seq_along(variant.error.probs), function(i)
                variant.error.probs[i] + genotype.error.probs[i, ])
            
            # If qual matrix simplified, restore to matrix.
            if ( ! is.matrix(qual.matrix) ) {
                qual.matrix <- matrix(qual.matrix, ncol=length(variant.error.probs))
            }
            
            # Replace NA values with maximum error probability.
            qual.matrix[ is.na(qual.matrix) ] <- const$qual$prob$range[2]
            
            # Clamp quality scores within range of error probabilities.
            qual.matrix <- sapply(seq_along(file.snps), function(i)
                clamp(qual.matrix[, i], const$qual$prob$range))
            
            # If qual matrix simplified, restore to matrix.
            if ( ! is.matrix(qual.matrix) ) {
                qual.matrix <- matrix(qual.matrix, ncol=length(file.snps))
            }
            
            geno.data <- c(geno.matrix, qual.matrix)
            
        } else { # ..otherwise create object with only genotypes.
            
            geno.data <- geno.matrix
        }
        
        # Prep file variant data.
        data.names <- list(file.samples, file.snps, slices)
        data.shape <- lengths(data.names)
        
        # Set file variant data.
        snps[[i]] <- array(geno.data, dim=data.shape, dimnames=data.names)
    }
    
    # If multiple relevant files, combine SNP genotypes,
    # taking the union of all variants in input files..
    if ( length(infiles) > 1 ) {
        
        # Prep combined variant data.
        combined.samples <- sort( unique( unlist( lapply(snps, rownames) ) ) )
        combined.snps <- sort( unique( unlist( lapply(snps, colnames) ) ) )
        combined.names <- list(combined.samples, combined.snps, slices)
        combined.shape <- lengths(combined.names)
        
        # Init combined variant data.
        result <- array(const$missing.value, dim=combined.shape, dimnames=combined.names)
        
        # Set combined variant data from each input file.
        # NB: we previously checked for duplicate samples and coinciding
        # variants, so we can assume that there will be no conflicts.
        for ( i in seq_along(snps) ) {
            for ( slice in slices ) {
                for ( sample.id in rownames(snps[[i]]) ) {
                    for ( snp.id in colnames(snps[[i]]) ) {
                        result[sample.id, snp.id, slice] <- snps[[i]][sample.id, snp.id, slice]
                    }
                }
            }
        }
        
    } else { # ..otherwise take single set of SNP genotypes.
        
        result <- snps[[1]]
    }
    
    # Get combined number of SNP variants.
    num.snps <- ncol(result)
    stopifnot( num.snps > 0 )
    
    # If filters specified, filter SNP variants.
    if ( require.all || require.any || require.polymorphic ) {
        
        mask <- rep(TRUE, num.snps)
        
        if (require.all) { # Remove incomplete variants.
            mask <- mask & sapply(1:num.snps, function(i)
                all(result[, i, 'geno'] != const$missing.value))
        }
        
        if (require.any) { # Remove variants without genotypes.
            mask <- mask & sapply(1:num.snps, function(i)
                any(result[, i, 'geno'] != const$missing.value))
        }
        
        if (require.polymorphic) { # Remove monomorphic variants.
            mask <- mask & sapply(1:num.snps, function(i) {
                geno.data <- result[, i, 'geno']
                g.symbols <- unique(geno.data)
                genotypes <- g.symbols[ g.symbols != const$missing.value ]
                return( length(genotypes) > 1 )
            })
        }
        
        # Apply filter mask.
        result <- result[, mask, , drop=FALSE]
    }
    
    return(result)
}

# readGenoVCF ------------------------------------------------------------------
#' Read genotype data from a VCF file.
#' 
#' This function reads SNP genotype data from one or more VCF files, and
#' returns these as an \pkg{R/qtl} \code{cross} \code{geno} object.
#' 
#' If no founder samples are specified, this function assigns enumerated
#' genotypes at each locus according to the observed raw SNP genotype. So
#' for example, if the raw SNP alleles at a given locus are \code{'A'} and
#' \code{'C'}, samples are assigned the genotypes \code{'1'} and \code{'2'},
#' respectively.
#' 
#' If founder samples are specified, this function assigns to each marker a
#' genotype symbol consisting of alleles that each correspond to a specific
#' founder.
#' 
#' If the \code{alleles} parameter is specified, this must be a mapping of
#' founder sample IDs to allele symbols (e.g.
#' \code{mapping( c(DBVPG6044 = 'W', YPS128 = 'S') )}). If the \code{alleles}
#' parameter is not specified, allele symbols are taken from the letters of the
#' alphabet (i.e. \code{'A'}, \code{'B'} etc.).
#' 
#' @param ... Input VCF file paths.
#' @param samples Cross sample IDs.
#' @param founders Founder sample IDs.
#' @param alleles Mapping of founder sample IDs to founder allele symbols.
#' 
#' @return An \pkg{R/qtl} \code{cross} \code{geno} object.
#' 
#' @template section-geno-raw
#' @template section-geno-enum
#' @template section-geno-founder
#' 
#' @export
#' @family VCF functions
#' @rdname readGenoVCF
readGenoVCF <- function(..., samples, founders=NULL, alleles=NULL) {
    
    # TODO: optimise.
    
    sample.data <- readSnpsVCF(..., samples=samples, require.any=TRUE,
        require.polymorphic=TRUE)
    
    if ( ! is.null(founders) ) {
        
        clashing <- intersect(samples, founders)
        if ( length(clashing) > 0 ) {
            stop("clashing sample IDs - ", toString(clashing), "'")
        }
        
        founder.data <- readSnpsVCF(..., samples=founders,
            require.all=TRUE, require.polymorphic=TRUE)
    
    } else {
        
        founder.data <- NULL
    }
    
    geno <- makeGeno(sample.data, founder.data, alleles=alleles)
    
    return(geno)
}

# End of vcf.R #################################################################