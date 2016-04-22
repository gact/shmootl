# Start of CrossInfo-class.R ###################################################

# CrossInfo --------------------------------------------------------------------
#' An S4 class to hold yeast cross information.
#' 
#' A CrossInfo object holds yeast cross information for a specific \code{cross}
#' object. The contents of its slots should match its corresponding object.
#'
#' @slot seq A character vector of sequence identifiers, with the name
#' of each element being the name of the given sequence. These must be
#' unique.
#' 
#' @slot pheno A vector of cross phenotypes, with the name of each element being 
#' the syntactically valid name of the phenotype ID (as produced by the 
#' function \code{make.names}). Phenotype IDs must be unique.
#' 
#' @slot markers A \code{data.frame} with one required column - 'marker' - 
#' containing the marker ID, and one optional column: 'seq'. If present, the 
#' 'seq' column contains the sequence corresponding to the given marker. 
#' Each row name contains the syntactically valid name of the marker ID (as 
#' produced by the function \code{make.names}). Marker IDs must be unique.
#' 
#' @slot samples A \code{data.frame} with one required column - 'sample.index' - 
#' that contains indices of the samples in the given cross dataset. Other optional 
#' columns can be present, including 'sample.id', 'sample.name', 'strain.index' 
#' and 'tetrad.index'. 
#' 
#' The 'sample.id' column contains sample IDs, which can be any valid item ID.
#' If this column is present, another column called 'sample.name' will
#' contain the syntactically valid name of the sample ID (as produced by the 
#' function \code{make.names}). Duplicate sample IDs are permissible, but only
#' if referring to replicate samples of the same strain. Different strains can 
#' have different numbers of replicates, but samples from a given strain must 
#' be in consecutive rows.
#'  
#' The 'strain.index' column, if present, contains strain indices denoting 
#' which strain each sample is taken from. This is only necessary if the cross
#' contains sample replicates. Replicate samples of a given strain must have 
#' the same genotype data, and their sample IDs must match, if present. 
#' 
#' If the samples are tetradic, this can also include a 'tetrad.index' column 
#' containing the tetrad index of each sample, indicating which tetrad that 
#' sample is from. Tetrads are assumed to be formed of each group of four 
#' consecutive strains, unless sample IDs are present that indicate tetrad 
#' membership. Samples can be missing from a tetrad, but in such cases the 
#' sample IDs should be such that it is possible to assign each sample to
#' a tetrad. In any case, samples from the same tetrad must be in consecutive 
#' rows of the cross data.
#' 
#' @slot alleles A vector of cross allele symbols. 
#'  
#' @docType class
#' @export
#' @keywords internal
#' @rdname CrossInfo-class
CrossInfo <- setClass('CrossInfo',
    
    slots = c( 
        seq = 'character', 
        pheno = 'character', 
        markers = 'data.frame',
        samples = 'data.frame',
        alleles = 'character'
    ),
    
    prototype=list( 
        seq = character(),
        pheno = character(),
        markers = data.frame(marker=character(), 
            stringsAsFactors=FALSE),
        samples = data.frame(sample.index=integer(), 
            stringsAsFactors=FALSE),
        alleles = character()
    )
)

# getAlleles -------------------------------------------------------------------
#' Get allele symbols.
#' 
#' @template param-CrossInfo
#'  
#' @return Character vector of allele symbols.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getAlleles-methods
setGeneric('getAlleles', function(cross.info) { 
    standardGeneric('getAlleles') })

# CrossInfo::getAlleles --------------------------------------------------------
#' @aliases getAlleles,CrossInfo-method
#' @export
#' @rdname getAlleles-methods
setMethod('getAlleles', signature='CrossInfo', 
    definition = function(cross.info) { 
        return(cross.info@alleles)
})

# getMarkerIndices -------------------------------------------------------------
#' Get marker indices.
#' 
#' @template param-CrossInfo
#' @template param-markers
#' 
#' @return Integer vector of marker indices.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getMarkerIndices-methods
setGeneric('getMarkerIndices', function(cross.info, markers=NULL) { 
  standardGeneric('getMarkerIndices') })

# CrossInfo::getMarkerIndices --------------------------------------------------
#' @aliases getMarkerIndices,CrossInfo-method
#' @export
#' @rdname getMarkerIndices-methods
setMethod('getMarkerIndices', signature='CrossInfo', 
  definition = function(cross.info, markers=NULL) { 

    if ( ! is.null(markers) ) {
        
        stopifnot( length(markers) > 0 )
        stopif( anyNA(markers) )
        
        if ( is.integer(markers) ) {
            
            exrange <- markers[ markers < 1 | markers > nrow(cross.info@markers) ]
            if ( length(exrange) > 0 ) {
                stop("marker indices out of range - '", toString(exrange), "'")
            }
            
            indices <- unname(markers)
            
        } else if ( is.logical(markers) ) {
            
            if ( length(markers) != nrow(cross.info@markers) ) {
                stop("marker logical vector length mismatch")
            }
            
            indices <- unname( which(markers) )
            
        } else if ( is.character(markers) ) {
            
            index.list <- lapply(markers, function(marker) 
                which( cross.info@markers$marker == marker | 
                rownames(cross.info@markers) == marker ) )
            
            unfound <- markers[ lengths(index.list) == 0 ]
            if ( length(unfound) > 0 ) {
                stop("markers not found - '", toString(unfound), "'")
            }
            
            indices <- unlist(index.list)
            
        } else {
            
            stop("marker vector must be of type logical, integer, or character")
        }

    } else {
        
        indices <- getRowIndices(cross.info@markers)
    }
    
    return(indices) 
})

# getMarkerNames ---------------------------------------------------------------
#' Get marker names.
#' 
#' @template param-CrossInfo
#' @template param-markers
#' 
#' @return Character vector of marker names.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getMarkerNames-methods
setGeneric('getMarkerNames', function(cross.info, markers=NULL) {
  standardGeneric('getMarkerNames') })

# CrossInfo::getMarkerNames ----------------------------------------------------
#' @aliases getMarkerNames,CrossInfo-method
#' @export
#' @rdname getMarkerNames-methods
setMethod('getMarkerNames', signature='CrossInfo', 
  definition = function(cross.info, markers=NULL) { 
    indices <- getMarkerIndices(cross.info, markers)
    return( rownames(cross.info@markers)[indices] )
})

# getMarkers -------------------------------------------------------------------
#' Get marker IDs.
#' 
#' @template param-CrossInfo
#' @template param-markers
#' 
#' @return Character vector of marker IDs.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getMarkers-methods
setGeneric('getMarkers', function(cross.info, markers=NULL) {
  standardGeneric('getMarkers') })

# CrossInfo::getMarkers --------------------------------------------------------
#' @aliases getMarkers,CrossInfo-method
#' @export
#' @rdname getMarkers-methods
setMethod('getMarkers', signature='CrossInfo', 
  definition = function(cross.info, markers=NULL) { 
    indices <- getMarkerIndices(cross.info, markers)
    return( cross.info@markers$marker[indices] )
})

# getMarkerSeqs ----------------------------------------------------------------
#' Get sequences by marker.
#' 
#' @template param-CrossInfo
#' @template param-markers
#' 
#' @return Vector of sequences corresponding to the given markers. Returns 
#' \code{NULL} if no sequence-marker information is available.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getMarkerSeqs-methods
setGeneric('getMarkerSeqs', function(cross.info, markers=NULL) { 
    standardGeneric('getMarkerSeqs') })

# CrossInfo::getMarkerSeqs -----------------------------------------------------
#' @aliases getMarkerSeqs,CrossInfo-method
#' @export
#' @rdname getMarkerSeqs-methods
setMethod('getMarkerSeqs', signature='CrossInfo', 
    definition = function(cross.info, markers=NULL) { 
              
    if ( 'seq' %in% colnames(cross.info@markers) ) {
        marker.indices <- getMarkerIndices(cross.info, markers)
        seqs <- cross.info@markers$seq[marker.indices]
    } else {
        seqs <- NULL
    }          
              
    return(seqs)
})

# getNumMarkers ----------------------------------------------------------------
#' Get the number of markers.
#' 
#' @template param-CrossInfo
#' 
#' @return Number of markers.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getNumMarkers-methods
setGeneric('getNumMarkers', function(cross.info) {
  standardGeneric('getNumMarkers') })

# CrossInfo::getNumMarkers -----------------------------------------------------
#' @aliases getNumMarkers,CrossInfo-method
#' @export
#' @rdname getNumMarkers-methods
setMethod('getNumMarkers', signature='CrossInfo', 
  definition = function(cross.info) { 
    return( nrow(cross.info@markers) )
})

# getNumPhenotypes -------------------------------------------------------------
#' Get the number of phenotypes.
#' 
#' @template param-CrossInfo
#' 
#' @return Number of phenotypes.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getNumPhenotypes-methods
setGeneric('getNumPhenotypes', function(cross.info) {
  standardGeneric('getNumPhenotypes') })

# CrossInfo::getNumPhenotypes --------------------------------------------------
#' @aliases getNumPhenotypes,CrossInfo-method
#' @export
#' @rdname getNumPhenotypes-methods
setMethod('getNumPhenotypes', signature='CrossInfo', 
  definition = function(cross.info) { 
    return( length(cross.info@pheno) )
})

# getNumSamples ----------------------------------------------------------------
#' Get the number of samples.
#' 
#' @template param-CrossInfo
#' 
#' @return Number of samples.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getNumSamples-methods
setGeneric('getNumSamples', function(cross.info) {
  standardGeneric('getNumSamples') })

# CrossInfo::getNumSamples -----------------------------------------------------
#' @aliases getNumSamples,CrossInfo-method
#' @export
#' @rdname getNumSamples-methods
setMethod('getNumSamples', signature='CrossInfo', 
  definition = function(cross.info) {  
    return( nrow(cross.info@samples) )
})

# getNumSeqs -------------------------------------------------------------------
#' Get the number of sequences.
#' 
#' @template param-CrossInfo
#' 
#' @return Number of sequences.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getNumSeqs-methods
setGeneric('getNumSeqs', function(cross.info) {
    standardGeneric('getNumSeqs') })

# CrossInfo::getNumSeqs --------------------------------------------------------
#' @aliases getNumSeqs,CrossInfo-method
#' @export
#' @rdname getNumSeqs-methods
setMethod('getNumSeqs', signature='CrossInfo', 
    definition = function(cross.info) { 
    return( length(cross.info@seq) )
})

# getPhenotypeIndices ----------------------------------------------------------
#' Get phenotype indices.
#' 
#' @template param-CrossInfo
#' @template param-phenotypes
#'  
#' @return Integer vector of phenotype indices.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getPhenotypeIndices-methods
setGeneric('getPhenotypeIndices', function(cross.info, phenotypes=NULL) { 
    standardGeneric('getPhenotypeIndices') })

# CrossInfo::getPhenotypeIndices -----------------------------------------------
#' @aliases getPhenotypeIndices,CrossInfo-method
#' @export
#' @rdname getPhenotypeIndices-methods
setMethod('getPhenotypeIndices', signature='CrossInfo', 
    definition = function(cross.info, phenotypes=NULL) { 
    
    if ( ! is.null(phenotypes) ) {
        
        stopifnot( length(phenotypes) > 0 )
        stopif( anyNA(phenotypes) )
    
        if ( is.integer(phenotypes) ) {
        
            exrange <- phenotypes[ phenotypes < 1 | phenotypes > length(cross.info@pheno) ]
            if ( length(exrange) > 0 ) {
                stop("phenotype indices out of range - '", toString(exrange), "'")
            }
            
            indices <- unname(phenotypes)

        } else if ( is.logical(phenotypes) ) {
             
            if ( length(phenotypes) != length(cross.info@pheno) ) {
                stop("phenotype logical vector length mismatch")
            }
      
            indices <- unname( which(phenotypes) )
            
        } else if ( is.character(phenotypes) ) {
            
            index.list <- lapply(phenotypes, function(phenotype) 
                which( cross.info@pheno == phenotype | 
                names(cross.info@pheno) == phenotype ) )
            
            unfound <- phenotypes[ lengths(index.list) == 0 ]
            if ( length(unfound) > 0 ) {
                stop("phenotypes not found - '", toString(unfound), "'")
            }
            
            indices <- unlist(index.list)
            
        } else {
        
            stop("phenotype vector must be of type logical, integer, or character")
        }
    
    } else {
        
        indices <- getIndices(cross.info@pheno)
    }
        
    return(indices)
})

# getPhenotypeNames ------------------------------------------------------------
#' Get phenotype names.
#' 
#' @template param-CrossInfo
#' @template param-phenotypes
#'  
#' @return Character vector of phenotype names.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getPhenotypeNames-methods
setGeneric('getPhenotypeNames', function(cross.info, phenotypes=NULL) { 
  standardGeneric('getPhenotypeNames') })

# CrossInfo::getPhenotypeNames -------------------------------------------------
#' @aliases getPhenotypeNames,CrossInfo-method
#' @export        
#' @rdname getPhenotypeNames-methods
setMethod('getPhenotypeNames', signature='CrossInfo', 
  definition = function(cross.info, phenotypes=NULL) { 
    indices <- getPhenotypeIndices(cross.info, phenotypes)
    return( names(cross.info@pheno)[indices] )
})

# getPhenotypes ----------------------------------------------------------------
#' Get phenotype IDs.
#' 
#' @template param-CrossInfo
#' @template param-phenotypes
#'  
#' @return Character vector of phenotype IDs.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getPhenotypes-methods
setGeneric('getPhenotypes', function(cross.info, phenotypes=NULL) { 
  standardGeneric('getPhenotypes') })

# CrossInfo::getPhenotypes -----------------------------------------------------
#' @aliases getPhenotypes,CrossInfo-method
#' @export
#' @rdname getPhenotypes-methods
setMethod('getPhenotypes', signature='CrossInfo', 
  definition = function(cross.info, phenotypes=NULL) { 
    indices <- getPhenotypeIndices(cross.info, phenotypes)
    return( unname(cross.info@pheno[indices]) )
})

# getSampleIndices -------------------------------------------------------------
#' Get sample indices.
#' 
#' Get indices of the specified samples. Samples can be specified at one of the 
#' following levels: sample, strain, or tetrad. If none are specified, all 
#' sample indices are returned.
#' 
#' @template param-CrossInfo
#' @template param-samples
#' @template param-strains
#' @template param-tetrads
#' @param simplify Simplify list return value to a vector.
#' 
#' @return If samples are specified by strain or tetrad, and if the option 
#' \code{simplify} is not TRUE, this method returns a list of integer vectors,
#' each containing the sample indices corresponding to a given strain/tetrad.
#' Otherwise, a vector of sample indices is returned.
#'  
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getSampleIndices-methods
setGeneric('getSampleIndices', function(cross.info, samples=NULL, strains=NULL, 
    tetrads=NULL, simplify=FALSE) { standardGeneric('getSampleIndices') })

# CrossInfo::getSampleIndices --------------------------------------------------
#' @aliases getSampleIndices,CrossInfo-method
#' @export
#' @rdname getSampleIndices-methods
setMethod('getSampleIndices', signature='CrossInfo', 
    definition = function(cross.info, samples=NULL, strains=NULL, tetrads=NULL, 
    simplify=FALSE) { 
    
    sample.args <- list(sample=samples, strain=strains, tetrad=tetrads)
    
    sample.args <- sample.args[ ! sapply(sample.args, is.null) ]
    
    if ( length(sample.args) > 1 ) {
        stop("samples can be specified at only one of the following levels - ", toString(const$sample.levels))
    }
    
    if ( length(sample.args) == 1 ) {
        
        k <- names(sample.args)
        
        X <- sample.args[[k]]
        
        stopifnot( length(X) > 0 )
        stopif( anyNA(X) )
        
        index.key <- const$sample.aspects[k, 'index']
        
        if ( ! index.key %in% colnames(cross.info@samples) ) {
            stop(k, " index column not found")
        }
        
        max.index <- max(cross.info@samples[, index.key])
        
        if ( is.integer(X) ) {
            
            exrange <- samples[ X < 1 | X > max.index ]
            if ( length(exrange) > 0 ) {
                stop(k, " indices out of range - '", toString(exrange), "'")
            }

            index.list <- as.list( unname(X) )
            
            index.list <- lapply(index.list, function(i) 
                which( cross.info@samples[, index.key] == i ) )
            
        } else if ( is.logical(X) ) {
            
            if ( length(X) != max.index ) {
                stop(k, " logical vector length mismatch")
            }
            
            index.list <- as.list( unname( which(X) ) )
            
            index.list <- lapply(index.list, function(i) 
                which( cross.info@samples[, index.key] == i ) )
            
        } else if ( is.character(X) ) {
            
            headings <- const$sample.aspects[k, c('id', 'name')]
  
            for ( h in headings ) {
                
                if ( is.na(h) ) {
                    stop(k, " ", names(h), " column not supported")
                }
                
                if ( ! h %in% colnames(cross.info@samples) ) {
                    stop(k, " ", names(h), " column not found")
                }
            }
            
            index.list <- lapply(X, function(x)
                which( cross.info@samples[[ headings['id'] ]] == x | 
                cross.info@samples[[ headings['name'] ]] == x ) )
            
            unfound <- X[ lengths(index.list) == 0 ]
            if ( length(unfound) > 0 ) {
                stop(k, " values not found - '", toString(unfound), "'")
            }
            
        } else {
            
            stop("sample vector must be of type logical, integer, or character")
        }        

        if( k == 'sample' || simplify ) {
            indices <- unlist(index.list)
        } else {
            indices <- index.list
        }
        
    } else {
        
        indices <- unname(cross.info@samples$sample.index)
    }
    
    return(indices) 
})

# getSampleNames ---------------------------------------------------------------
#' Get sample names.
#' 
#' Get names of the specified samples. Samples can be specified at one of the 
#' following levels: sample, strain, or tetrad. If none are specified, all 
#' sample names are returned.
#' 
#' @template param-CrossInfo
#' @template param-samples
#' @template param-strains
#' @template param-tetrads
#' @param simplify Simplify list return value to a vector.
#' 
#' @return If samples are specified by strain or tetrad, and if the option 
#' \code{simplify} is not TRUE, this method returns a list of character vectors,
#' each containing the sample names corresponding to a given strain/tetrad.
#' Otherwise, a vector of sample names is returned.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getSampleNames-methods
setGeneric('getSampleNames', function(cross.info, samples=NULL, strains=NULL, 
    tetrads=NULL, simplify=FALSE) { standardGeneric('getSampleNames') })

# CrossInfo::getSampleNames ----------------------------------------------------
#' @aliases getSampleNames,CrossInfo-method
#' @export
#' @rdname getSampleNames-methods
setMethod('getSampleNames', signature='CrossInfo', 
    definition = function(cross.info, samples=NULL, strains=NULL, tetrads=NULL, 
    simplify=FALSE) { 
      
    sample.indices <- getSampleIndices(cross.info, samples=samples, 
        strains=strains, tetrads=tetrads)
    
    if ( hasSampleIDs(cross.info) ) {
        
        if ( is.list(sample.indices) ) {
            sample.names <- lapply(sample.indices, function(I)
                cross.info@samples$sample.name[I])
        } else {
            sample.names <- sapply(sample.indices, function(I)
                cross.info@samples$sample.name[I])
        }
        
    } else {
        sample.names <- NULL
    }
    
    if (simplify) {
        sample.names <- unlist(sample.names)
    }
    
    return(sample.names) 
})

# getSamples -------------------------------------------------------------------
#' Get sample IDs.
#' 
#' Get IDs of the specified samples. Samples can be specified at one of the 
#' following levels: sample, strain, or tetrad. If none are specified, all 
#' sample IDs are returned.
#' 
#' @template param-CrossInfo
#' @template param-samples
#' @template param-strains
#' @template param-tetrads
#' @param simplify Simplify list return value to a vector.
#' 
#' @return If samples are specified by strain or tetrad, and if the option 
#' \code{simplify} is not TRUE, this method returns a list of character vectors,
#' each containing the sample IDs corresponding to a given strain/tetrad.
#' Otherwise, a vector of sample IDs is returned.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getSamples-methods
setGeneric('getSamples', function(cross.info, samples=NULL, strains=NULL, 
    tetrads=NULL, simplify=FALSE) { standardGeneric('getSamples') })

# CrossInfo::getSamples --------------------------------------------------------
#' @aliases getSamples,CrossInfo-method
#' @export
#' @rdname getSamples-methods
setMethod('getSamples', signature='CrossInfo', 
    definition = function(cross.info, samples=NULL, strains=NULL, tetrads=NULL, 
    simplify=FALSE) { 
              
    sample.indices <- getSampleIndices(cross.info, samples=samples, 
        strains=strains, tetrads=tetrads)
              
    if ( hasSampleIDs(cross.info) ) {
        
        if ( is.list(sample.indices) ) {
            sample.ids <- lapply(sample.indices, function(I)
                cross.info@samples$sample.id[I])
        } else {
            sample.ids <- sapply(sample.indices, function(I)
                cross.info@samples$sample.id[I])
        }
        
    } else {
        sample.ids <- NULL
    }

    if (simplify) {
        sample.ids <- unlist(sample.ids)
    }    
         
    return(sample.ids)
})

# getSeqIndices ----------------------------------------------------------------
#' Get sequence indices.
#' 
#' @template param-CrossInfo
#' @template param-sequences
#' 
#' @return Integer vector of sequence indices.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getSeqIndices-methods
setGeneric('getSeqIndices', function(cross.info, sequences=NULL) { 
    standardGeneric('getSeqIndices') })

# CrossInfo::getSeqIndices -----------------------------------------------------
#' @aliases getSeqIndices,CrossInfo-method
#' @export
#' @rdname getSeqIndices-methods
setMethod('getSeqIndices', signature='CrossInfo', 
    definition = function(cross.info, sequences=NULL) { 
  
    if ( ! is.null(sequences) ) {
      
        stopifnot( length(sequences) > 0 )
        stopif( anyNA(sequences) )
      
        if ( is.integer(sequences) ) {
          
            exrange <- sequences[ sequences < 1 | sequences > length(cross.info@seq) ]
            if ( length(exrange) > 0 ) {
                stop("sequence indices out of range - '", toString(exrange), "'")
            }
            
            indices <- unname(sequences)
          
        } else if ( is.logical(sequences) ) {
          
            if ( length(sequences) != length(cross.info@seq) ) {
                stop("sequence logical vector length mismatch")
            }
          
            indices <- unname( which(sequences) )
          
        } else if ( is.character(sequences) ) {
          
            norm.seqs <- normSeq(sequences)
            
            index.list <- lapply(norm.seqs, function(norm.seq) 
                which( cross.info@seq == norm.seq ) )
          
            unfound <- sequences[ lengths(index.list) == 0 ]
            if ( length(unfound) > 0 ) {
                stop("sequences not found - '", toString(unfound), "'")
            }
          
            indices <- unlist(index.list)
          
        } else {
          
            stop("sequence vector must be of type logical, integer, or character")
        }
      
    } else {
      
        indices <- getIndices(cross.info@seq)
    }
  
    return(indices) 
})

# getSeqMarkers ----------------------------------------------------------------
#' Get markers by sequence.
#' 
#' @template param-CrossInfo
#' @template param-sequences
#' @param simplify Simplify list return value to a vector.
#' 
#' @return List of character vectors, each containing the marker IDs
#' corresponding to a given sequence. If option \code{simplify} is TRUE, 
#' this is simplified to a vector. Returns \code{NULL} if no sequence-marker 
#' information is available.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getSeqMarkers-methods
setGeneric('getSeqMarkers', function(cross.info, sequences=NULL, 
    simplify=FALSE) { standardGeneric('getSeqMarkers') })

# CrossInfo::getSeqMarkers -----------------------------------------------------
#' @aliases getSeqMarkers,CrossInfo-method
#' @export
#' @rdname getSeqMarkers-methods
setMethod('getSeqMarkers', signature='CrossInfo', 
    definition = function(cross.info, sequences=NULL, simplify=FALSE) { 
            
    if ( 'seq' %in% colnames(cross.info@markers) ) {
        sequences <- getSequences(cross.info, sequences)
        markers <- lapply(sequences, function(s) 
            cross.info@markers$marker[ cross.info@markers$seq == s ] )
    } else {
        markers <- NULL
    }
            
    if (simplify) {
        markers <- unlist(markers)
    }
            
    return(markers)   
})

# getSeqNames ------------------------------------------------------------------
#' Get vector of sequence names.
#' 
#' @template param-CrossInfo
#' @template param-sequences
#' 
#' @return Vector of sequence names.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getSeqNames-methods
setGeneric('getSeqNames', function(cross.info, sequences=NULL) { 
    standardGeneric('getSeqNames') })

# CrossInfo::getSeqNames -------------------------------------------------------
#' @aliases getSeqNames,CrossInfo-method
#' @export
#' @rdname getSeqNames-methods
setMethod('getSeqNames', signature='CrossInfo', 
    definition = function(cross.info, sequences=NULL) { 
    indices <- getSeqIndices(cross.info, sequences)
    return( names(cross.info@seq[indices]) )
})

# getSequences -----------------------------------------------------------------
#' Get vector of normalised sequence labels.
#' 
#' @template param-CrossInfo
#' @template param-sequences
#' 
#' @return Character vector of normalised sequence labels.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getSequences-methods
setGeneric('getSequences', function(cross.info, sequences=NULL) { 
    standardGeneric('getSequences') })

# CrossInfo::getSequences ------------------------------------------------------
#' @aliases getSequences,CrossInfo-method
#' @export
#' @rdname getSequences-methods
setMethod('getSequences', signature='CrossInfo', 
    definition = function(cross.info, sequences=NULL) { 
    indices <- getSeqIndices(cross.info, sequences)
    return(cross.info@seq[indices])
})

# getStrainIndices -------------------------------------------------------------
#' Get strain indices for the given samples.
#' 
#' @template param-CrossInfo
#' @template param-samples
#' 
#' @return Integer vector of strain indices for the given samples. If no samples 
#' are specified, all strain indices are returned. Returns \code{NULL} if strain 
#' indices are not present.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getStrainIndices-methods
setGeneric('getStrainIndices', function(cross.info, samples=NULL) { 
    standardGeneric('getStrainIndices') })

# CrossInfo::getStrainIndices --------------------------------------------------
#' @aliases getStrainIndices,CrossInfo-method
#' @export
#' @rdname getStrainIndices-methods
setMethod('getStrainIndices', signature='CrossInfo', 
    definition = function(cross.info, samples=NULL) { 
              
    sample.indices <- getSampleIndices(cross.info, samples=samples)
              
    if ( hasStrainIndices(cross.info) ) {
        strain.indices <- cross.info@samples$strain.index[sample.indices]
    } else {
        strain.indices <- cross.info@samples$sample.index[sample.indices]
    }
              
    return(strain.indices)              
})   

# getTetradIndices -------------------------------------------------------------
#' Get tetrad indices for the given samples.
#' 
#' @template param-CrossInfo
#' @template param-samples
#' 
#' @return Integer vector of tetrad indices for the given samples. If no samples 
#' are specified, all tetrad indices are returned. Returns \code{NULL} if tetrad 
#' indices are not present.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname getTetradIndices-methods
setGeneric('getTetradIndices', function(cross.info, samples=NULL) { 
  standardGeneric('getTetradIndices') })

# CrossInfo::getTetradIndices --------------------------------------------------
#' @aliases getTetradIndices,CrossInfo-method
#' @export
#' @rdname getTetradIndices-methods
setMethod('getTetradIndices', signature='CrossInfo', 
  definition = function(cross.info, samples=NULL) { 
    
    sample.indices <- getSampleIndices(cross.info, samples=samples)
    
    if ( hasTetradIndices(cross.info) ) {
        tetrad.indices <- cross.info@samples$tetrad.index[sample.indices]
    } else {
        tetrad.indices <- NULL
    }
    
    return(tetrad.indices)   
})

# hasMarkerSeqs ----------------------------------------------------------------
#' Test if \code{CrossInfo} object has a marker-sequence mapping.
#' 
#' @template param-CrossInfo
#' 
#' @return TRUE if marker-sequence mapping is present; FALSE otherwise.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname hasMarkerSeqs-methods
setGeneric('hasMarkerSeqs', function(cross.info) { 
    standardGeneric('hasMarkerSeqs') })

# CrossInfo::hasMarkerSeqs -----------------------------------------------------
#' @aliases hasMarkerSeqs,CrossInfo-method
#' @export
#' @rdname hasMarkerSeqs-methods
setMethod('hasMarkerSeqs', signature='CrossInfo', 
    definition = function(cross.info) {
    return( 'seq' %in% colnames(cross.info@markers) )
})

# hasSampleIDs -----------------------------------------------------------------
#' Test if \code{CrossInfo} object has sample IDs.
#' 
#' @template param-CrossInfo
#' 
#' @return TRUE if sample IDs are present; FALSE otherwise.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname hasSampleIDs-methods
setGeneric('hasSampleIDs', function(cross.info) { 
    standardGeneric('hasSampleIDs') })

# CrossInfo::hasSampleIDs ------------------------------------------------------
#' @aliases hasSampleIDs,CrossInfo-method
#' @export
#' @rdname hasSampleIDs-methods
setMethod('hasSampleIDs', signature='CrossInfo', 
    definition = function(cross.info) {
    return( 'sample.id' %in% colnames(cross.info@samples) )
})

# hasStrainIndices -------------------------------------------------------------
#' Test if \code{CrossInfo} object has strain indices.
#' 
#' @template param-CrossInfo
#' 
#' @return TRUE if strain indices are present; FALSE otherwise.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname hasStrainIndices-methods
setGeneric('hasStrainIndices', function(cross.info) { 
    standardGeneric('hasStrainIndices') })

# CrossInfo::hasStrainIndices --------------------------------------------------
#' @aliases hasStrainIndices,CrossInfo-method
#' @export
#' @rdname hasStrainIndices-methods
setMethod('hasStrainIndices', signature='CrossInfo', 
    definition = function(cross.info) {
    return( 'strain.index' %in% colnames(cross.info@samples) )
})

# hasTetradIndices -------------------------------------------------------------
#' Test if \code{CrossInfo} object has sample tetrad indices.
#' 
#' @template param-CrossInfo
#' 
#' @return TRUE if sample tetrad indices are present; FALSE otherwise.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname hasTetradIndices-methods
setGeneric('hasTetradIndices', function(cross.info) { 
    standardGeneric('hasTetradIndices') })

# CrossInfo::hasTetradIndices --------------------------------------------------
#' @aliases hasTetradIndices,CrossInfo-method
#' @export
#' @rdname hasTetradIndices-methods
setMethod('hasTetradIndices', signature='CrossInfo', 
    definition = function(cross.info) {
    return( 'tetrad.index' %in% colnames(cross.info@samples) )
})

# setAlleles -------------------------------------------------------------------
#' Set allele symbols.
#' 
#' @template param-CrossInfo
#' @param alleles Vector of allele symbols.
#' 
#' @return Input \code{CrossInfo} object with the given alleles.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname setAlleles-methods
setGeneric('setAlleles', function(cross.info, alleles) { 
    standardGeneric('setAlleles') })

# CrossInfo::setAlleles --------------------------------------------------------
#' @aliases setAlleles,CrossInfo-method
#' @export
#' @rdname setAlleles-methods
setMethod('setAlleles', signature='CrossInfo', 
    definition = function(cross.info, alleles) { 
    cross.info@alleles <- alleles
    validateAlleles(cross.info)
    return(cross.info)
})

# setMarkers -------------------------------------------------------------------
#' Set markers.
#' 
#' @template param-CrossInfo
#' @param markers Vector of marker IDs.
#' 
#' @return Input \code{CrossInfo} object with the given markers.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname setMarkers-methods
setGeneric('setMarkers', function(cross.info, markers) 
    { standardGeneric('setMarkers') })

# CrossInfo::setMarkers --------------------------------------------------------
#' @aliases setMarkers,CrossInfo-method
#' @export
#' @rdname setMarkers-methods
setMethod('setMarkers', signature='CrossInfo', definition = 
  function(cross.info, markers) { 
  
    stopif( missing(markers) )
      
    cross.info@markers <- data.frame( marker=as.character(markers), 
        row.names=make.names(markers), stringsAsFactors=FALSE)
    
    validateMarkers(cross.info)
    
    return(cross.info)
})

# setMarkerSeqs ----------------------------------------------------------------
#' Set marker sequences.
#' 
#' @template param-CrossInfo
#' @param markers Vector containing markers to assign sequences.
#' @param sequences Vector containing sequences for the given markers.
#' 
#' @return Input \code{CrossInfo} object with the given marker-sequence info.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname setMarkerSeqs-methods
setGeneric('setMarkerSeqs', function(cross.info, markers=NULL, sequences=NULL) { 
    standardGeneric('setMarkerSeqs') })

# CrossInfo::setMarkerSeqs -----------------------------------------------------
#' @aliases setMarkerSeqs,CrossInfo-method
#' @export
#' @rdname setMarkerSeqs-methods
setMethod('setMarkerSeqs', signature='CrossInfo', definition = 
    function(cross.info, markers=NULL, sequences=NULL) {
    
    stopif( is.null(sequences) )
        
    marker.indices <- getMarkerIndices(cross.info, markers)
      
    norm.seqs <- normSeq(sequences)
      
    if ( length(norm.seqs) != length(marker.indices) ) {
        stop("cannot assign ", length(norm.seqs), " sequences to ", 
            length(marker.indices)," markers")
    }
      
    if ( ! hasMarkerSeqs(cross.info) ) {
        cross.info@markers$seq <- NA
    }
      
    cross.info@markers$seq[marker.indices] <- norm.seqs
      
    return(cross.info)
})        

# setPhenotypes ----------------------------------------------------------------
#' Set phenotypes.
#' 
#' @template param-CrossInfo
#' @param phenotypes Vector of phenotype IDs.
#' 
#' @return Input \code{CrossInfo} object with the given phenotypes.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname setPhenotypes-methods
setGeneric('setPhenotypes', function(cross.info, phenotypes) {
    standardGeneric('setPhenotypes') })

# CrossInfo::setPhenotypes -----------------------------------------------------
#' @aliases setPhenotypes,CrossInfo-method
#' @export
#' @rdname setPhenotypes-methods
setMethod('setPhenotypes', signature='CrossInfo', 
    definition = function(cross.info, phenotypes) { 
              
    cross.info@pheno <- as.character(phenotypes)
    names(cross.info@pheno) <- make.names(cross.info@pheno)
    
    validatePhenotypes(cross.info)
    
    return(cross.info)
})

# setSamples -------------------------------------------------------------------
#' Set samples by index or sample ID.
#' 
#' @template param-CrossInfo
#' @param samples Integer vector of sample indices or character vector of 
#' sample IDs.
#' 
#' @return Input \code{CrossInfo} object with the given samples.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname setSamples-methods
setGeneric('setSamples', function(cross.info, samples) { 
    standardGeneric('setSamples') })

# CrossInfo::setSamples --------------------------------------------------------
#' @aliases setSamples,CrossInfo-method
#' @export
#' @rdname setSamples-methods
setMethod('setSamples', signature='CrossInfo', definition = 
    function(cross.info, samples) { 
                  
    stopif( missing(samples) )
    
    indices <- getIndices(samples)
        
    if ( is.character(samples) ) {
      
        ids <- samples
        
    } else if ( is.integer(samples) ) {
      
        ids <- NULL
      
        if ( any(samples != indices) ) {
            stop("integer sample vector must contain sample indices")
        }
      
    } else {
      
        stop("sample vector must be of type integer or character")
    }
    
    cross.info@samples <- data.frame(sample.index=indices, 
        stringsAsFactors=FALSE)
    
    if ( ! is.null(ids) ) {
        cross.info@samples$sample.id <- ids
        cross.info@samples$sample.name <- make.names(ids)
    }
    
    validateSamples(cross.info)
    
    return(cross.info)
})

# setSequences -----------------------------------------------------------------
#' Set sequences.
#' 
#' @template param-CrossInfo
#' @param sequences Vector of sequence labels.
#' 
#' @return Input \code{CrossInfo} object with the given sequences.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname setSequences-methods
setGeneric('setSequences', function(cross.info, sequences) 
    { standardGeneric('setSequences') })

# CrossInfo::setSequences ----------------------------------------------------
#' @aliases setSequences,CrossInfo-method
#' @export
#' @rdname setSequences-methods
setMethod('setSequences', signature='CrossInfo', 
    definition = function(cross.info, sequences) { 
              
    norm.seqs <- normSeq(sequences)
    fmt.seqs <- formatSeq(norm.seqs)
              
    cross.info@seq <- structure(norm.seqs, names=fmt.seqs)
              
    validateSequences(cross.info)
              
    return(cross.info)
})

# setStrainIndices -------------------------------------------------------------
#' Set strain indices.
#' 
#' @template param-CrossInfo
#' @param samples Vector of samples to assign strain indices.
#' @param strains Vector of strain indices.
#' 
#' @return Input \code{CrossInfo} object with the given strain indices.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname setStrainIndices-methods
setGeneric('setStrainIndices', function(cross.info, samples=NULL, strains=NULL) { 
    standardGeneric('setStrainIndices') })

# CrossInfo::setStrainIndices --------------------------------------------------
#' @aliases setStrainIndices,CrossInfo-method
#' @export
#' @rdname setStrainIndices-methods
setMethod('setStrainIndices', signature='CrossInfo', definition = function(cross.info, 
    samples=NULL, strains=NULL) {
    
    sample.indices <- getSampleIndices(cross.info, samples)
    
    if ( length(strains) != length(sample.indices) ) {
        stop("cannot assign ", length(strains), " strain indices to ", 
            length(sample.indices)," samples")
    }
    
    if ( ! hasStrainIndices(cross.info) ) {
        cross.info@samples$strain.index <- NA
    }
    
    cross.info@samples$strain.index[sample.indices] <- strains
    
    validateSamples(cross.info)
    
    return(cross.info)
})

# setTetradIndices -------------------------------------------------------------
#' Set sample tetrad indices.
#' 
#' @template param-CrossInfo
#' @param samples Vector of samples for which a tetrad index should be set.
#' @param tetrads Vector of tetrad indices.
#' 
#' @return Input \code{CrossInfo} object with the given tetrad indices.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname setTetradIndices-methods
setGeneric('setTetradIndices', function(cross.info, samples=NULL, 
    tetrads=NULL) { standardGeneric('setTetradIndices') })

# CrossInfo::setTetradIndices --------------------------------------------------
#' @aliases setTetradIndices,CrossInfo-method
#' @export
#' @rdname setTetradIndices-methods
setMethod('setTetradIndices', signature='CrossInfo', definition = 
    function(cross.info, samples=NULL, tetrads=NULL) { 
    
    sample.indices <- getSampleIndices(cross.info, samples)
    
    if ( length(tetrads) != length(sample.indices) ) {
        stop("cannot assign ", length(tetrads), " tetrad indices to ", 
            length(sample.indices)," samples")
    }
    
    if ( ! hasTetradIndices(cross.info) ) {
        cross.info@samples$tetrad.index <- NA
    }
    
    cross.info@samples$tetrad.index[sample.indices] <- tetrads
    
    validateSamples(cross.info)
    
    return(cross.info)    
})

# validateAlleles --------------------------------------------------------------
#' Validate allele information.
#' 
#' @template param-CrossInfo
#' 
#' @return TRUE if alleles are valid; otherwise, returns first error.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname validateAlleles-methods
setGeneric('validateAlleles', function(cross.info) { 
    standardGeneric('validateAlleles') })

# CrossInfo::validateAlleles ---------------------------------------------------
#' @aliases validateAlleles,CrossInfo-method
#' @export
#' @rdname validateAlleles-methods
setMethod('validateAlleles', signature='CrossInfo', 
    definition = function(cross.info) { 
      
    stopifnot( is.character(cross.info@alleles) )
      
    if ( anyNA(cross.info@alleles) ) {
        stop("incomplete allele info")
    }  
    
    dup.geno <- cross.info@alleles[ duplicated(cross.info@alleles) ]
    if ( length(dup.geno) > 0 ) {
        stop("duplicate alleles - '", toString(dup.geno), "'")
    }
    
    err.geno <- cross.info@alleles[ ! grepl('^[[:alnum:]]+$', cross.info@alleles) ]
    if ( length(err.geno) > 0 ) {
        stop("invalid allele values - '", toString(err.geno), "'")
    }
    
    return(TRUE)
})

# validateMarkers --------------------------------------------------------------
#' Validate marker information.
#' 
#' @template param-CrossInfo
#' 
#' @return TRUE if marker information is valid; otherwise, returns first error.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname validateMarkers-methods
setGeneric('validateMarkers', function(cross.info) { 
    standardGeneric('validateMarkers') })

# CrossInfo::validateMarkers ---------------------------------------------------
#' @aliases validateMarkers,CrossInfo-method
#' @export
#' @rdname validateMarkers-methods
setMethod('validateMarkers', signature='CrossInfo', 
    definition = function(cross.info) { 
    
    stopifnot( is.data.frame(cross.info@markers) )
    
    if ( anyNA(cross.info@markers) ) {
        stop("incomplete marker info")
    }    
    
    headings <- colnames(cross.info@markers)
    
    uh <- headings[ ! headings %in% const$marker.headings ]
    if ( length(uh) > 0 ) {
        stop("unknown marker headings - '", toString(uh), "'")
    }  
    
    dh <- headings[ duplicated(headings) ]
    if ( length(dh) > 0 ) {
        stop("duplicate marker headings - '", toString(dh), "'")
    }   
    
    if ( headings[1] != 'marker' ) {
        stop("marker ID column not found")
    }
    
    if ( nrow(cross.info@markers) == 0 ) {
        ih <- headings[ headings != 'marker' ]
        stop("headings invalid in CrossInfo object with zero markers - '", toString(ih), "'")
    }
    
    marker.ids <- cross.info@markers$marker
    
    dup.mkr <- marker.ids[ duplicated(marker.ids) ]
    if ( length(dup.mkr) > 0 ) {
        stop("duplicate marker IDs - '", toString(dup.mkr), "'")
    }
    
    invalid.ids <- marker.ids[ ! isValidID(marker.ids) ]
    if ( length(invalid.ids) > 0 ) {
        stop("invalid marker IDs - '", toString(invalid.ids), "'")
    }
    
    mkr.names <- make.names(cross.info@markers$marker)
      
    err.names <- rownames(cross.info@markers)[ rownames(cross.info@markers) != mkr.names ]
    if ( length(err.names) > 0 ) {
        stop("invalid marker names - '", toString(err.names), "'")
    }
      
    dup.names <- rownames(cross.info@markers)[ duplicated( rownames(cross.info@markers) ) ]
    if ( length(dup.names) > 0 ) {
        stop("duplicate marker names - '", toString(dup.names), "'")
    }
    
    if ( length(headings) > 1 ) {
        
        if ( length(headings) > 2 || headings[2] != 'seq' ) {
            stop("invalid marker info headings")
        }
        
        indices <- which( ! is.na(cross.info@markers$seq) & 
            ! isNormSeq(cross.info@markers$seq) )
        err.seqs <- unique(cross.info@markers$seq[indices])
        if ( length(err.seqs) > 0 ) {
            stop("invalid marker sequences - '", toString(err.seqs), "'")
        }
    }
    
    return(TRUE)
})

# validatePhenotypes -----------------------------------------------------------
#' Validate phenotype information.
#' 
#' @template param-CrossInfo
#' 
#' @return TRUE if phenotypes are valid; otherwise, returns first error.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname validatePhenotypes-methods
setGeneric('validatePhenotypes', function(cross.info) { 
    standardGeneric('validatePhenotypes') })

# CrossInfo::validatePhenotypes ------------------------------------------------
#' @aliases validatePhenotypes,CrossInfo-method
#' @export
#' @rdname validatePhenotypes-methods
setMethod('validatePhenotypes', signature='CrossInfo', 
    definition = function(cross.info) { 
              
    stopifnot( is.character(cross.info@pheno) )
    
    if ( anyNA(cross.info@pheno) ) {
        stop("incomplete phenotype info")
    }

    dup.pheno <- cross.info@pheno[ duplicated(cross.info@pheno) ]
    if ( length(dup.pheno) > 0 ) {
        stop("duplicate phenotypes - '", toString(dup.pheno), "'")
    }
    
    reserved.ids <- cross.info@pheno[ tolower(cross.info@pheno) %in% 
        const$disallowed.phenotypes ]
    if ( length(reserved.ids) > 0 ) {
        stop("disallowed reserved phenotype IDs - '", toString(reserved.ids), "'")
    }
    
    invalid.ids <- cross.info@pheno[ ! isValidID(cross.info@pheno) ]
    if ( length(invalid.ids) > 0 ) {
        stop("invalid phenotype IDs - '", toString(invalid.ids), "'")
    }
    
    phenames <- make.names(cross.info@pheno)
    
    err.names <- names(cross.info@pheno)[ names(cross.info@pheno) != phenames ]
    if ( length(err.names) > 0 ) {
        stop("invalid phenotype names - '", toString(err.names), "'")
    }
    
    dup.names <- names(cross.info@pheno)[ duplicated( names(cross.info@pheno) ) ]
    if ( length(dup.names) > 0 ) {
        stop("duplicate phenotype names - '", toString(dup.names), "'")
    }
    
    return(TRUE)
})

# validateSamples --------------------------------------------------------------
#' Validate sample information.
#' 
#' @template param-CrossInfo
#' 
#' @return TRUE if sample information is valid; otherwise, returns first error.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname validateSamples-methods
setGeneric('validateSamples', function(cross.info) { 
    standardGeneric('validateSamples') })

# CrossInfo::validateSamples ---------------------------------------------------
#' @aliases validateSamples,CrossInfo-method
#' @export
#' @rdname validateSamples-methods
setMethod('validateSamples', signature='CrossInfo', 
    definition = function(cross.info) {
    
    stopifnot( is.data.frame(cross.info@samples) )
    
    if ( anyNA(cross.info@samples) ) {
        stop("incomplete sample info")
    }
    
    headings <- colnames(cross.info@samples)
    
    uh <- headings[ ! headings %in% const$sample.headings ]
    if ( length(uh) > 0 ) {
        stop("unknown sample headings - '", toString(uh), "'")
    }

    dh <- headings[ duplicated(headings) ]
    if ( length(dh) > 0 ) {
        stop("duplicate sample headings - '", toString(dh), "'")
    }   
    
    if ( headings[1] != 'sample.index' ) {
        stop("sample indices not found")
    }
    
    if ( nrow(cross.info@samples) == 0 ) {
        ih <- headings[ headings != 'sample.index' ]
        stop("headings invalid in CrossInfo object with zero samples - '", toString(ih), "'")
    }
    
    if ( any( cross.info@samples$sample.index != getRowIndices(cross.info@samples) ) ) {
        stop("invalid sample indices")
    }
    
    if ( 'sample.id' %in% headings ) {
        
        sample.ids <- cross.info@samples$sample.id
        sample.runs <- rle(sample.ids)
        if ( anyDuplicated(sample.runs$values) ) {
            stop("non-consecutive identical sample IDs")
        }
        
        invalid.ids <- sample.ids[ ! isValidID(sample.ids) ]
        if ( length(invalid.ids) > 0 ) {
            stop("invalid sample IDs - '", toString(invalid.ids), "'")
        }
        
        sample.names <- make.names(sample.ids)
        err.names <- cross.info@samples$sample.name[ cross.info@samples$sample.name != sample.names ]
        if ( length(err.names) > 0 ) {
            stop("invalid sample names - '", toString(err.names), "'")
        }
    }
    
    if ( 'strain.index' %in% headings ) {
        
        strain.indices <- cross.info@samples$strain.index
        
        if ( length(strain.indices) > 0 ) {
            
            strain.start.valid <- strain.indices[1] == 1
            strain.steps.valid <- all( diff(strain.indices) %in% 0:1 )
            
            if ( ! ( strain.start.valid && strain.steps.valid ) ) {
                stop("invalid strain indices")
            }
            
            if ( 'sample.id' %in% headings ) {
                
                strain.runs <- rle(strain.indices)
                
                if ( ! ( length(strain.runs$lengths) == length(sample.runs$lengths) && 
                    all(strain.runs$lengths == sample.runs$lengths) ) ) {
                    stop("mismatch between sample IDs and strain indices")
                } 
            }
        }
    }

    if ( 'tetrad.index' %in% headings ) {
        
        tetrad.indices <- cross.info@samples$tetrad.index
        
        if ( length(tetrad.indices) > 0 ) {
        
            if ( ! 'strain.index' %in% headings ) {
                strain.indices <- getRowIndices(cross.info@samples)
            }
            
            exemplar.sindices <- sapply(unique(strain.indices), match, strain.indices)
            exemplar.tindices <- tetrad.indices[ exemplar.sindices ]
            
            tetrad.start.valid <- tetrad.indices[1] == 1
            tetrad.steps.valid <- all( diff(tetrad.indices) %in% 0:1 )
            tetrad.sizes.valid <- all( table(exemplar.tindices) <= 4 )
            
            if ( ! ( tetrad.start.valid && tetrad.steps.valid && tetrad.sizes.valid ) ) {
                stop("invalid tetrad indices")
            }
        }
    }
    
    return(TRUE)
})

# validateSequences ----------------------------------------------------------
#' Validate sequence information.
#' 
#' @template param-CrossInfo
#' 
#' @return TRUE if sequences are valid; otherwise, returns first error.
#' 
#' @docType methods
#' @export
#' @keywords internal
#' @rdname validateSequences-methods
setGeneric('validateSequences', function(cross.info) { 
    standardGeneric('validateSequences') })

# CrossInfo::validateSequences -----------------------------------------------
#' @aliases validateSequences,CrossInfo-method
#' @export
#' @rdname validateSequences-methods
setMethod('validateSequences', signature='CrossInfo', 
    definition = function(cross.info) { 
      
    stopifnot( is.character(cross.info@seq) )
      
    if ( anyNA(cross.info@seq) ) {
        stop("incomplete sequences info")
    }  
      
    norm.seqs <- normSeq(cross.info@seq)
    err.ids <- cross.info@seq[ cross.info@seq != norm.seqs ]
    if ( length(err.ids) > 0 ) {
        stop("invalid sequence labels - '", toString(err.ids), "'")
    }
    
    fmt.seqs <- formatSeq(norm.seqs)
    err.names <- names(cross.info@seq)[ names(cross.info@seq) != fmt.seqs ]
    if ( length(err.names) > 0 ) {
        stop("invalid sequence labels - '", toString(err.names), "'")
    }
    
    dup.seqs <- cross.info@seq[ duplicated(cross.info@seq) ]
    if ( length(dup.seqs) > 0 ) {
        stop("duplicate sequence labels - '", toString(dup.seqs), "'")
    }
      
    return(TRUE)
})

# CrossInfo::setValidity -------------------------------------------------------
#' Validate \code{CrossInfo} object.
#' 
#' @param object A \code{CrossInfo} object.
#' 
#' @return TRUE if object is valid; otherwise, a character vector of errors.
#' 
#' @aliases setValidity,CrossInfo-method
#' @docType methods
#' @keywords internal
#' @name setValidity
#' @rdname setValidity-methods
setValidity('CrossInfo', function(object) { 

    errors <- vector('character')
    
    validators <- c(validateAlleles, validateSequences,
        validateMarkers, validatePhenotypes, validateSamples)
    
    for ( validator in validators ) {
        tryCatch({
            validator(object)
        }, error=function(e) {
            errors <- c(errors, e[['message']])
        })
    }
    
    if ( hasMarkerSeqs(object) ) {
        
        orphans <- object@markers$marker[ ! object@markers$seq %in% object@seq ]
        if ( length(orphans) > 0 ) {
            e <- paste0("no sequence found for markers - '", toString(orphans), "'")
            errors <- c(errors, e)
        }
        
        empties <- object@seq[ ! object@seq %in% object@markers$seq ]
        if ( length(empties) > 0 ) {
            e <- paste0("no markers found for sequences - '", toString(empties), "'")
            errors <- c(errors, e)
        }
    }
    
    return( if ( length(errors) == 0 ) {TRUE} else {errors}  )
})

# End of CrossInfo-class.R #####################################################