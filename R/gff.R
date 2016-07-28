# Start of gff.R ###############################################################

# readFeaturesGFF --------------------------------------------------------------
#' Read annotation from GFF file.
#' 
#' @param annofile Input annotation GFF file.
#' 
#' @return A \code{data.frame} of primary genome features, such that each
#' row contains at least the reference sequence (\code{'chr'}), start and
#' end positions (\code{'start'} and \code{'end'} respectively), and \code{'ID'}
#' of a feature. Additional feature annotation is included, if available.
#' 
#' @export
#' @importFrom rtracklayer readGFF
#' @importFrom S4Vectors unstrsplit
#' @include const.R
#' @include seq.R
#' @include util.R
#' @rdname readFeaturesGFF
readFeaturesGFF <- function(annofile) {
    
    stopifnot( isSingleString(annofile) )
    stopifnot( file.exists(annofile) )
    
    # Get info on annotation headings.
    supp.hdgs <- const$anno$supported.headings
    req.hdgs <- const$anno$required.headings
    
    # Read annotation DataFrame from GFF file.
    anno <- rtracklayer::readGFF(annofile, columns=const$anno$columns,
        tags=const$anno$tags)
    
    # Identify and remove unwanted features.
    irrelevant <- anno$type %in% const$anno$irrelevant
    secondary <- lengths(anno$Parent) > 0
    unidentified <- is.na(anno$ID)
    anno <- anno[ ! ( unidentified | secondary | irrelevant ), ]
    
    # Create features data-frame (S3) from annotation DataFrame (S4).
    features <- as.data.frame( matrix( nrow=nrow(anno), ncol=length(supp.hdgs),
        dimnames=list(NULL, supp.hdgs) ), stringsAsFactors=FALSE )
    
    features$chr <- normSeq(anno$seqid)
    features$start <- anno$start
    features$end <- anno$end
    features$strand <- anno$strand
    features$ID <- anno$ID
    features$Alias <- S4Vectors::unstrsplit(anno$Alias, '; ')
    features$type <- as.character(anno$type)
    features$orf_classification <- anno$orf_classification
    features$source <- as.character(anno$source)
    features$dbxref <- anno$dbxref
    features$Ontology_term <- S4Vectors::unstrsplit(anno$Ontology_term, ' ')
    features$Note <- S4Vectors::unstrsplit(anno$Note, '...')
    
    # Remove empty columns.
    nonempty <- sapply( getColIndices(features),
        function(i) ! allNA(features[, i]) )
    features <- features[, nonempty]
    
    unfound <- req.hdgs[ ! req.hdgs %in% colnames(features) ]
    if ( length(unfound) > 0 ) {
        stop("required annotation headings not found - '", unfound, "'")
    }
    
    # Check feature positions are within range of reference sequence.
    seqtab <- getSeqTable()
    for ( anno.seq in unique(features$chr) ) {
        seq.pos <- unname( unlist(
            features[ features$chr == anno.seq, c('start', 'end') ] ) )
        seq.length <- seqtab[ seqtab$seqids == anno.seq, 'seqlengths' ]
        if ( any( seq.pos < 1 | seq.pos > seq.length ) ) {
            stop("features out of range of reference sequence")
        }
    }
    
    return(features)
}

# End of gff.R #################################################################