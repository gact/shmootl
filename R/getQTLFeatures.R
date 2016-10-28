# Start of getQTLFeatures.R ####################################################

# getQTLFeatures ---------------------------------------------------------------
#' Get annotation of features within QTL intervals.
#'
#' @param qtl.intervals A \code{qtlintervals} object.
#' @param features A \code{data.frame} of primary genome features, as returned
#' by \code{\link{readFeaturesGFF}}. For each feature, this \code{data.frame}
#' must include the reference sequence (\code{'chr'}), 1-offset base-pair
#' positions of the start and end of the feature (\code{'start'} and
#' \code{'end'} respectively), and the feature \code{'ID'}.
#' 
#' @return An object of class \code{qtlfeatures}, which is essentially
#' a list of \code{data.frame} objects, each containing the set of genome
#' features in a given QTL interval.
#' 
#' @export
#' @family QTL functions
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom IRanges restrict
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom methods slotNames
#' @rdname getQTLFeatures
getQTLFeatures <- function(qtl.intervals, features) {
    
    stopifnot( hasPhysicalPositions(qtl.intervals) )
    stopifnot( is.data.frame(features) )
    
    # Get info on annotation headings.
    supp.hdgs <- const$anno$supported.headings
    req.hdgs <- const$anno$required.headings
    
    # Identify and remove unwanted features.
    irrelevant <- features$type %in% const$anno$irrelevant
    unidentified <- is.na(features$ID)
    features <- features[ ! ( unidentified | irrelevant ), ]
    
    # Remove empty feature columns.
    nonempty <- sapply( getColIndices(features),
        function(i) ! allNA(features[, i]) )
    features <- features[, nonempty]
    
    unfound <- req.hdgs[ ! req.hdgs %in% colnames(features) ]
    if ( length(unfound) > 0 ) {
        stop("required annotation headings not found - '", toString(unfound), "'")
    }
    
    unknown <- colnames(features)[ ! colnames(features) %in% supp.hdgs ]
    if ( length(unknown) > 0 ) {
        stop("unsupported annotation headings - '", toString(unknown), "'")
    }
    
    # Init QTL features.
    qtl.features <- vector('list', length(qtl.intervals))
    class(qtl.features) <- c('qtlfeatures', 'list')
    names(qtl.features) <- names(qtl.intervals)
    
    # If there are QTL intervals, get annotated features for each QTL interval.
    if ( length(qtl.intervals) > 0 ) {
        
        # Get Seqinfo object for current genome.
        seqinfo <- getSeqinfo()
        
        # Get ranges of annotated features.
        feature.ranges <- GenomicRanges::makeGRangesFromDataFrame(features,
            seqnames.field='chr', start.field='start', end.field='end',
            strand.field='strand', seqinfo=seqinfo)
        
        # Get QTL interval ranges.
        qtl.seqs <- sapply( qtl.intervals, function(x) unique(x[, 'chr']) )
        qtl.starts <- sapply(qtl.intervals, function(x) x[1, 'pos (bp)'])
        qtl.ends <- sapply(qtl.intervals, function(x) x[3, 'pos (bp)'])
        qtl.iranges <- IRanges::IRanges(start=qtl.starts, end=qtl.ends)
        
        # Restrict QTL interval ranges to within genome sequences.
        seq.names <- GenomeInfoDb::seqnames(seqinfo)
        seq.starts <- structure(rep(1, length(seqinfo)), names=seq.names)
        seq.ends <- GenomeInfoDb::seqlengths(seqinfo)
        qtl.iranges <- IRanges::restrict(qtl.iranges, start=seq.starts[qtl.seqs],
            end=seq.ends[qtl.seqs], keep.all.ranges=TRUE)
        
        # Get QTL genomic ranges.
        qtl.ranges <- GenomicRanges::GRanges(seqnames=qtl.seqs,
            ranges=qtl.iranges, seqinfo=seqinfo)
        
        # Find overlaps between QTL intervals and genome features.
        qtl.overlaps <- GenomicRanges::findOverlaps(qtl.ranges,
            feature.ranges, ignore.strand=TRUE)
        
        # Assign QTL features overlapping each QTL interval.
        for ( i in seq_along(qtl.intervals) ) {
            
            if ( 'queryHits' %in% methods::slotNames(qtl.overlaps) ) {
                row.indices <- qtl.overlaps@subjectHits[ qtl.overlaps@queryHits == i ]
            } else {
                row.indices <- qtl.overlaps@to[ qtl.overlaps@from == i ]
            }
            
            qtl.features[[i]] <- features[row.indices, ]
            
            rownames(qtl.features[[i]]) <- NULL
        }
    }
    
    return(qtl.features)
}

# End of getQTLFeatures.R ######################################################