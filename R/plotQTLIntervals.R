# Start of plotQTLIntervals.R ##################################################

# plotQTLIntervals -------------------------------------------------------------
#' Plot QTL intervals.
#'
#' @usage 
#' ## Generic method.
#' plotQTLIntervals(x, ...)
#'
#' @param x A \code{scanone} or \code{qtlintervals} object.
#' 
#' @param ... Further arguments (see below).
#' @param lod.profile For a \code{qtlintervals} object, the \code{scanone}
#' object from which the QTL intervals were obtained.
#' @param chr Vector indicating which sequences to plot. If no sequences are
#' specified, all sequences are plotted.
#' @param drop For a \code{scanone} object, LOD units that the LOD profile must
#' drop to form a QTL interval.
#' @param expandtomarkers For a \code{scanone} object, expand the LOD interval
#' to the nearest flanking markers, or to the respective terminal loci.
#' @param threshold For a \code{scanone} object, this indicates 
#' the LOD significance threshold for identifying QTL intervals.
#' @param alpha For a \code{scanone} object, this contains
#' the significance level of the given threshold.
#' @template param-lodcolumn
#' @param phenotype Name of the phenotype, to be used in plot title.
#' 
#' @export
#' @importFrom graphics abline
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics strwidth
#' @importFrom graphics text
#' @importFrom graphics yinch
#' @importFrom grDevices graphics.off
#' @importFrom utils tail
#' @rdname plotQTLIntervals
plotQTLIntervals <- function(x, ...) {
    UseMethod('plotQTLIntervals', x)
}

# plotQTLIntervals.qtlintervals ------------------------------------------------
#' @export
#' @rdname plotQTLIntervals
plotQTLIntervals.qtlintervals <- function(x, lod.profile, chr=NULL,
    lodcolumn=NULL, phenotype=NULL, ...) {
    
    stopifnot( 'scanone' %in% class(lod.profile) )
    
    # Get LOD profile column indices.
    seqcol.index <- getSeqColIndex(lod.profile)
    poscol.index <- getPosColIndex(lod.profile)
    lodcol.index <- getDatColIndices(lod.profile, datcolumns=lodcolumn)
    
    if ( length(lodcol.index) > 1 ) {
        stop("cannot plot QTL intervals for multiple LOD columns - please choose one")
    } else if ( length(lodcol.index) == 0 ) {
        stop("no LOD column found")
    }
    
    # Get LOD threshold.
    threshold <- attr(x, 'threshold')
    
    # Get alpha value corresponding to LOD threshold.
    alpha <- attr(x, 'alpha')
    
    # Ensure LOD profile has normalised sequence IDs.
    lod.profile <- normSeq(lod.profile)
    
    # Get specified sequences.
    lodprof.seqs <- unique( pullLocusSeq(lod.profile) )
    chr <- lodprof.seqs <- subsetBySeq(lodprof.seqs, chr)
    stopifnot( length(chr) > 0 )
    
    # Subset QTL intervals by specified sequences.
    x <- subsetBySeq(x, chr)
    
    # Get sequences corresponding to QTL intervals.
    interval.seqs <- sapply(x, function(obj) unique(obj[, 'chr']))
    
    # Subset LOD profile by specified sequences.
    lod.profile <- subsetBySeq(lod.profile, chr)
    
    # Get sequences and positions of individual LOD profile loci.
    locus.seqs <- pullLocusSeq(lod.profile)
    
    # If phenotype specified, add it to plot title..
    if ( ! is.null(phenotype) ) {
        
        stopifnot( isSingleString(phenotype) )
        plot.title <- paste0("LOD Plot (", phenotype, ")")
        
    } else {
        
        # ..otherwise get the name of the given LOD column.
        profile.phename <- colnames(lod.profile)[lodcol.index]
        
        # If LOD column name is a phenotype, add it to plot title..
        if ( profile.phename != 'lod' ) {
            
            plot.title <- paste0("LOD Plot (", profile.phename, ")")
            
        } else { # ..otherwise set default plot title.
            
            plot.title <- "LOD Plot"
        }
    }
    
    # For single-sequence plots, add sequence label to plot title.
    if ( length(chr) == 1 ) {
        plot.title <- paste('Chromosome', chr, plot.title)
    }
    
    # Include markers if plotting one chromosome with fewer than 100 markers.
    incl.markers = length(chr) == 1 && nrow(lod.profile) < 100
    
    # Set maximum y-value, ensure greater than or equal to one.
    max.lod <- max( lod.profile[, lodcol.index],
        threshold, 1.0, na.rm=TRUE )
    
    # Set inter-sequence gap width (value from R/qtl docs).
    inter.seq.gap <- 25 # cM
    
    # Set width of gaps between sequences.
    gap.width <- ifelse(length(chr) > 1, inter.seq.gap, 0)
    half.gap <- 0.5 * gap.width
    
    # Get last sequence.
    last.seq <- tail(chr, 1)
    
    # Init gap position vector.
    gap.positions <- numeric()
    
    # Init sequence offset vector.
    seq.offsets <- numeric()
    
    # Init cumulative plot width.
    cum.plot.width <- 0.0
    
    # Get horizontal positions of gaps between the
    # sequences, for all but the last sequence.   
    for ( plot.seq in chr ) {
        
        # Get indices of cross result for this sequence.
        row.indices <- which(locus.seqs == plot.seq)
        
        # Get index of last locus for this sequence.
        last.locus <- tail(row.indices, 1)
        
        # Get sequence length from last locus.
        seq.length <- lod.profile[last.locus, poscol.index]
        
        # Get sequence offset from cumulative plot width.
        seq.offsets <- c(seq.offsets, cum.plot.width)
        
        # Add sequence length to cumulative plot width.
        cum.plot.width <- cum.plot.width + seq.length
        
        # If this isn't the last sequence, add gap marker position.
        if ( plot.seq != last.seq ) {
            gap.positions <- c(gap.positions, cum.plot.width + half.gap)
        }
        
        # Add gap width to cumulative plot width.
        cum.plot.width <- cum.plot.width + gap.width
    }
    
    # Subtract gap width from final cumulative plot width.
    cum.plot.width <- cum.plot.width - gap.width
    
    # Set names of sequence offsets.
    names(seq.offsets) <- chr
    
    # Ensure graphics parameters will be reset.
    on.exit({ par(op) })
    
    # Set font settings.
    op <- par(cex=1.0, family='sans')
    
    # Set xlim (plot space horizontal range); values taken from R/qtl source code.
    xlim <- c(-half.gap, cum.plot.width + half.gap)
    
    # Set y-axis range from range of LOD values.
    ylim <- c( 0.0, ceiling(max.lod + 0.5) )
    
    # Draw LOD plot.
    plot(lod.profile, xlim=xlim, ylim=ylim, chr=chr, main=plot.title,
         ylab='LOD', show.marker.names=FALSE, incl.markers=incl.markers)
    
    # Get user space parameters.
    usr <- par('usr')
    
    # Draw dashed lines in gap positions between chromosomes.
    for ( gap.pos in gap.positions ) {
        segments(gap.pos, usr[3], gap.pos, usr[4], col='gray', lwd=0.5, lty='dashed')
    }
    
    # Init interval regions matrix.
    interval.regions <- matrix(nrow=length(x), ncol=2)
    
    # Set vertical offset of 1.5-LOD interval from LOD peak.
    y.offset <- yinch(0.25)
    
    # Display each QTL interval, remember the plot regions it occupies.
    for ( i in getIndices(x) ) {
        
        qtl.interval <- x[[i]]

        # Get the horizontal offset for this sequence (i.e. cumulative 
        # length of previous sequences and gaps between.).
        x.offset <- seq.offsets[ interval.seqs[[i]] ]
        
        # Get horizontal positions of QTL interval parts.
        xpos <- qtl.interval[, 'pos'] + x.offset
        
        # Get peak LOD value for this QTL interval.
        peak.lod <- qtl.interval[2, 'lod']
        
        # Get interval line vertical position from peak LOD and vertical offset.
        iline <- peak.lod + y.offset
        
        # Get vertical positions of QTL interval parts.
        ypos <- c(iline - 0.5 * y.offset, iline, iline + 0.5 * y.offset)
        
        # Get QTL interval line segment endpoints.
        x0 <- c(xpos[1], xpos[1], xpos[3])
        x1 <- c(xpos[1], xpos[3], xpos[3])
        y0 <- c(ypos[1], ypos[2], ypos[1])
        y1 <- c(ypos[3], ypos[2], ypos[3])
        
        # Draw QTL interval line segments.
        segments(x0, y0, x1, y1, col='black', lwd=0.5, lty='solid')
        
        # If this LOD plot is of one chromosome, draw point at LOD peak.
        if ( length(chr) == 1 ) {
            points(xpos[2], iline, col='black', lwd=0.5, pch=20)
        }
        
        # Add this region to the matrix of interval regions.
        interval.regions[i, ] <- c(xpos[1], xpos[3])
    }
    
    if ( ! is.null(threshold) ) {
        
        # Draw horizontal red dashed line to indicate LOD threshold.
        abline(threshold, 0, col='gray32', lwd=0.5, lty='dashed')
        
        if ( ! is.null(alpha) ) {
        
            # Set threshold label text.
            thresh.label.text <- bquote( alpha ~ '=' ~ .(alpha) )
            
            # Set threshold label width.
            thresh.label.width <- strwidth(thresh.label.text)
            
            # Set threshold label offset from rightmost available position.
            thresh.label.offset <- 0.025 * cum.plot.width
            
            # Set threshold label default position.
            default.position <- cum.plot.width - thresh.label.offset - thresh.label.width
            
            # Set threshold label position from default.
            thresh.label.pos <- default.position
            
            # If interval regions are present, try to adjust threshold label to avoid these.
            if ( nrow(interval.regions) > 0 ) {
                
                # Check label against each interval in reverse order. Shift label left for 
                # each coinciding interval. Stop if label does not coincide with any interval.
                for ( r in nrow(interval.regions):1 ) {
                    
                    region <- interval.regions[r,]
                    
                    if ( (thresh.label.pos + thresh.label.width) >= region[1] && region[2] >= thresh.label.pos ) {
                        thresh.label.pos <- region[1] - thresh.label.offset - thresh.label.width
                    } else {
                        break
                    }
                }
                
                # If, after shifting threshold label left to avoid QTL intervals, we have
                # moved the label off the plot, just put it back in the default position.
                if (thresh.label.pos < 0.0) {
                    thresh.label.pos <- default.position
                }
            }
            
            # Draw LOD significance threshold label.
            text(thresh.label.pos, threshold + yinch(0.025), thresh.label.text, 
                col='gray32', adj=c(0, 0), cex=0.8)
        }
    }
    
    return( invisible() )
}

# plotQTLIntervals.scanone -----------------------------------------------------
#' @export
#' @rdname plotQTLIntervals
plotQTLIntervals.scanone <- function(x, chr=NULL, drop=1.5, expandtomarkers=FALSE,
    threshold=NULL, alpha=NULL, lodcolumn=NULL, phenotype=NULL, ...) {
    
    intervals <- getQTLIntervals(x, chr=chr, drop=drop, threshold=threshold,
        alpha=alpha, expandtomarkers=expandtomarkers, lodcolumn=lodcolumn)
    
    plotQTLIntervals(intervals, x, chr=chr, phenotype=phenotype)
    
    return( invisible() )
}

# End of plotQTLIntervals.R ####################################################