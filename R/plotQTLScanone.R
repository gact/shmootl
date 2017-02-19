# Start of plotQTLScanone.R ####################################################

# plotQTLScanone ---------------------------------------------------------------
#' Plot \pkg{R/qtl} \code{scanone} result.
#' 
#' Plotting function based on \pkg{R/qtl} \code{plot.scanone} (Broman
#' \emph{et al.} 2003) and \pkg{qqman} \code{manhattan} (Turner 2014).
#' Plots a \code{scanone} result as a LOD curve (with \code{type} set
#' to \code{'l'}) or as a Manhattan plot (with \code{type} set to
#' \code{'p'}). If no plot type is specified, this is set automatically
#' based on marker density.
#' 
#' @param x An \pkg{R/qtl} \code{scanone} object.
#' @param chr Vector indicating which sequences to plot. If no sequences are
#' specified, all are plotted.
#' @template param-lodcolumn
#' @param qtl.intervals A \code{qtlintervals} object created from the same data
#' as the \code{scanone} object. These intervals are added to the plot, unless
#' the median interval size is too small to be distinguishable. Threshold
#' information is taken from this object, if available.
#' @param threshold A single \code{numeric} LOD significance threshold, or an
#' object (e.g. \code{summary.scanoneperm}) containing one such threshold and
#' its associated significance level.
#' @param col Analogous to the standard \code{'col'} plotting parameter. This is
#' recycled to match the number of sequences being plotted. As in the package
#' \pkg{qqman}, this defaults to two alternating monochrome shades.
#' @param gap Gap (in cM) between chromosomes in multi-chromosome plots.
#' @param phenotype Name of the phenotype, to be shown in plot information.
#' @param type Type of plot. Set to \code{'l'} to output a stanard LOD curve, or
#' to \code{'p'} for a Manhattan plot of individual LOD scores.
#' @param ... Unused arguments.
#' 
#' @template ref-broman-2003
#' @template ref-turner-2014
#' 
#' @export
#' @family plot functions
#' @family QTL functions
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics box
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics strwidth
#' @importFrom graphics text
#' @importFrom graphics xinch
#' @importFrom graphics yinch
#' @importFrom grDevices graphics.off
#' @importFrom stats median
#' @importFrom utils head
#' @importFrom utils tail
#' @rdname plotQTLScanone
plotQTLScanone <- function(x, chr=NULL, lodcolumn=NULL, qtl.intervals=NULL,
    threshold=NULL, col=c('gray10', 'gray60'), gap=25, phenotype=NULL,
    type=NULL, ...) {
    
    opar <- graphics::par(no.readonly=TRUE)
    on.exit( graphics::par(opar) )
    
    stopifnot( 'scanone' %in% class(x) )
    stopifnot( isSinglePositiveWholeNumber(gap) )
    stopifnot( emptyArgs(...) )
    
    # Assume no threshold info.
    tinfo <- NULL
    
    # Get threshold info, if specified directly.
    if ( ! is.null(threshold) ) {
        tinfo <- getScanoneThresholdInfo(threshold)
    }
    
    # Get threshold info from QTL intervals, if specified.
    if ( ! is.null(qtl.intervals) ) {
        
        stopifnot( 'qtlintervals' %in% class(qtl.intervals) )
        qinfo <- getScanoneThresholdInfo(qtl.intervals)
        
        if ( is.null(threshold) ) {
            tinfo <- qinfo
        } else if ( ! identical(qinfo, tinfo) ) {
            stop("threshold information does not match that of QTL intervals")
        }
    }
    
    if ( ! is.null(type) && ! type %in% c('l', 'p') ) {
        stop("unsupported plot type - '", toString(type), "'")
    }
    
    # Get scanone result column indices.
    poscol.index <- getPosColIndex(x)
    lodcol.index <- getLodColIndex(x, lodcolumn=lodcolumn)
    
    # Ensure scanone result has normalised sequence IDs.
    x <- normSeq(x)
    
    # Get specified sequences.
    x.seqs <- pull.chr(x)
    chr <- x.seqs <- subsetBySeq(x.seqs, chr)
    stopifnot( length(chr) > 0 )
    
    # Subset scanone result by specified sequences.
    x <- subsetBySeq(x, chr)
    
    # Get sequences and positions of individual scanone result loci.
    locus.seqs <- pullLocusSeq(x)
    
    # Replace any NA values with zero.
    # TODO: consider removing these
    x[ is.na(x[, lodcol.index]), lodcol.index] <- 0
    
    # Get list containing scanone row indices for each sequence.
    seq.indices <- lapply(chr, function(s) which(locus.seqs == s) )
    names(seq.indices) <- chr
    
    # Set maximum y-value, ensure greater than or equal to one.
    max.lod <- max( x[, lodcol.index], tinfo[['threshold']], 1.0, na.rm=TRUE )
    
    # Assemble sequence plotting info ------------------------------------------
    
    # Set width of gaps between sequences.
    gap.width <- ifelse(length(chr) > 1, gap, 0)
    
    # Init cumulative plot width.
    cum.plot.width <- 0.0
    
    # Init sequence plotting info.
    seq.par <- matrix( NA_real_, nrow=length(chr), ncol=3,
        dimnames=list(chr, c('offset', 'midpoint', 'length') ) )
    
    for ( i in seq_along(chr) ) {
        
        if ( length(chr) > 1 ) {
            seq.start <- x[head(seq.indices[[i]], 1), poscol.index]
        } else {
            seq.start <- 0
        }
        
        seq.end <- x[tail(seq.indices[[i]], 1), poscol.index]
        
        seq.par[i, 'offset'] <- cum.plot.width
        
        seq.par[i, 'length'] <- seq.end - seq.start
        
        seq.par[i, 'midpoint'] <- seq.par[i, 'offset'] + (0.5 * seq.par[i, 'length'])
        
        cum.plot.width <- cum.plot.width + seq.par[i, 'length']
        
        if ( i < length(chr) ) {
            cum.plot.width <- cum.plot.width + gap.width
        }
    }
    
    # Assemble general plot info -----------------------------------------------
    
    plot.info <- mapping()
    
    # Get phenotype from argument or LOD column, if possible.
    if ( ! is.null(phenotype) ) {
        stopifnot( isSingleString(phenotype) )
        plot.info['Phenotype'] <- phenotype
    } else {
        profile.phename <- colnames(x)[lodcol.index]
        if ( profile.phename != 'lod' ) {
            plot.info['Phenotype'] <- profile.phename
        }
    }
    
    # For single-sequence plots, add sequence label to plot info.
    if ( length(chr) == 1 ) {
        plot.info['Chromosome'] <- chr
    }
    
    # --------------------------------------------------------------------------
    
    # Set top margin line counts, accounting for any plot info lines.
    adj.top.mar <- length(plot.info) + 4 # NB: minimum 4 margin lines
    graphics::par( mar=(c(5, 4, adj.top.mar, 2) + 0.1) )
    
    # Set plot info margin indices.
    plot.info.indices <- rev( seq_along(plot.info) )
    plot.title.index <- adj.top.mar - 2
    
    # Set x-axis plotting parameters; values taken from R/qtl.
    if ( length(chr) > 1 ) {
        xlim <- c(-(0.5 * gap.width), cum.plot.width + (0.5 * gap.width))
        xlab <- 'Chromosome'
        xaxt <- 'n'
    } else {
        xlim <- c(0, cum.plot.width)
        xlab <- 'Map position (cM)'
        xaxt <- 's'
    }
    
    # Set y-axis plotting parameters.
    headspace <- 0.2
    ylim <- c( 0.0, max.lod / (1.0 - headspace) )
    ylab <- 'LOD'
    
    # Set fixed plotting arguments.
    fixed.args <- list(family='sans', xaxs='i', xpd=FALSE, yaxs='i')
    
    # Set default plotting arguments; defaults taken from R/qtl and R/qqman.
    default.args <- list(bg='white', cex=1, las=1, lty=1,
        xaxt=xaxt, xlab=xlab, xlim=xlim, ylab=ylab, ylim=ylim)
    
    # Set args from fixed and default values.
    # TODO: allow user to set some arguments
    args <- c(fixed.args, default.args)
    
    # Set palette for sequences.
    col <- rep_len(col, length(chr))
    
    # Start new plot.
    do.call('plot', c(NA, args))
    
    # If plot type not specified, choose based on marker density.
    if ( is.null(type) ) {
        
        num.markers <- length( locus.seqs[ isMarkerID(locus.seqs) ] )
        
        if ( num.markers > 10000 ) { # TODO: improve
            type <- 'p'
        } else {
            type <- 'l'
        }
    }
    
    # Plot scanone result for each sequence.
    for ( i in seq_along(chr) ) {
        
        seq.x <- x[seq.indices[[i]], poscol.index] + seq.par[i, 'offset']
        seq.y <- x[seq.indices[[i]], lodcol.index]
        
        if (type == 'l') {
            lines(seq.x, seq.y, col=col[i], lwd=2, lty='solid')
        } else if (type == 'p') {
            points(seq.x, seq.y, col=col[i], lwd=0.5, pch=20)
        }
    }
    
    # If QTL intervals available, add to plot if they will be visually distinguishable.
    if ( ! is.null(qtl.intervals) ) {
        
        # Subset QTL intervals by specified sequences.
        qtl.intervals <- subsetBySeq(qtl.intervals, chr)
        
        if ( length(qtl.intervals) > 0 ) {
            
            # Get size of QTL interval in centiMorgans.
            interval.widths <- sapply(qtl.intervals, function(qtl.interval)
                diff(qtl.interval[c(1,3), 'pos (cM)']) )
            
            # Try to get size of bullet symbol (pch=20) in user-space..
            bullet.width <- tryCatch({
                result <- suppressWarnings( strwidth('\u2022',
                    cex=args$cex, family ='sans') )
            }, error=function(e) { # ..otherwise get default bullet width.
                bullet.inches <- 0.06944444 # default bullet width in inches (family='sans', cex=1)
                plot.inches <- graphics::par('pin')[1] # plot width in inches
                plot.space <- diff(xlim) # plot width in user-space
                result <- (bullet.inches / plot.inches) * plot.space
            })
            
            # Plot QTL intervals if plotting LOD curves and typical
            # QTL intervals would be visually distinguishable.
            if ( type == 'l' && median(interval.widths) >= bullet.width ) {
                
                # Get sequences corresponding to QTL intervals.
                interval.seqs <- sapply(qtl.intervals,
                    function(obj) unique(obj[, 'chr']))
                
                # Set vertical offset of 1.5-LOD interval from
                # LOD peak in terms of the overall plot height.
                y.offset <- 0.1 * diff(ylim)
                
                # Display each QTL interval, remember the plot regions it occupies.
                for ( i in seq_along(qtl.intervals) ) {
                    
                    qtl.interval <- qtl.intervals[[i]]
                    
                    # Get the horizontal offset for this sequence (i.e. cumulative
                    # length of previous sequences and gaps between.).
                    interval.offset <- seq.par[interval.seqs[[i]], 'offset']
                    
                    # Get horizontal positions of QTL interval parts.
                    interval.xpos <- qtl.interval[, 'pos (cM)'] + interval.offset
                    
                    # Get peak LOD value for this QTL interval.
                    peak.lod <- qtl.interval[2, 'lod']
                    
                    # Get interval line vertical position from peak LOD and vertical offset.
                    iline <- peak.lod + y.offset
                    
                    # Get vertical positions of QTL interval parts.
                    interval.ypos <- c(
                        iline - 0.25 * y.offset,
                        iline,
                        iline + 0.25 * y.offset
                    )
                    
                    # Get QTL interval line segment endpoints.
                    x0 <- interval.xpos[c(1,1,3)]
                    x1 <- interval.xpos[c(1,3,3)]
                    y0 <- interval.ypos[c(1,2,1)]
                    y1 <- interval.ypos[c(3,2,3)]
                    
                    # Draw QTL interval line segments.
                    segments(x0, y0, x1, y1, col='black', lwd=0.5, lty='solid')
                    
                    # Draw point at position of LOD peak.
                    points(interval.xpos[2], iline, col='black', lwd=0.5, pch=20)
                }
            }
        }
    }
    
    # Plot LOD threshold if available.
    if ( ! is.null(tinfo[['threshold']]) ) {
        
        # Draw horizontal red dashed line to indicate LOD threshold.
        abline(tinfo[['threshold']], 0, col='red', lwd=1.0, lty='dotted')
        
        # If significance level or false-discovery rate available, add to threshold line.
        if ( ! is.null(tinfo[['alpha']]) || ! is.null(tinfo[['fdr']]) ) {
            
            # Set threshold label text.
            if ( ! is.null(tinfo[['alpha']]) ) {
                thresh.label.text <- bquote( bold(alpha ~ '=' ~ .(tinfo[['alpha']])) )
            } else { # ! is.null(fdr)
                thresh.label.text <- paste0('FDR=', tinfo[['fdr']])
            }
            
            # Set threshold label width.
            thresh.label.width <- strwidth(thresh.label.text)
            
            # Set threshold label position.
            thresh.label.pos <- cum.plot.width - thresh.label.width - xinch(0.025)
            
            # Draw LOD significance threshold label.
            text( thresh.label.pos, tinfo[['threshold']] + yinch(0.025),
                thresh.label.text, col='red', adj=c(0, 0), cex=0.8)
        }
    }
    
    # If multiple sequences, add sequence ticks and labels.
    if ( length(chr) > 1 ) {
        
        # Adjust size of sequence labels to ensure that they don't overlap.
        seq.label.size <- max( strwidth(chr) )
        seq.plot.size <- min( seq.par[, 'length'] )
        if ( seq.plot.size < seq.label.size ) {
            cex.axis <- seq.plot.size / seq.label.size
        } else {
            cex.axis <- 1
        }
        
        # Plot x-axis with given ticks and labels.
        for ( i in seq_along(chr) ) {
            axis(side=1, at=seq.par[i, 'midpoint'], labels=chr[i], cex.axis=cex.axis)
        }
    }
    
    # Plot box around graph.
    box(lwd=3)
    
    # Write plot title.
    graphics::title(main='Scanone', line=plot.title.index, cex=1.2, col='black', family='sans')
    
    # If any plot info, add to top margin.
    if ( length(plot.info) > 0 ) {
        
        plot.info.lines <- sapply( mappingKeys(plot.info), function(k)
            paste0(k, ': ', plot.info[k]) )
        
        for ( i in seq_along(plot.info.lines) ) {
            graphics::mtext(plot.info.lines[i], line=plot.info.indices[i], side=3, adj=0,
                cex=1.0, col='black', family='sans', font=1)
        }
    }
    
    return( invisible() )
}

# End of plotQTLScanone.R ######################################################