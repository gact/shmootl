# Start of plotQTLScantwo.R ####################################################

# plotQTLScantwo ---------------------------------------------------------------
#' Plot \pkg{R/qtl} \code{scantwo} result.
#' 
#' This function is a restricted version of \pkg{R/qtl} function
#' \code{plot.scantwo} (Broman \emph{et al.} 2003), plotting the
#' input \code{scantwo} result as a heatmap. The set of palettes
#' are the same as in \pkg{R/qtl} \code{plot.scantwo}, except for
#' \code{'viridis'} which uses the recently-developed \code{R}
#' package of the same name.
#' 
#' @param x An \pkg{R/qtl} \code{scantwo} object.
#' @param chr Vector indicating which sequences to plot.
#' If none are specified, all are plotted.
#' @param incl.markers When \code{scantwo} has been run with a fixed non-zero
#' step width, this indicates if markers (as distinct from pseudomarkers)
#' should be included in the plot.
#' @template param-lodcolumn
#' @param lower The \code{scantwo} LOD type to
#' output in the lower triangle of the plot.
#' @param upper The \code{scantwo} LOD type to
#' output in the upper triangle of the plot.
#' @param zscale Display a z-scale in the right margin of the plot.
#' @param palette Palette to use in the plot. This is equivalent to the
#' \pkg{R/qtl} \code{plot.scantwo} parameter \code{'col.scheme'}.
#' @param phenotype Name of the phenotype, to be shown in plot information.
#' @param ... Unused arguments.
#' 
#' @template ref-broman-2003
#' 
#' @seealso Function \code{plot.scantwo} in the
#' \href{http://www.rqtl.org/manual/qtl-manual.pdf}{R/qtl manual}.
#' 
#' @export
#' @family plot functions
#' @importFrom graphics image
#' @rdname plotQTLScantwo
plotQTLScantwo <- function(x, chr=NULL, incl.markers=FALSE, lodcolumn=NULL,
    lower=c('full', 'add', 'cond-int', 'cond-add', 'int'),
    upper=c('int', 'cond-add', 'cond-int', 'add', 'full'), zscale=TRUE,
    palette=c('viridis', 'redblue', 'cm', 'gray', 'heat', 'terrain', 'topo'),
    phenotype=NULL, ...) {
    
    opar <- graphics::par(no.readonly=TRUE)
    on.exit( graphics::par(opar), add=TRUE )
    
    stopifnot( 'scantwo' %in% class(x) )
    stopifnot( emptyArgs(...) )
    
    palette <- match.arg(palette)
    
    if ( isSingleString(lower) ) {
        lower <- resolveScantwoLodtypes(lower, to='plot.scantwo')
    }
    lower <- match.arg(lower)
    
    if ( isSingleString(upper) ) {
        upper <- resolveScantwoLodtypes(upper, to='plot.scantwo')
    }
    upper <- match.arg(upper)
    
    # Generate requested palette.
    palette.size <- 256
    gamma <- 0.6
    if ( palette == 'viridis' ) {
        
        if ( ! 'viridis' %in% rownames(utils::installed.packages()) ) {
            stop("cannot get 'viridis' palette without R package 'viridis'")
        }
        
        col <- viridis::viridis(palette.size)
        
    } else  if ( palette == 'redblue' ) {
        
        hex.col <- rev( grDevices::rainbow(palette.size, start=0, end=0.666) )
        rgb.col <- ( grDevices::col2rgb(hex.col) / (palette.size - 1) ) ^ gamma
        col <- grDevices::rgb(rgb.col[1, ], rgb.col[2, ], rgb.col[3, ])
        
    } else  if ( palette == 'gray' ) {
        
        if ( isPositiveNumber(gamma) ) {
            col <- rev( grDevices::gray( log( seq( from=1, to=exp(gamma),
                len=palette.size ) ) / gamma ) )
        } else {
            col <- rev( grDevices::gray( seq(0.0, 1.0, len=palette.size) ) )
        }
        
    } else  if ( palette == 'cm' ) {
        col <- grDevices::cm.colors(palette.size)
    } else  if ( palette == 'heat' ) {
        col <- grDevices::heat.colors(palette.size)
    } else  if ( palette == 'terrain' ) {
        col <- grDevices::terrain.colors(palette.size)
    } else  if ( palette == 'topo' ) {
        col <- grDevices::topo.colors(palette.size)
    }
    
    # Get scantwo result LOD column indices.
    lodcol.index <- getLodColIndex(x, lodcolumn=lodcolumn)
    
    # Ensure scantwo result has normalised sequence IDs.
    x <- normSeq(x)
    
    # Get specified sequences.
    x.seqs <- pull.chr(x$map)
    chr <- x.seqs <- subsetBySeq(x.seqs, chr)
    stopifnot( length(chr) > 0 )
    
    # Subset scantwo map by specified sequences.
    x.map <- subsetBySeq(x$map, chr)
    
    # Remove markers from map, if requested.
    if ( ! incl.markers && any( x.map[, 'eq.spacing'] == 0 ) ) {
        x.map <- x.map[ x.map[, 'eq.spacing'] == 1, ]
    }
    
    # Get the LOD matrix of the scantwo object, subsetting by the
    # specified sequences, and with markers removed if requested.
    lod <- getLODMatrix(x, chr=chr, incl.markers=incl.markers,
        lodcolumn=lodcolumn, lower=lower, upper=upper)
    
    L <- lower.tri(lod)
    U <- upper.tri(lod)
    
    diag(lod) <- 0 # NB: single-QTL LOD values are not currently shown.
    
    # Get LOD value ranges of matrix.
    min.lod <- 0
    upper.max.lod <- max( lod[U][ is.finite(lod[U]) ] )
    lower.max.lod <- max( lod[L][ is.finite(lod[L]) ] )
    max.lod <- max( lod[ is.finite(lod) ] )
    
    # Get z-axis ranges of heatmap.
    zlim.lower <- c(min.lod, lower.max.lod)
    zlim.upper <- c(min.lod, upper.max.lod)
    zlim <- c(min.lod, max.lod)
    
    missing.lods <- is.na(lod)
    if( any(missing.lods) ) {
        warning("setting ", sum(missing.lods), " missing LOD scores to zero")
        lod[missing.lods] <- 0
    }
    
    negative.lods <- isNegativeNumber(lod)
    if( any(negative.lods) ) {
        warning("setting ", sum(negative.lods), " negative LOD scores to zero")
        lod[negative.lods] <- 0
    }
    
    infinite.lods <- is.infinite(lod)
    if( any(infinite.lods) ) {
        warning("setting ", sum(infinite.lods), " infinite LOD scores to",
            " maximum finite LOD value (", max.lod, ")")
        lod[infinite.lods] <- max.lod
    }
    
    # Rescale upper- and lower-triangle LOD scores so
    # that they both use the full range of the palette.
    lod[L] <- lod[L] * ( max.lod / lower.max.lod )
    lod[U] <- lod[U] * ( max.lod / upper.max.lod )
    
    # Assemble sequence plotting info ------------------------------------------
    
    seq.par <- matrix( NA_real_, nrow=length(chr), ncol=3,
        dimnames=list(chr, c('offset', 'midpoint', 'length') ) )
    
    if ( length(chr) > 1 ) {
        
        seq.indices <- matchSeqRowIndices(x.map, chr)
        seq.nmar <- unlist( lengths(seq.indices) )
        seq.bounds <- c(0.5, cumsum(seq.nmar) + 0.5)
        seq.par[, 'offset'] <- seq.bounds[ -length(seq.bounds) ]
        seq.par[, 'length'] <- diff(seq.bounds)
        seq.par[, 'midpoint'] <- seq.par[, 'offset'] + (0.5 * seq.par[, 'length'])
        
    } else {
        
        seq.par[, 'offset'] <- 0
        seq.par[, 'length'] <- x.map[nrow(x.map), 'pos']
        seq.par[, 'midpoint'] <- seq.par[, 'offset'] + (0.5 * seq.par[, 'length'])
    }
    
    # Assemble general plot info -----------------------------------------------
    
    plot.info <- mapping()
    
    # Get phenotype from argument if specified, or LOD column otherwise.
    if ( ! is.null(phenotype) ) {
        stopifnot( isSingleString(phenotype) )
        plot.info['Phenotype'] <- phenotype
    } else {
        plot.info['Phenotype'] <- attr(x, 'phenotypes')[lodcol.index]
    }
    
    # For single-sequence plots, add sequence label to plot info.
    if ( length(chr) == 1 ) {
        plot.info['Chromosome'] <- chr
    }
    
    # Add LOD types for upper and lower triangles.
    plot.info['Upper Triangle'] <- resolveScantwoLodtypes(upper, to='title')
    plot.info['Lower Triangle'] <- resolveScantwoLodtypes(lower, to='title')
    
    # --------------------------------------------------------------------------
    
    # Set margin line counts, adjusting for plot info lines and z-scale.
    top.mar <- length(plot.info) + 4 # NB: minimum 4 margin lines
    right.mar <- if (zscale) { 5 } else { 2 }
    mar <- (c(5, 4, top.mar, right.mar) + 0.1)
    
    # Set plot info margin indices.
    plot.info.indices <- rev( seq_along(plot.info) )
    plot.title.index <- top.mar - 2
    
    # Set plot dimensions.
    xlim <- ylim <- c(seq.par[1, 'offset'], seq.par[nrow(seq.par), 'offset'] +
        seq.par[nrow(seq.par), 'length'])
    
    # Set plotting parameters.
    if ( length(chr) > 1 ) {
        gridlines <- getRowIndices(lod)
        xlab <- ylab <- 'Chromosome'
        xaxt <- yaxt <- 'n'
        
    } else {
        gridlines <- pullLocusPos(x.map)
        xlab <- ylab <- 'Location (cM)'
        xaxt <- yaxt <- 's'
    }
    
    # Set fixed plotting arguments.
    fixed.args <- list(x=gridlines, y=gridlines, z=lod, col=col, family='sans',
        fg='black', xaxs='i', xaxt=xaxt, xlim=xlim, yaxs='i', yaxt=yaxt,
        ylim=ylim, zlim=zlim)
    
    # Set default plotting arguments.
    default.args <- list(xlab=xlab, ylab=ylab)
    
    # Set args from fixed and default values.
    # TODO: allow user to set some arguments
    args <- c(fixed.args, default.args)
    
    # Set plot graphics parameters.
    graphics::par(
        las = 1,   # horizontal axis labels
        mar = mar, # margin line counts
        pty = 's'  # square plot
    )
    
    # Plot heatmap of two-QTL LOD matrix.
    # NB: in R, a matrix is usually printed from top-left to bottom-right, but
    # the 'image' function plots a matrix heatmap from bottom-left to top-right.
    do.call('image', args)
    
    # If multiple sequences, add sequence ticks, labels, and gridlines.
    if( length(chr) > 1 ) {
        
        # TODO: adjust size of sequence labels so that they don't overlap.
        
        # Plot x- and y-axes with given ticks and labels.
        for ( side in 1:2 ) {
            for ( i in seq_along(chr) ) {
                graphics::axis(side=side, at=seq.par[i, 'midpoint'],
                    labels=chr[i], cex.axis=1.0, col.axis='black', font.axis=1)
            }
        }
        
        # Plot gridlines between sequences.
        graphics::abline(h=seq.par[, 'offset'], xpd=FALSE)
        graphics::abline(v=seq.par[, 'offset'], xpd=FALSE)
    }
    
    # Plot diagonal line separating upper and lower triangles.
    graphics::segments(xlim[1], xlim[3], xlim[2], xlim[4], lwd=3)
    
    # Plot box around graph.
    graphics::box(lwd=3)
    
    # Write plot title.
    graphics::title(main='Scantwo', line=plot.title.index, cex.main=1.2,
        col.main='black', font.main=2)
    
    # If any plot info, add to top margin.
    if ( length(plot.info) > 0 ) {
        
        plot.info.lines <- sapply( mappingKeys(plot.info), function(k)
            paste0(k, ': ', plot.info[k]) )
        
        for ( i in seq_along(plot.info.lines) ) {
            graphics::mtext(plot.info.lines[i], line=plot.info.indices[i], side=3, adj=0,
                cex=1.0, col='black', font=1)
        }
    }
    
    # Plot z-scale in right margin, if requested.
    if (zscale) {
        
        # Set z-scale graphical parameters.
        zpar <- list(cex=0.8, col='black', font=1, lty=1, lwd=1)
        
        # Get one hundredth of plot size in user-space units,
        # which is used to define z-scale position and size.
        cp <- diff(xlim) / 100
        
        # Get x- and y-coordinates of z-scale rectangle.
        x1 <- xlim[2] + 10*cp
        x2 <- x1 + 5*cp
        y1 <- ylim[1]
        y2 <- ylim[2]
        
        # Plot z-scale gradient.
        breaks <- seq(y1, y2, length.out=256)
        for ( i in seq( length(breaks) - 1 ) ) {
            graphics::rect(x1, breaks[i], x2, breaks[i+1],
                border=NA, col=col[i], xpd=TRUE)
        }
        
        # Plot z-scale border.
        graphics::rect(x1, y1, x2, y2, border=zpar$col, col=NA,
            lty=zpar$lty, lwd=zpar$lwd, xpd=TRUE)
        
        # Plot lower- and upper-triangle labels.
        midpoint <- y1 + 0.5 * (y2 - y1)
        graphics::text(x2 + 2.5*cp, y=midpoint, labels=plot.info['Lower Triangle'],
            adj=0.5, cex=zpar$cex, col=zpar$col, font=zpar$font, xpd=TRUE, srt=270)
        graphics::text(x1 - 2.5*cp, y=midpoint, labels=plot.info['Upper Triangle'],
            adj=0.5, cex=zpar$cex, col=zpar$col, font=zpar$font, xpd=TRUE, srt=90)
        
        # Plot tick marks for LOD score minimum and maximum values.
        t1 <- x1 - cp
        t2 <- x2 + cp
        graphics::segments(
            x0 = c(t1, t1), y0 = c(y1, y2),
            x1 = c(t2, t2), y1 = c(y1, y2),
        col=zpar$col, lty=zpar$lty, lwd=zpar$lwd, xpd=TRUE)
        
        # Plot lower-triangle LOD score minimum and maximum values.
        graphics::text(x2, y1, labels=round(zlim.lower[1], 1), pos=4,
            cex=zpar$cex, col=zpar$col, font=zpar$font, xpd=TRUE)
        graphics::text(x2, y2, labels=round(zlim.lower[2], 1), pos=4,
            cex=zpar$cex, col=zpar$col, font=zpar$font, xpd=TRUE)
        
        # Plot upper-triangle LOD score minimum and maximum values.
        graphics::text(x1, y1, labels=round(zlim.upper[1], 1), pos=2,
            cex=zpar$cex, col=zpar$col, font=zpar$font, xpd=TRUE)
        graphics::text(x1, y2, labels=round(zlim.upper[2], 1), pos=2,
            cex=zpar$cex, col=zpar$col, font=zpar$font, xpd=TRUE)
    }
    
    return( invisible() )
}

# End of plotQTLScantwo.R ######################################################