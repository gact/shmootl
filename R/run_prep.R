# Start of run_prep.R ##########################################################

# run_prep ---------------------------------------------------------------------
#' Prepare data in \pkg{R/qtl} CSV file.
#' 
#' @description This script will take an \pkg{R/qtl} input data file, and
#' prepare it for input to \pkg{shmootl}. Possible prep actions include
#' resolving map sequence labels (\code{'normseq'}), estimating a new
#' genetic map with the \pkg{R/qtl} function \code{'est.map'} (Broman
#' \emph{et al.} 2003), jittering existing map positions with the
#' \pkg{R/qtl} function \code{'jittermap'} (Broman \emph{et al.} 2003),
#' and replacing empty cells with an explicit missing value symbol
#' (\code{'replace.missing'}).
#' 
#' @description If no prep actions are specified, map sequence labels are
#' normalised; a genetic map is estimated if not found, otherwise the existing
#' map is jittered; and empty cells are replaced with an explicit missing
#' value. If one or more prep actions are specified, only the specified
#' actions will be applied. It is the responsibility of the user to ensure
#' that any changes are valid for a given input file.
#' 
#' @description Aside from prep actions, other options include
#' \code{'require.mapunit'}, which indicates if the input file
#' must contain explicit map units; and \code{'include.mapunit'},
#' which indicates if explicit map units should be included in
#' the output file. Options \code{'error.prob'} and \code{'map.function'}
#' are passed to the \pkg{R/qtl} function \code{'est.map'}, while option
#' \code{'amount'} is passed to the \pkg{R/qtl} \code{'jittermap'} function.
#' 
#' @param datafile \pkg{R/qtl} CSV file [required]
#' @param normseq normalise map sequence labels
#' @param estimap estimate genetic map
#' @param error.prob genotyping error rate
#' @param map.function genetic map function
#' @param jittermap jitter map positions
#' @param amount jitter amount
#' @param replace.missing make missing values explicit
#' @param require.mapunit require map units in input
#' @param include.mapunit include map units in output
#' 
#' @template ref-broman-2003
#' 
#' @concept shmootl:preparation
#' @export
#' @family pipeline functions
#' @importFrom utils read.csv
#' @importFrom utils write.table
#' @rdname run_prep
run_prep <- function(datafile=NA_character_, normseq=FALSE, estimap=FALSE,
    error.prob=0.0001, map.function=c('haldane','kosambi','c-f','morgan'),
    jittermap=FALSE, amount=1e-6, replace.missing=FALSE, require.mapunit=FALSE,
    include.mapunit=TRUE) {
    
    stopifnot( isSingleString(datafile) )
    stopifnot( file.exists(datafile) )
    stopifnot( isBOOL(normseq) )
    stopifnot( isBOOL(estimap) )
    stopifnot( isSingleProbability(error.prob) )
    stopifnot( isBOOL(jittermap) )
    stopifnot( isSinglePositiveNumber(amount) )
    stopifnot( isBOOL(replace.missing) )
    stopifnot( isBOOL(require.mapunit) )
    stopifnot( isBOOL(include.mapunit) )
    
    map.function <- match.arg(map.function)
    
    # Set prep actions from corresponding options.
    actions = c(
        normseq         = normseq,
        estimap         = estimap,
        jittermap       = jittermap,
        replace.missing = replace.missing
    )
    
    # If no prep actions were specified, do default prep actions.
    default.prep <- all( actions == FALSE )
    if (default.prep) {
        actions <- replace(actions, seq_along(actions), TRUE)
    }
    
    # Read input CSV file.
    x <- utils::read.csv(datafile, header=FALSE, check.names=FALSE, quote='',
        na.strings=const$missing.value, colClasses='character',
        strip.white=TRUE, stringsAsFactors=FALSE)
    
    # Trim any blank rows/columns from the bottom/right, respectively.
    x <- bstripBlankRows( rstripBlankCols(x) )
    
    params <- getMetadataCSV(x)
    
    if ( ! params$class %in% c('cross', 'geno', 'pheno') ) {
        stop("cannot prep - ", params$class," input data")
    }
    
    # If specified, replace empty cells in data with NA values.
    # NB: this effectively inserts the explicit missing value symbol.
    if (actions['replace.missing']) {
        x.data <- x[params$dat.rows, ]
        x.data[ x.data == '' ] <- NA
        x[params$dat.rows, ] <- x.data
    }
    
    # If this is a cross/geno data file, prep genetic map..
    if ( params$class %in% c('cross', 'geno') ) {
        
        # Get normalised copy of locus sequence labels.
        norm.locus.seqs <- normSeq( as.character(x[2, params$geno.cols]) )
        
        # If specified, ensure sequence labels are normalised.
        if (actions['normseq']) {
            x[2, params$geno.cols] <- norm.locus.seqs
        }
        
        # If map estimation specified, or applied default prep steps to a data
        # file without a genetic map, estimate a jittered genetic map..
        if ( actions['estimap'] || ( default.prep && ! params$map.present ) ) {
            
            head.rows <- if (params$map.present) {1:3} else {1:2}
            
            # Get genotype data-frame, replacing empty cells with NA values.
            geno.data <- x[params$dat.rows, params$geno.cols]
            geno.data[ geno.data == '' ] <- NA
            geno.frame <- rbind(x[head.rows, params$geno.cols], geno.data)
            
            # Insert placeholder 'id' column in genotype data-frame.
            num.samples <- length(params$dat.rows)
            digits <- ceiling( log10(num.samples) )
            suffixes <- formatC(1:num.samples, width=digits, flag=0)
            placeholder.column <- c( 'id', rep_len('', length(head.rows) - 1 ),
                paste0('S', suffixes) )
            geno.frame <- insertColumn(geno.frame, 1, data=placeholder.column)
            
            # Convert genotype data-frame to genotype object.
            geno <- as.geno(geno.frame, require.mapunit=FALSE)
            
            # If genotype object has founder genotypes, estimate genetic map..
            if ( hasFounderGenotypes(geno) ) {
                
                # Create placeholder phenotype object.
                pheno <- makePlaceholderPheno(samples=pull.ind(geno))
                
                # Create temporary cross object.
                cross <- makeCross(geno, pheno)
                
                # Estimate genetic map from temporary cross object.
                gmap <- qtl::est.map(cross, error.prob=error.prob,
                    map.function=map.function, offset=0)
                
                # If specified, ensure no coinciding markers.
                if (actions['jittermap']) {
                    gmap <- qtl::jittermap(gmap)
                }
                
                # Set map positions in initial rows, making sure
                # to include explicit map units, if specified.
                map.blanks <- rep_len('', ncol(x) - length(params$geno.cols))
                output.mapunit <- if (include.mapunit) { 'cM' } else { NULL }
                locus.pos <- setPosColDataMapUnit(pullLocusPos(gmap), output.mapunit)
                map.row <- c(map.blanks, locus.pos)
                x <- rbind(x[1:2, ], map.row, x[params$dat.rows, ])
                
            } else if( ! default.prep ) {
                # ..otherwise check that map estimation was not explicitly specified.
                stop("cannot estimate map from non-founder genotype data")
            }
            
        } else if (params$map.present) {
            
            # Create map table from initial rows.
            map.table <- as.data.frame( t(x[1:3, params$geno.cols]),
                stringsAsFactors=FALSE)
            colnames(map.table) <- const$maptable.colnames[1:3]
            map.table <- setRownamesFromColumn(map.table, col.name='id')
            
            # Get map unit from map positions.
            map.unit <- getPosColDataMapUnit(map.table)
            
            if ( ! is.na(map.unit) ) {
                if ( map.unit != 'cM' ) {
                    stop("map positions must be in centiMorgans (e.g. '47 cM')")
                }
            } else {
                if (require.mapunit) {
                    stop("map positions must include map units (e.g. '47 cM')")
                }
            }
            
            # Convert map table to genetic map object.
            gmap <- as.map(map.table, map.unit='cM')
            
            # If specified, ensure no coinciding markers.
            if (actions['jittermap']) {
                gmap <- qtl::jittermap(gmap)
            }
            
            # Replace map positions in initial rows, making
            # sure to include explicit map units, if specified.
            output.mapunit <- if (include.mapunit) { 'cM' } else { NULL }
            locus.pos <- setPosColDataMapUnit(pullLocusPos(gmap), output.mapunit)
            x[3, params$geno.cols] <- locus.pos
        }
    }
    
    # Create temp file.
    tmp <- tempfile()
    on.exit( file.remove(tmp) )
    
    # Write prepped data to temp file.
    utils::write.table(x, file=tmp, na=const$missing.value, sep=',',
        quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    # Move temp file to final output file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, datafile, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_prep.R ############################################################