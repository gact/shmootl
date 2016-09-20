# Start of run_prep.R ##########################################################

# run_prep ---------------------------------------------------------------------
#' Prepare data in \pkg{R/qtl} CSV file.
#' 
#' This script will take an \pkg{R/qtl} input data file, and prepare it
#' for input to \pkg{shmootl}. Possible prep actions include resolving
#' map sequence labels (\code{'normseq'}), jittering map positions (with
#' the \pkg{R/qtl} function \code{jittermap}), and replacing empty cells
#' with an explicit missing value symbol (\code{'replace.missing'}). If
#' no prep actions are specified, all are applied. Other options include
#' \code{'require.mapunit'}, which indicates if the input file must contain
#' explicit map units; and \code{'include.mapunit'}, which indicates if
#' explicit map units should be included in the output file. It is the
#' responsibility of the user to ensure that these changes are valid for
#' a given input file.
#' 
#' @param datafile \pkg{R/qtl} CSV file
#' @param jittermap jitter map positions
#' @param normseq normalise map sequence labels
#' @param replace.missing make missing values explicit
#' @param require.mapunit require map units in input
#' @param include.mapunit include map units in output
#' 
#' @concept shmootl:utilities
#' @export
#' @family pipeline functions
#' @importFrom utils read.csv
#' @importFrom utils write.table
#' @rdname run_prep
run_prep <- function(datafile, jittermap=FALSE, normseq=FALSE,
    replace.missing=FALSE, require.mapunit=FALSE, include.mapunit=TRUE) {
    
    stopifnot( isSingleString(datafile) )
    stopifnot( file.exists(datafile) )
    stopifnot( isBOOL(jittermap) )
    stopifnot( isBOOL(normseq) )
    stopifnot( isBOOL(replace.missing) )
    stopifnot( isBOOL(require.mapunit) )
    stopifnot( isBOOL(include.mapunit) )
    
    # Set prep actions from corresponding options.
    actions = c(
        jittermap       = jittermap,
        normseq         = normseq,
        replace.missing = replace.missing
    )
    
    # If none of the prep actions were specified, perform all prep actions.
    if ( all( actions == FALSE ) ) {
        actions <- replace(actions, seq_along(actions), TRUE)
    }
    
    # Read input CSV file.
    x <- utils::read.csv(datafile, header=FALSE, check.names=FALSE, quote='',
        stringsAsFactors=FALSE, strip.white=TRUE, na.strings=const$missing.value)
    
    # Trim any blank rows/columns from the bottom/right, respectively.
    x <- bstripBlankRows( rstripBlankCols(x) )
    
    params <- getMetadataCSV(x)
    
    if ( ! params$class %in% c('cross', 'geno', 'pheno') ) {
        stop("cannot prep - ", params$class," input data")
    }
    
    # If this is a cross/geno data file, ensure genetic map units present..
    if ( params$class %in% c('cross', 'geno') ) {
        
        # TODO: estimate map if not present ('estimap' option)
        
        # If specified, ensure sequence labels are normalised.
        if (actions['normseq']) {
            x[2, params$geno.cols] <- normSeq( as.character(x[2, params$geno.cols]) )
        }
        
        if (params$map.present) {
            
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
            
            # Replace map markers/positions in initial rows,
            # including explicit map units if specified.
            output.mapunit <- if (include.mapunit) { 'cM' } else { NULL }
            map.table <- setPosColDataMapUnit(as.data.frame(gmap), output.mapunit)
            map.table <- setColumnFromRownames(map.table, col.name='id')
            x[c(1,3), params$geno.cols] <- t(map.table)[c(1,3), ]
        }
    }
    
    # If specified, replace empty cells in data with NA values.
    # NB: this effectively inserts the explicit missing value symbol.
    if (actions['replace.missing']) {
        x.data <- x[params$dat.rows, ]
        x.data[ x.data == '' ] <- NA
        x[params$dat.rows, ] <- x.data
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