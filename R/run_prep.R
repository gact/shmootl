# Start of run_prep.R ##########################################################

# run_prep ---------------------------------------------------------------------
#' Prepare data in \pkg{R/qtl} CSV file.
#' 
#' @description This script will take an \pkg{R/qtl} input data file, and
#' prepare it for input to \pkg{shmootl}. Marker sequence IDs are resolved,
#' and empty cells are replaced with a missing value symbol. If a genetic
#' map is present, map positions are jittered and map units are added. It
#' is the responsibility of the user to ensure that these changes are valid
#' for a given input file.
#' 
#' @param datafile \pkg{R/qtl} CSV file
#' 
#' @export
#' @importFrom utils read.csv
#' @importFrom utils write.table
#' @rdname run_prep
run_prep <- function(datafile) {
    
    stopifnot( isSingleString(datafile) )
    stopifnot( file.exists(datafile) )
    
    # Read input CSV file.
    x <- utils::read.csv(datafile, header=FALSE, check.names=FALSE, quote='',
        stringsAsFactors=FALSE, strip.white=TRUE, na.strings=const$missing.value)
    
    # Trim any blank rows/columns from the bottom/right, respectively.
    x <- bstripBlankRows( rstripBlankCols(x) )
    
    params <- getMetadataCSV(x)
    
    if ( is.null(params$class) ) {
        stop("cannot prep - unknown input data")
    }
    
    if ( ! params$class %in% c('cross', 'geno', 'pheno') ) {
        stop("cannot prep ", params$class," data")
    }
    
    # If file has genotype columns, ensure marker sequence IDs are normalised.
    if ( ! is.null(params$geno.cols) ) {
        x[2, params$geno.cols] <- normSeq( as.character(x[2, params$geno.cols]) )
    }
    
    # If this is a cross/geno data file, ensure genetic map units present..
    if ( params$class %in% c('cross', 'geno') ) {
        
        # TODO: estimate map if not present
        
        if (params$map.present) {
            
            head.rows <- 1:3
            
            # Create map table from initial rows.
            map.table <- as.data.frame( t(x[head.rows, params$geno.cols]),
                stringsAsFactors=FALSE)
            colnames(map.table) <- const$maptable.colnames[head.rows]
            map.table <- setRownamesFromColumn(map.table, col.name='id')
            
            # Get map unit from map positions.
            map.unit <- getPosColDataMapUnit(map.table)
            
            # Check map unit if present.
            if ( ! is.na(map.unit) && map.unit != 'cM' ) {
                stop("map positions must be in centiMorgans (e.g. '47 cM')")
            }
            
            # Convert map table to genetic map object.
            gmap <- as.map(map.table, map.unit='cM')
            
            # Ensure no coinciding markers.
            gmap <- qtl::jittermap(gmap)
            
            # Replace map in initial rows.
            map.table <- setPosColDataMapUnit(as.data.frame(gmap), 'cM')
            map.table <- setColumnFromRownames(map.table, col.name='id')
            x[head.rows, params$geno.cols] <- t(map.table)
        }
    }
    
    # Replace empty cells in data with NA values.
    # NB: this will effectively insert the missing value symbol.
    x.data <- x[params$dat.rows, ]
    x.data[ x.data == '' ] <- NA
    x[params$dat.rows, ] <- x.data
    
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