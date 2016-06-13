# Start of csv.R ###############################################################

# CSV Utilities ----------------------------------------------------------------
#' CSV input/output utilities.
#' 
#' Functions for validated input and output of \pkg{R/qtl} data in CSV format.
#' These include functions for reading (\code{\link{readCrossCSV}}) or writing 
#' (\code{\link{writeCrossCSV}}) an \pkg{R/qtl} \code{cross} object. There are 
#' also functions for reading (\code{\link{readMapCSV}}) or writing 
#' (\code{\link{writeMapCSV}}) an \pkg{R/qtl} \code{map} object.
#' 
#' @docType package
#' @name CSV Utilities
NULL

# getMetadataCSV -------------------------------------------------------------
#' Get parameters of \pkg{R/qtl} input data.
#' 
#' @description Given a data frame as loaded from an \pkg{R/qtl} input CSV file,
#' this function determines the parameters of that data (e.g. phenotype columns).
#' 
#' @details The returned list includes information on the following parameters:
#' 
#' \itemize{
#'  \item{\code{pheno.cols}}{Indices of phenotype columns, or NA if none found.}
#'  \item{\code{id.col}}{Index of ID column, or NA if not found.}
#'  \item{\code{geno.cols}}{Indices of genotype columns, or NA if none found.}
#'  \item{\code{dat.rows}}{Row indices of data.}
#'  \item{\code{map.present}}{TRUE if the data contains a map; FALSE otherwise.}
#'  \item{\code{class}}{Expected class of the input data.}
#' }
#' 
#' This information can then guide further processing of the input data.
#' 
#' @param x An \pkg{R/qtl} input data frame as read from a CSV file.
#'
#' @return A list containing parameters of \pkg{R/qtl} input data.
#' 
#' @keywords internal
#' @rdname getMetadataCSV
getMetadataCSV <- function(x) {
    
    stopifnot( is.data.frame(x) )
    
    if ( nrow(x) < 2 ) {
        stop("data not found")
    }
    
    if ( anyNA(x[1:2, ]) ) {
        stop("invalid headings found")
    }
    
    # Get logical vector indicating which columns are blank
    # in the first row after the initial heading row.
    seq.blanks <- ! is.na(x[2, ]) & x[2, ] == ''
    
    # If third row present, get indices of columns with blanks
    # (where the optional genetic map data is placed)..
    if ( nrow(x) >= 3 ) {
        map.blanks <- ! is.na(x[3, ]) & x[3, ] == ''
    } else { # ..otherwise assume no map data.
        map.blanks <- rep.int(FALSE, ncol(x))
    }
    
    # Verify genetic map correctly formatted, if present.
    if ( any(map.blanks) && ( length(map.blanks) != length(seq.blanks) ||
        any(map.blanks != seq.blanks) ) ) {
        stop("map data row is invalid")
    }
    
    # Set phenotype columns from those with blank sequence row.
    # This includes the ID column for now.
    pheno.cols <- which(seq.blanks)
    
    # Get index of 'id' columns.
    id.col <- which(tolower(x[1, ]) == 'id')
    
    if ( length(id.col) > 1 ) {
        stop("multiple ID columns found")
    }
    
    # If there are phenotype columns, this is either cross or geno data..
    if ( length(pheno.cols) > 0 ) {
        
        if ( pheno.cols[1] != 1 || any(diff(pheno.cols) != 1) ) {
            stop("sequence data row is invalid")
        } 
        
        if ( length(id.col) > 0 ) {
            pheno.cols <- pheno.cols[-id.col]
        }
        
        geno.cols <- which( ! seq.blanks )
        
        if ( length(geno.cols) == 0 ) {
            stop("genotype data not found")
        }
        
        map.present <- any(map.blanks)
        head.rows <- if (map.present) {1:3} else {1:2}
        
        if ( length(pheno.cols) == 0 ) {
            
            if ( length(id.col) == 0 ) {
                stop("ID column not found")
            }
            
            data.class <- 'geno'
            
        } else {
            
            data.class <- 'cross'
        }
        
    } else { # ..otherwise it's either map or pheno data.
        
        if ( length(id.col) == 0 ) {
            stop("ID column not found")
        }
        
        geno.cols <- NA
        head.rows <- 1
        
        seq.col <- which( x[1, ] == 'chr' )
        pos.col <- which(grepl(const$pattern$poscol, x[1, ], ignore.case=TRUE))
        
        if ( length(seq.col) == 1 && length(pos.col) == 1 &&
             id.col == 1 && seq.col == 2 && pos.col == 3 ) {
            
            data.class <- 'map'
            map.present <- TRUE
            
        } else {
            
            data.class <- 'pheno'
            map.present <- FALSE
            pheno.cols <- getColIndices(x) # every column is a pheno column
        }
    }
    
    first.data.row <- length(head.rows) + 1
    last.data.row <- nrow(x)
    
    if ( last.data.row < first.data.row ) {
        stop("data not found")
    }
    
    dat.rows <- first.data.row : last.data.row
    
    if ( length(pheno.cols) == 0 ) {
        pheno.cols <- NA
    }
    
    if ( length(id.col) == 0 ) {
        id.col <- NA
    }
    
    params <- list(
        pheno.cols  = pheno.cols,
        id.col      = id.col,
        geno.cols   = geno.cols,
        dat.rows    = dat.rows,
        map.present = map.present,
        class       = data.class
    )
    
    return(params)
}

# hasMapCSV --------------------------------------------------------------------
#' Test if CSV file contains map data.
#' 
#' @param infile Input CSV file path.
#'
#' @return TRUE if CSV file contains map data in a known format;
#' FALSE otherwise.
#'
#' @export
#' @importFrom utils read.csv
#' @keywords internal
#' @rdname hasMapCSV
hasMapCSV <- function(infile) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    
    x <- utils::read.csv(infile, header=FALSE, nrows=4, check.names=FALSE, quote='',
        stringsAsFactors=FALSE, strip.white=TRUE, na.strings=const$missing.value)
    
    x <- rstripBlankCols(x)
    
    status <- FALSE
    
    tryCatch({
        params <- getMetadataCSV(x)
        status <- params$map.present
    }, error=function(e) {})
    
    return(status)
}

# readCrossCSV -----------------------------------------------------------------
#' Read yeast \code{cross} from a CSV file.
#' 
#' This function reads yeast cross data from an \pkg{R/qtl} CSV file and returns 
#' an \pkg{R/qtl} \code{cross} object with an attribute \code{'info'} of type 
#' \code{CrossInfo}. 
#' 
#' @param infile Input CSV file path.
#' @param error.prob Genotyping error rate (ignored unless estimating genetic map).
#' @param map.function Genetic map function (ignored unless estimating genetic map).
#' @param require.mapunit Require map unit information in map positions.
#'  
#' @return An \pkg{R/qtl} \code{cross} object with an attribute \code{'info'} of 
#' type \code{CrossInfo}.
#' 
#' @export
#' @family csv utilities
#' @importFrom methods new
#' @importFrom utils read.csv
#' @importFrom utils write.table
#' @rdname readCrossCSV
readCrossCSV <- function(infile, error.prob=0.0001,
    map.function=c('haldane', 'kosambi', 'c-f', 'morgan'),
    require.mapunit=TRUE) {
    
    crosstype <- 'bc' # TODO: support other cross types
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    error.prob <- as.numeric(error.prob)
    stopifnot( isSingleProbability(error.prob) )
    stopifnot( isBOOL(require.mapunit) )
    
    map.function <- match.arg(map.function)
    
    # Read cross input data as CSV file.
    cross.table <- utils::read.csv(infile, header=FALSE, check.names=FALSE, quote='',
        stringsAsFactors=FALSE, strip.white=TRUE, na.strings=const$missing.value)
    
    # Trim any blank rows/columns from the bottom/right, respectively.
    cross.table <- bstripBlankRows( rstripBlankCols(cross.table) )
    
    # Get logical vector indicating which columns are blank 
    # in the first row after the initial heading row.
    seq.is.blank <- cross.table[2, ] == ''
    
    # Set phenotype columns from those with blank sequence row.   
    pheno.cols <- which( seq.is.blank )
    stopifnot( length(pheno.cols) > 0 )
    
    # Set marker columns from those with nonempty sequence row.   
    geno.cols <- which( ! seq.is.blank )
    stopifnot( length(geno.cols) > 0 )
    
    # Verify that phenotype columns form a contiguous block at left of table.
    if ( pheno.cols[1] != 1 || any(diff(pheno.cols) != 1) ) {
        stop("sequence data row is invalid")
    } 
    
    # Get indices of columns with blanks in second row after initial  
    # headings (where the optional genetic map data is placed).
    map.blanks <- which( cross.table[3, ] == '' )
    
    # Check whether genetic map data appears to be present.
    map.present <- length(map.blanks) > 0
    
    # Verify genetic map correctly formatted, if present.
    if ( map.present && ! all(map.blanks == pheno.cols) ) {
        stop("map data row is invalid")
    }
    
    # Get indices of initial rows.
    head.rows <- if (map.present) {1:3} else {1:2}
    
    # Get offset of main body of cross data.
    dat.offset = length(head.rows)
    
    # Get index of first and last data rows.
    first.data.row <- dat.offset + 1
    last.data.row <- nrow(cross.table)
    
    stopifnot( last.data.row >= first.data.row )
    
    # Get vector of data row indices.
    dat.rows <- first.data.row : last.data.row
    
    # Verify that there are no blank phenotype values.
    if ( any(cross.table[dat.rows, pheno.cols] == '') ) {
        stop("blank phenotype values found")
    }
    
    # Get set of unique symbols used in marker genotype data.
    geno.symbols <- sort( unique( as.character( 
        unlist(cross.table[dat.rows, geno.cols]) ) ) )
    
    # Verify that there are no blank genotype entries.
    if ( any(geno.symbols == '') ) {
        stop("blank genotype values found")
    }
    
    # Get putative genotypes: genotype symbols that are not 'missing data' symbols.
    genotypes <- geno.symbols[ ! geno.symbols %in% const$missing.value ]
    
    # Verify that genotypes are haploid.
    # TODO: handle other ploidies.
    if ( any( nchar(genotypes) > 1 ) ) {
        stop("unsupported genotype ploidy")
    }
        
    # Verify that there are exactly two genotypes.
    # TODO: handle cross with more than two genotypes.
    if ( length(genotypes) != 2 ) { 
        stop("unsupported number of genotypes - '", length(genotypes), "'")
    }
    
    # Get allele symbols from characters in genotype symbols.
    alleles <- unique( unlist( strsplit(genotypes, '') ) )
    
    # Get phenotype names.
    phenotypes <- as.character(cross.table[1, pheno.cols])
    
    # Get index of phenotype columns with R/qtl special heading 'id'. If present,
    # this contains identifiers of sampled individuals. Search is case-insensitive.
    id.col <- which( tolower(phenotypes) == 'id' )
    
    if ( length(id.col) > 0 ) {
        # If 'id' column present, set vector of sample IDs 
        # and remove column from phenotype column indices..
        stopifnot( length(id.col) == 1 )
        samples <- cross.table[dat.rows, id.col]
        pheno.cols <- pheno.cols[-id.col]
        phenotypes <- as.character(cross.table[1, pheno.cols])
    
    } else {
        # ..otherwise set vector of sample indices.
        samples <- getIndices(dat.rows)
    }
    
    # Get locus info.
    locus.ids <- as.character(cross.table[1, geno.cols])
    locus.seqs <- as.character(cross.table[2, geno.cols])
    
    # Create CrossInfo object.
    cross.info <- methods::new('CrossInfo')
    cross.info <- setMarkers(cross.info, markers=locus.ids)
    cross.info <- setMarkerSeqs(cross.info, sequences=locus.seqs)
    cross.info <- setPhenotypes(cross.info, phenotypes)
    cross.info <- setAlleles(cross.info, alleles)
    cross.info <- setGenotypes(cross.info, genotypes)
    cross.info <- setCrossType(cross.info, crosstype)
    cross.info <- setSamples(cross.info, samples)
    
    # Get normalised marker sequences, replace these in original table.
    locus.seqs <- getMarkerSeqs(cross.info)
    cross.table[2, geno.cols] <- locus.seqs
    
    # Set sorted sequences in CrossInfo object.
    sequences <- sortSeq( unique(locus.seqs) )
    cross.info <- setSequences(cross.info, sequences)
    
    # If map is present, check for map unit info.
    if (map.present) {
        
        # Get cross map positions.
        map.pos <- as.character(cross.table[3, geno.cols])
        
        # Get map unit from map positions.
        map.unit <- getMapUnitSuffix(map.pos)
        
        if ( ! is.na(map.unit) ) {
            if ( map.unit != 'cM' ) {
                stop("cross map positions must be in centiMorgans (e.g. '47 cM')")
            }
        } else {
            if (require.mapunit) {
                stop("cross map positions must include map units (e.g. '47 cM')")
            }
        }
        
        # Strip map unit information.
        map.pos <- setPosColDataMapUnit(map.pos, NULL)
        
        # Replace original map positions.
        cross.table[3, geno.cols] <- map.pos
    }
    
    # Create temp file for adjusted cross data.
    temp.file <- tempfile(fileext='csv')
    utils::write.table(cross.table, file=temp.file, na='', sep=',',
        quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    # Read adjusted cross data.
    cross <- qtl::read.cross('csv', '', temp.file, genotypes=genotypes, 
        alleles=alleles, na.strings=const$missing.value, error.prob=error.prob,
        map.function=map.function)
    
    # Infer strain indices..keep if there are replicate samples.
    strain.indices <- inferStrainIndices(cross)
    if ( max( table(strain.indices) ) > 1 ) {
        cross.info <- setStrainIndices(cross.info, strains=strain.indices)
    }
    
    # Infer tetrad indices..keep if samples are tetradic.
    tetrad.indices <- inferTetradIndices(cross)
    if ( ! is.null(tetrad.indices) ) {
        cross.info <- setTetradIndices(cross.info, tetrads=tetrad.indices)
    }
    
    attr(cross, 'info') <- cross.info
    
    return(cross)
}

# readGenoCSV ------------------------------------------------------------------
#' Read yeast genotype data from a CSV file.
#' 
#' @param infile Input CSV file path.
#' @param require.mapunit Require map unit information in map positions.
#' 
#' @export
#' @family csv utilities
#' @importFrom utils read.csv
#' @rdname readGenoCSV
readGenoCSV <- function(infile, require.mapunit=TRUE) {
    
    stopifnot( isSingleString(infile) )
    
    # Read genotype input data as CSV file.
    geno.table <- utils::read.csv(infile, header=FALSE, check.names=FALSE, quote='',
        stringsAsFactors=FALSE, strip.white=TRUE, colClasses='character',
        na.strings=const$missing.value)
    
    # Trim any blank rows/columns from the bottom/right, respectively.
    geno.table <- bstripBlankRows( rstripBlankCols(geno.table) )
    
    if ( nrow(geno.table) < 3 || ncol(geno.table) < 2 ) {
        stop("invalid genotype data file")
    }
    
    # Make geno table column names from first row of input table.
    colnames(geno.table) <- make.names(geno.table[1, ])
    
    # Check for ID heading in first column.
    id.col <- which( tolower( colnames(geno.table) ) == 'id' )
    if ( length(id.col) == 0 || id.col[1] != 1 ) {
        stop("ID column not found in genotype data file - '", infile, "'")
    }
    
    # Check for sample IDs - required in genotype file.
    head.rows <- if (geno.table[3, id.col] == '') {1:3} else {1:2}
    dat.rows <- getRowIndices(geno.table)[-head.rows]
    if ( any(geno.table[dat.rows, id.col] == '') ) {
        stop("ID column is incomplete in genotype data file - '", infile, "'")
    }
    
    return( as.geno(geno.table, require.mapunit=require.mapunit) )
}

# readMapCSV -------------------------------------------------------------------
#' Read \code{map} from a CSV file.
#' 
#' @param infile Input CSV file path.
#' @param require.mapunit Require map unit information.
#'     
#' @return An \pkg{R/qtl} \code{map} object.
#' 
#' @export
#' @family csv utilities
#' @rdname readMapCSV
readMapCSV <- function(infile, require.mapunit=TRUE) {
    return( as.map( readMapframeCSV(infile, require.mapunit=require.mapunit) ) )
}

# readMapframeCSV --------------------------------------------------------------
#' Read \code{mapframe} from a CSV file.
#' 
#' @param infile Input CSV file path.
#' @param require.mapunit Require map unit information.
#'     
#' @return A \code{mapframe} object.
#' 
#' @export
#' @family csv utilities
#' @importFrom utils read.csv
#' @rdname readMapframeCSV
readMapframeCSV <- function(infile, require.mapunit=TRUE) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    
    # Read mapframe from CSV file.
    x <- utils::read.csv(infile, check.names=FALSE, quote='', strip.white=TRUE,
        comment.char='', stringsAsFactors=FALSE, colClasses='character',
        na.strings='')
    
    # Validate map unit information.
    map.unit <- getMapUnit(x)
    
    if ( is.na(map.unit) ) {
        
        if (require.mapunit) {
            stop("input map must include map unit information (e.g. 'pos (cM)', '30 kb')")
        }
        
        map.unit <- 'cM'
    }
    
    # Set locus IDs from an input 'id' column, if present.
    if ( 'id' %in% colnames(x) ) {
        x <- setRownamesFromColumn(x, col.name='id')
    }
    
    return( as.mapframe(x, map.unit=map.unit) )
}

# readPhenoCSV -----------------------------------------------------------------
#' Read yeast phenotype data from a CSV file.
#'
#' @param infile Input CSV file path.
#'
#' @export
#' @family csv utilities
#' @importFrom utils read.csv
#' @rdname readPhenoCSV
readPhenoCSV <- function(infile) {
    
    stopifnot( isSingleString(infile) )
    
    # Read phenotype input data as CSV file.
    pheno.table <- utils::read.csv(infile, header=FALSE, check.names=FALSE, quote='',
        stringsAsFactors=FALSE, strip.white=TRUE, colClasses='character',
        na.strings=const$missing.value)
    
    # Trim any blank rows/columns from the bottom/right, respectively.
    pheno.table <- bstripBlankRows( rstripBlankCols(pheno.table) )
    
    # Make pheno table column names from first row of input table.
    colnames(pheno.table) <- make.names(pheno.table[1, ])
    
    # Check for ID heading.
    id.col <- which( tolower( colnames(pheno.table) ) == 'id' )
    if ( length(id.col) == 0 ) {
        stop("ID column not found in phenotype data file - '", infile, "'")
    }
    
    return( as.pheno(pheno.table) )
}

# recodeCSV --------------------------------------------------------------------
#' Recode \pkg{R/qtl} data in CSV file.
#' 
#' @param infile Input CSV file path.
#' @param outfile Output CSV file path.
#' @param geno Named vector of genotype symbols, with vector names containing
#' existing genotype symbols, and vector elements containing their replacement
#' values.
#'
#' @importFrom utils read.csv
#' @importFrom utils write.table
#' @keywords internal
#' @rdname recodeCSV
recodeCSV <- function(infile, outfile, geno=NULL) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    stopifnot( isSingleString(outfile) )
    
    # At least one recoding parameter must be specified.
    if ( is.null(geno) ) {
        stop("cannot recode - no recoding specified")
    }
    
    # Read input CSV file.
    x <- utils::read.csv(infile, header=FALSE, check.names=FALSE, quote='', 
        stringsAsFactors=FALSE, strip.white=TRUE, na.strings=const$missing.value)
    
    # Trim any blank rows/columns from the bottom/right, respectively.
    x <- bstripBlankRows( rstripBlankCols(x) )
    
    params <- getMetadataCSV(x)
    
    if ( ! is.null(geno) ) {
        
        stopifnot( is.character(geno) )
        stopifnot( hasNames(geno) )
        
        src.geno <- names(geno)
        dest.geno <- unname(geno)
        
        if ( length(params$geno.cols) == 0 ) {
            stop("cannot recode - no genotype data")
        }
        
        if ( anyNA(src.geno) || anyNA(dest.geno) ) {
            stop("cannot recode - incomplete genotype recoding")
        }
        
        if ( anyDuplicated(src.geno) || anyDuplicated(dest.geno) ) {
            stop("cannot recode - ambiguous genotype recoding")
        }
        
        # Existing genotypes can be any string,
        # but replacements must be valid genotypes.
        validateGenotypeSet(dest.geno)
        
        # Get genotype data.
        geno.data <- x[params$dat.rows, params$geno.cols]
        
        # Get set of symbols in genotype data.
        curr.geno <- unique( as.character( unlist(geno.data) ) )
        curr.geno <- curr.geno[ ! is.na(curr.geno) ]
        
        err.geno <- curr.geno[ ! curr.geno %in% src.geno ]
        if ( length(err.geno) > 0 ) {
            stop("no recoding defined for genotype symbols - '", toString(err.geno), "'")
        }
        
        # Recode genotypes.
        for ( i in getIndices(geno) ) {
            geno.data[ geno.data == src.geno[i] ] <- dest.geno[i]
        }
        
        # Replace genotype data.
        x[params$dat.rows, params$geno.cols] <- geno.data
    } 
    
    # Write recoded data to file.
    utils::write.table(x, file=outfile, na=const$missing.value, sep=',',
        quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# sniffCSV ---------------------------------------------------------------------
#' Identify type of \pkg{R/qtl} data in CSV file.
#' 
#' @param infile Input CSV file path.
#'
#' @return A string describing the type of data in the input CSV file. This can
#' be \code{'cross'}, \code{'geno'}, \code{'pheno'}, or \code{'map'}. Returns
#' \code{NULL} if the data could not be identified.
#'
#' @export
#' @importFrom utils read.csv
#' @keywords internal
#' @rdname sniffCSV
sniffCSV <- function(infile) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    
    x <- utils::read.csv(infile, header=FALSE, nrows=4, check.names=FALSE, quote='',
        stringsAsFactors=FALSE, strip.white=TRUE, na.strings=const$missing.value)
    
    x <- rstripBlankCols(x)
    
    data.class <- NULL
    
    tryCatch({
        params <- getMetadataCSV(x)
        data.class <- params$class
    }, error=function(e) {})
    
    return(data.class)
}

# writeCrossCSV ----------------------------------------------------------------
#' Write yeast \code{cross} to a CSV file.
#' 
#' @description Function \code{writeCrossCSV} writes a yeast \code{cross} to a 
#' CSV file. Phenotype, genotype, and sample IDs are taken from the \code{'info'}
#' attribute of the \code{cross}, if present.
#'  
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param outfile Output CSV file path.
#' @param chr Vector of sequences for which genotype data should be included in 
#' the output file. If none are specified, genotype data are output for all 
#' sequences.
#' @param digits If specified, round genetic map positions and numeric phenotype 
#' values to the specified number of digits.
#' @param include.mapunit Include map unit information in map positions.
#'  
#' @export
#' @family csv utilities
#' @importFrom utils write.table
#' @rdname writeCrossCSV
writeCrossCSV <- function(cross, outfile, chr=NULL, digits=NULL, 
    include.mapunit=TRUE) {
    
    stopifnot( isSingleString(outfile) )
    stopifnot( isBOOL(include.mapunit) )
    
    map.unit <- 'cM'
    
    # Get indices of phenotype columns.
    pheno.col <- getPhenoColIndices(cross)
    
    # Get CrossInfo object.
    cross.info <- attr(cross, 'info')
    
    # Get specified sequences.
    chr <- subsetBySeq(pull.chr(cross), chr)
    stopifnot( length(chr) > 0 )

    # Take cross info from CrossInfo, if available..
    if ( ! is.null(cross.info) ) {
        compareCrossInfo(cross, cross.info)
        crosstype <- getCrossType(cross.info)
        phenotypes <- getPhenotypes(cross.info)
        alleles <- getAlleles(cross.info)
        markers <- getSeqMarkers(cross.info, normSeq(chr), simplify=TRUE)
        sample.ids <- getSamples(cross.info)
    } else { # ..otherwise take directly from cross.
        crosstype <- pull.crosstype(cross)
        phenotypes <- qtl::phenames(cross)[pheno.col]
        alleles <- pull.alleles(cross)
        markers <- qtl::markernames(cross, chr)
        sample.ids <- pull.ind(cross)
    }
    
    stopifnot( crosstype == 'bc' ) # TODO: support other cross types
    stopifnot( length(phenotypes) > 0 )
    stopifnot( length(alleles) > 0 )
    stopifnot( length(markers) > 0 )
    stopifnot( length(sample.ids) > 0 )
    
    # Get phenotypes, map, genotypes.
    pheno.table <- cross$pheno[, pheno.col, drop=FALSE]
    map.table <- as.data.frame(qtl::pull.map(cross, chr), map.unit=map.unit)
    geno.table <- qtl::pull.geno(cross)[, markers]
    
    # As cross is haploid (though marked as 'bc'),
    # we can set genotypes directly from alleles.
    # TODO: handle other ploidies.
    genotypes <- alleles
    
    # If digits specified, round numeric phenotype values and map positions.
    if ( ! is.null(digits) ) {
        stopifnot( isSinglePositiveWholeNumber(digits) )
        mask <- sapply(pheno.table, is.numeric)
        pheno.table[, mask] <- round(pheno.table[, mask], digits=digits)
        map.table$pos <- round(map.table$pos, digits=digits)
    }
    
    # Replace encoded genotypes with actual genotype values.
    for ( i in getIndices(genotypes) ) {
        geno.table[ geno.table == i ] <- genotypes[i]
    }
    
    # Replace missing genotype data with missing value symbol.
    geno.table[ is.na(geno.table) ] <- const$missing.value
    
    # If writing map unit, add map unit info to map table.
    if (include.mapunit) {
        map.table <- setPosColDataMapUnit(map.table, map.unit)
    }
    
    # Bind map and genotype tables into one output table.
    output.table <- rbind(t(map.table), geno.table)
    
    # Pad phenotypes with two blank rows at top, corresponding 
    # to the location of the map above the genotype data.
    pheno.padding <- data.frame(matrix('', nrow=2, ncol=ncol(pheno.table), 
        dimnames=list(NULL, colnames(pheno.table))), stringsAsFactors=FALSE)
    pheno.table <- rbind(pheno.padding, pheno.table)
    
    # Bind phenotype table to output table.
    output.table <- cbind(pheno.table, output.table)
    output.table <- data.frame( lapply(output.table, as.character) )
    
    # Set output table headings as first row of output table.
    # NB: we want table headings to be written unmodified.
    output.headings <- data.frame( t( c(phenotypes, markers) ), 
        stringsAsFactors=FALSE)
    colnames(output.headings) <- colnames(output.table)
    output.table <- rbind(output.headings, output.table)
    
    # If sample IDs defined, insert sample ID  
    # column between phenotypes and genotypes.
    if ( ! is.null(sample.ids) ) {
        
        # Get ID column name.
        id.col.name <- colnames(cross$pheno)[ getIdColIndex(cross) ]
        
        # Pad sample IDs with column heading and with two blank values, 
        # corresponding to the location of the map above the genotype data.
        padded.ids <- c(id.col.name, rep_len('', 2), sample.ids)
        
        # Set ID column index to first column after phenotypes.
        id.col <- ncol(pheno.table) + 1
        
        # Insert sample ID column into output table.
        output.table <- insertColumn(output.table, col.index=id.col, 
            col.name=id.col.name, data=padded.ids)
    } 
    
    # Write cross data to file.
    utils::write.table(output.table, file=outfile, na=const$missing.value, sep=',',
        quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    return( invisible() )
}

# writeGenoCSV -----------------------------------------------------------------
#' Write yeast genotype data to a CSV file.
#'   
#' @param geno An \pkg{R/qtl} \code{cross} \code{geno} object.
#' @param outfile Output CSV file path.
#' @param chr Vector of sequences for which genotype data should be included in 
#' the output file. If none are specified, a genotype file is output for all 
#' available sequences.
#' @param digits If specified, round genetic map positions to the specified 
#' number of digits.
#' @param include.mapunit Include map unit information in map positions.
#'  
#' @export
#' @family csv utilities
#' @importFrom utils write.table
#' @rdname writeGenoCSV
writeGenoCSV <- function(geno, outfile, chr=NULL, digits=NULL, 
    include.mapunit=TRUE) {
        
    stopifnot( isSingleString(outfile) )
    
    # Convert geno data to a data frame for output.
    geno.table <- as.data.frame(geno, chr=chr, digits=digits, 
        include.mapunit=include.mapunit)
    
    # Check for sample IDs - required in genotype file.
    if ( any(geno.table[4:nrow(geno.table), 'id'] == '') ) {
        stop("cannot output genotype data without sample IDs")
    }
    
    # Write cross geno data to CSV file.
    utils::write.table(geno.table, file=outfile, na=const$missing.value, sep=',',
        quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    return( invisible() )
}

# writeMapCSV ------------------------------------------------------------------
#' Write \code{map} to a CSV file.
#' 
#' @param map An \pkg{R/qtl} \code{map} object.
#' @param outfile Output CSV file path.
#' @param include.mapunit Include map unit information in map positions.
#'  
#' @export
#' @family csv utilities
#' @rdname writeMapCSV
writeMapCSV <- function(map, outfile, include.mapunit=TRUE) {
    stopifnot( 'map' %in% class(map) )
    writeMapframeCSV(as.mapframe(map), outfile, include.mapunit=include.mapunit)
    return( invisible() )
}

# writeMapframeCSV -------------------------------------------------------------
#' Write \code{mapframe} to a CSV file.
#' 
#' @param x A \code{mapframe} object.
#' @param outfile Output CSV file path.
#' @param include.mapunit Include map unit information in map positions.
#'  
#' @export
#' @family csv utilities
#' @importFrom utils write.csv
#' @rdname writeMapframeCSV
writeMapframeCSV <- function(x, outfile, include.mapunit=TRUE) {
    
    stopifnot( 'mapframe' %in% class(x) )
    stopifnot( isSingleString(outfile) )
    stopifnot( isBOOL(include.mapunit) )
    
    x <- as.data.frame(x)
    
    # If including map unit, append a map 
    # unit suffix to each map position.
    if (include.mapunit) {
        x <- setPosColNameMapUnit(x, getMapUnit(x))
    }
    
    # Set 'id' column from locus IDs, if present.
    if ( hasRownames(x) ) {
        x <- setColumnFromRownames(x, col.name='id')
    }
    
    # Write mapframe to CSV file.
    utils::write.csv(x, file=outfile, quote=FALSE, row.names=FALSE)
    
    return( invisible() )
}

# writePhenoCSV ----------------------------------------------------------------
#' Write yeast pheno data to a CSV file.
#'   
#' @param pheno An \pkg{R/qtl} \code{cross} \code{pheno} object.
#' @param outfile Output CSV file path.
#' @param digits If specified, round numeric phenotype values to the specified
#' number of digits.
#'  
#' @export
#' @family csv utilities
#' @importFrom utils write.table
#' @rdname writePhenoCSV
writePhenoCSV <- function(pheno, outfile, digits=NULL) {
    
    stopifnot( isSingleString(outfile) )
    
    # Convert pheno data to a data frame for output.
    pheno.table <- as.data.frame(pheno, digits=digits)
    
    # Write cross pheno data to CSV file.
    utils::write.table(pheno.table, file=outfile, na=const$missing.value, sep=',',
        quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    return( invisible() )
}

# End of csv.R #################################################################