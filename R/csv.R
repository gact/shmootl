# Start of csv.R ###############################################################

# getMetadataCSV -------------------------------------------------------------
#' Get parameters of \pkg{R/qtl} input data.
#' 
#' Given the path to an \pkg{R/qtl} input CSV file, or a \code{data.frame}
#' loaded from such a file, this function determines the parameters of that
#' data (e.g. phenotype columns).
#' 
#' @details The returned list includes information on the following parameters:
#' 
#' \itemize{
#'  \item{\code{pheno.cols}}{Indices of phenotype columns,
#'  or \code{NULL} if none found.}
#'  \item{\code{id.col}}{Index of ID column, or \code{NA} if not found.}
#'  \item{\code{geno.cols}}{Indices of genotype columns,
#'  or \code{NULL} if none found.}
#'  \item{\code{dat.rows}}{Row indices of data.}
#'  \item{\code{map.present}}{\code{TRUE} if the data contains a map;
#'  \code{FALSE} otherwise.}
#'  \item{\code{class}}{Expected class of the input data.
#'  This is \code{'unknown'} by default.}
#' }
#' 
#' This information can then guide further processing of the input data.
#' 
#' @param x An \pkg{R/qtl} input CSV file or
#' an \pkg{R/qtl} input \code{data.frame}.
#'
#' @return A list containing parameters of \pkg{R/qtl} input data.
#' 
#' @keywords internal
#' @rdname getMetadataCSV
getMetadataCSV <- function(x) {
    
    if ( is.data.frame(x) ) {
        
        data.class <- 'unknown'
        
        if ( nrow(x) < 2 ) {
            stop("data not found")
        }
        
        if ( anyNA(x[1, ]) ) {
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
                
                if ( ! id.col %in% pheno.cols ) {
                    stop("ID column not in phenotype columns")
                }
                
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
            
            geno.cols <- NULL
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
                pheno.cols <- getColIndices(x)[-id.col] # every column (except ID column) is a pheno column
            }
        }
        
        first.data.row <- length(head.rows) + 1
        last.data.row <- nrow(x)
        
        if ( last.data.row < first.data.row ) {
            stop("data not found")
        }
        
        dat.rows <- first.data.row : last.data.row
        
        if ( length(pheno.cols) == 0 ) {
            pheno.cols <- NULL
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
        
    } else if ( isSingleString(x) ) {
        
        stopifnot( file.exists(x) )
        
        x <- utils::read.csv(x, header=FALSE, nrows=4, check.names=FALSE, quote='',
            stringsAsFactors=FALSE, strip.white=TRUE, na.strings=const$missing.value)
        
        x <- bstripBlankRows( rstripBlankCols(x) )
        
        params <- getMetadataCSV(x)
        
    } else {
        
        stop("cannot get CSV metadata from object of class - '", class(x), "'")
    }
    
    return(params)
}

# hasMapCSV --------------------------------------------------------------------
#' Test if CSV file contains map data.
#' 
#' @param infile Input CSV file path.
#'
#' @return \code{TRUE} if CSV file contains map data in a known format;
#' \code{FALSE} otherwise.
#'
#' @export
#' @family CSV functions
#' @family map utility functions
#' @rdname hasMapCSV
hasMapCSV <- function(infile) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    
    status <- FALSE
    
    tryCatch({
        params <- getMetadataCSV(infile)
        status <- params$map.present
    }, error=function(e) {})
    
    return(status)
}

# readCovarCSV -----------------------------------------------------------------
#' Read covariate matrix from a CSV file.
#' 
#' This function reads a table of covariates (with optional \code{'ID'} column)
#' from a CSV file and returns a numeric \code{matrix} of covariate data. Sample
#' IDs are not included in the returned matrix, as \pkg{R/qtl} matches the
#' covariate matrix row-by-row to the sample rows of the \code{cross} object.
#' The contents of the input table do not need to be numeric, but they will be
#' converted to numeric values when being read from file. For example, columns
#' containing character values are read as factors, which are converted to their
#' numeric representation.
#' 
#' @param infile Input CSV file path.
#' 
#' @return A numeric \code{matrix} of covariate data.
#' 
#' @export
#' @family CSV functions
#' @importFrom utils read.csv
#' @importFrom utils type.convert
#' @rdname readCovarCSV
readCovarCSV <- function(infile) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    
    # Read covariate data as character data-frame from CSV file.
    covar.table <- utils::read.csv(infile, header=TRUE, check.names=FALSE,
        quote='', stringsAsFactors=TRUE, strip.white=TRUE,
        colClasses='character', na.strings=const$missing.value)
    stopif( anyDuplicated( colnames(covar.table) ) )
    
    # Trim any blank rows/columns from the bottom/right, respectively.
    covar.table <- bstripBlankRows( rstripBlankCols(covar.table) )
    
    # If 'id' column present, remove it from covariate table.
    id.col <- which( tolower( colnames(covar.table) ) == 'id' )
    if ( length(id.col) > 0 ) {
        covar.table <- deleteColumn(covar.table, col.index=id.col)
    }
    
    # Convert columns of covariate data-frame
    # from character to appropriate datatypes.
    covar.table <- do.call(cbind.data.frame,
        lapply(covar.table, utils::type.convert))
    stopifnot( all( sapply(covar.table, class) %in%
        c('double', 'factor', 'integer', 'logical', 'numeric') ) )
    
    # Convert covariate table to numeric matrix.
    # NB: converts factors to their numeric encoding.
    covar.matrix <- data.matrix(covar.table, rownames.force=FALSE)
    
    return(covar.matrix)
}

# readCrossCSV -----------------------------------------------------------------
#' Read yeast \code{cross} from a CSV file.
#' 
#' This function reads yeast cross data from an \pkg{R/qtl} CSV file and
#' returns an \pkg{R/qtl} \code{cross} object. The returned \code{cross}
#' has an attribute \code{'info'} of type \code{\linkS4class{CrossInfo}},
#' which contains information about the input \code{cross}.
#' 
#' If a genetic map is not present in the input file, this is estimated when
#' loading the cross. In such cases, parameters \code{error.prob} and
#' \code{map.function} are passed to \pkg{R/qtl} for map estimation.
#' 
#' By default, any genetic map positions must include map units (i.e.
#' \code{'cM'}). To read a cross without requiring map units, set
#' parameter \code{require.mapunit} to \code{FALSE}. If map units
#' are not required and not found, map positions are assumed to be
#' in centiMorgans.
#' 
#' @param infile Input CSV file path.
#' @param error.prob Genotyping error rate.
#' @param map.function Genetic map function.
#' @param require.mapunit Require map unit information in map positions.
#'  
#' @return An \pkg{R/qtl} \code{cross} object with an attribute
#' \code{'info'} of type \code{\linkS4class{CrossInfo}}.
#' 
#' @export
#' @family CSV functions
#' @family cross object functions
#' @importFrom methods new
#' @importFrom stats na.omit
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
    
    # Get vector of non-missing values in phenotype data.
    pheno.values <- stats::na.omit( unlist(cross.table[dat.rows, pheno.cols]) )
    
    # Verify that there are no blank phenotype values.
    if ( any( pheno.values == '' ) ) {
        stop("blank phenotype values found")
    }
    
    # Get set of unique symbols used in marker genotype data.
    geno.symbols <- sort( unique( as.character( # NB: sort removes NA values
        unlist(cross.table[dat.rows, geno.cols]) ) ) )
    
    # Verify that there are no blank genotype entries.
    if ( any( geno.symbols == '' ) ) {
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
    cross.info <- setCrosstype(cross.info, crosstype)
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
#' This function reads yeast cross genotype data from an \pkg{R/qtl} CSV file
#' and returns a \code{geno} object. The returned \code{geno} object has an
#' attribute \code{'info'} of type \code{\linkS4class{CrossInfo}}. This
#' \code{\linkS4class{CrossInfo}} object contains information about the
#' genotype data, including sample IDs, which must be present in the
#' input genotype data file.
#' 
#' By default, any genetic map positions must include map units (i.e.
#' \code{'cM'}). To read a genotype file without requiring map units,
#' set parameter \code{require.mapunit} to \code{FALSE}. If map units
#' are not required and not found, map positions are assumed to be in
#' centiMorgans.
#' 
#' @param infile Input CSV file path.
#' @param require.mapunit Require map unit information in map positions.
#' 
#' @export
#' @family CSV functions
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
    
    # Normalise sequence IDs.
    geno.table[2, -id.col] <- normSeq( as.character(geno.table[2, -id.col]) )
    
    return( as.geno(geno.table, require.mapunit=require.mapunit) )
}

# readMapCSV -------------------------------------------------------------------
#' Read \code{map} from a CSV file.
#' 
#' This function reads a yeast cross genetic or physical map from any
#' \pkg{R/qtl} CSV file that contains map data, whether that file
#' contains a map table, genotype data, or a full cross.
#' 
#' By default, the input map must include map unit information (e.g.
#' \code{'cM'}, \code{'bp'}). In the case of a map table file, this information
#' can be indicated in the map headings (e.g. \code{'pos (cM)'}) or in the map
#' positions (e.g. \code{'47 cM'}). For a cross or genotype file, map units
#' should be included with map positions.
#' 
#' To read an input file without requiring map units, set parameter
#' \code{require.mapunit} to \code{FALSE}. If map units are not
#' required and not found, map positions are assumed to be in
#' centiMorgans.
#' 
#' @param infile Input CSV file path.
#' @param require.mapunit Require map unit information.
#'     
#' @return An \pkg{R/qtl} \code{map} object.
#' 
#' @export
#' @family CSV functions
#' @family map utility functions
#' @rdname readMapCSV
readMapCSV <- function(infile, require.mapunit=TRUE) {
    
    stopifnot( isBOOL(require.mapunit) )
    
    params <- getMetadataCSV(infile)
    
    if ( ! params$map.present ) {
        stop("no map data found in file - '", infile,"'")
    }
    
    if ( params$class == 'cross' ) {
        
        cross <- readCrossCSV(infile, require.mapunit=require.mapunit)
        cross.map <- qtl::pull.map(cross)
        
    } else if ( params$class == 'geno' ) {
        
        geno <- readGenoCSV(infile, require.mapunit=require.mapunit)
        cross.map <- pullMap(geno)
        
    } else if ( params$class == 'map' ) {
        
        cross.map <- as.map( readMapframeCSV(infile,
            require.mapunit=require.mapunit) )
        
    } else {
        
        stop("cannot read map from ", params$class, " data file - '", infile,"'")
    }
    
    return(cross.map)
}

# readMapframeCSV --------------------------------------------------------------
#' Read \code{mapframe} from a CSV file.
#' 
#' This function reads a yeast cross genetic or physical mapframe from any
#' \pkg{R/qtl} CSV file that contains map data, whether that file contains
#' a map table, genotype data, or a full cross.
#' 
#' By default, the input mapframe must include map unit information (e.g.
#' \code{'cM'}, \code{'bp'}). In the case of a map table file, this information
#' can be indicated in the map headings (e.g. \code{'pos (cM)'}) or in the map
#' positions (e.g. \code{'47 cM'}). For a cross or genotype file, map units
#' should be included with map positions.
#' 
#' To read an input file without requiring map units, set parameter
#' \code{require.mapunit} to \code{FALSE}. If map units are not
#' required and not found, map positions are assumed to be in
#' centiMorgans.
#' 
#' @param infile Input CSV file path.
#' @param require.mapunit Require map unit information.
#'     
#' @return A \code{mapframe} object.
#' 
#' @export
#' @family CSV functions
#' @family map utility functions
#' @importFrom utils read.csv
#' @rdname readMapframeCSV
readMapframeCSV <- function(infile, require.mapunit=TRUE) {
    
    stopifnot( isBOOL(require.mapunit) )
    
    params <- getMetadataCSV(infile)
    
    if ( ! params$map.present ) {
        stop("no map data found in file - '", infile,"'")
    }
    
    if ( params$class %in% c('cross', 'geno') ) {
        
        cross.map <- as.mapframe( readMapCSV(infile,
            require.mapunit=require.mapunit) )
        
    } else if ( params$class == 'map' ) {
        
        # Read mapframe from CSV file.
        x <- utils::read.csv(infile, check.names=FALSE, quote='',
            strip.white=TRUE, comment.char='', stringsAsFactors=FALSE,
            colClasses='character', na.strings='')
        
        # Trim any blank rows/columns from the bottom/right, respectively.
        x <- bstripBlankRows( rstripBlankCols(x) )
        
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
        
        cross.map <- as.mapframe(x, map.unit=map.unit)
        
    } else {
        
        stop("cannot read mapframe from ", params$class, " data file - '", infile,"'")
    }
    
    return(cross.map)
}

# readPhenoCSV -----------------------------------------------------------------
#' Read yeast phenotype data from a CSV file.
#'
#' This function reads yeast cross phenotype data from an \pkg{R/qtl} CSV file
#' and returns a \code{pheno} object. The returned \code{pheno} object has an
#' attribute \code{'info'} of type \code{\linkS4class{CrossInfo}}. This
#' \code{\linkS4class{CrossInfo}} object contains information about the
#' phenotype data, including sample IDs, which must be present in the
#' input phenotype data file.
#' 
#' @param infile Input CSV file path.
#'
#' @export
#' @family CSV functions
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
#' This function recodes the data in an \pkg{R/qtl} cross or genotype CSV file.
#' Genotypes can be recoded by passing to the \code{geno} parameter a simple
#' \code{\link{mapping}} of old to new genotypes. Alternatively, genotype data
#' can be converted to enumerated genotypes by setting parameter
#' \code{enum.geno} to \code{TRUE}.
#' 
#' Any \code{\link{mapping}} of old to new symbols must be complete and
#' unambiguous. All existing values must be mapped to a new value, and
#' data with different input values cannot have the same output value.
#' 
#' @param infile Input CSV file path.
#' @param outfile Output CSV file path.
#' @param geno A simple \code{\link{mapping}} of old to new genotype symbols,
#' with names containing existing genotype symbols, and elements containing
#' their replacement values (incompatible with \code{enum.geno}).
#' @param enum.geno Option indicating if genotype data should be recoded as
#' enumerated genotypes (incompatible with \code{geno}).
#'
#' @export
#' @family CSV functions
#' @importFrom utils read.csv
#' @importFrom utils write.table
#' @rdname recodeCSV
recodeCSV <- function(infile, outfile, geno=NULL, enum.geno=FALSE) {
    
    # TODO: implement phenotype recoding
    # TODO: implement sample recoding ?
    # TODO: implement marker recoding ?
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    stopifnot( isSingleString(outfile) )
    stopifnot( isBOOL(enum.geno) )
    
    # Check if genotype recoding specified.
    recode.geno <- ! is.null(geno) && length(geno) > 0
    
    # At least one recoding parameter must be specified.
    if ( ! ( recode.geno || enum.geno ) ) {
        stop("cannot recode - no recoding option specified")
    }
    
    # Check no more than one genotype recoding option specified.
    if ( recode.geno && enum.geno ) {
        stop("multiple genotype recoding options specified")
    }
    
    # Read input CSV file.
    x <- utils::read.csv(infile, header=FALSE, check.names=FALSE, quote='', 
        stringsAsFactors=FALSE, strip.white=TRUE, na.strings=const$missing.value)
    
    # Trim any blank rows/columns from the bottom/right, respectively.
    x <- bstripBlankRows( rstripBlankCols(x) )
    
    params <- getMetadataCSV(x)
    
    if (recode.geno) {
        
        stopifnot( is.mapping(geno) )
        
        src.geno <- names(geno)
        dest.geno <- unlist( values(geno) )
        
        if ( is.null(params$geno.cols) ) {
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
        recoded.geno <- original.geno <- x[params$dat.rows, params$geno.cols]
        
        # Get set of symbols in genotype data.
        curr.geno <- unique( as.character( unlist(original.geno) ) )
        curr.geno <- curr.geno[ ! is.na(curr.geno) ]
        
        err.geno <- curr.geno[ ! curr.geno %in% src.geno ]
        if ( length(err.geno) > 0 ) {
            stop("no recoding defined for genotype symbols - '", toString(err.geno), "'")
        }
        
        # Recode genotypes.
        for ( i in getIndices(geno) ) {
            recoded.geno[ original.geno == src.geno[i] ] <- dest.geno[i]
        }
        
        # Replace genotype data.
        x[params$dat.rows, params$geno.cols] <- recoded.geno
        
    } else if (enum.geno) {
        
        # Get original genotype data.
        original.geno <- x[params$dat.rows, params$geno.cols]
        
        # Init recoded genotype data.
        recoded.geno <- matrix(NA_integer_, nrow=nrow(original.geno),
            ncol=ncol(original.geno), dimnames=dimnames(original.geno) )
        
        # Recode each column independently.
        for ( i in getColIndices(original.geno) ) {
            
            # Get sample symbols and genotypes for this locus.
            locus.symbols <- original.geno[, i]
            locus.levels <- unique(locus.symbols)
            locus.genotypes <- locus.levels[ locus.levels != const$missing.value ]
            
            # Skip loci with more genotypes than can be represented.
            if ( length(locus.genotypes) > length(const$enum.geno.charset) ) {
                next
            }
            
            # Set genotype numbers for this locus.
            recoded.geno[, i] <- match(locus.symbols, locus.genotypes)
        }
        
        # Get number of null loci in recoded genotype data.
        null.count <- sum( apply(recoded.geno, 2, function(column)
            allNA( as.vector(column) ) ) )
        
        if ( null.count > 0 ) {
            warning("enumerated genotype data contains ", null.count, " null loci")
        }
        
        # Replace genotype data.
        x[params$dat.rows, params$geno.cols] <- recoded.geno
    }
    
    # Write recoded data to file.
    utils::write.table(x, file=outfile, na=const$missing.value, sep=',',
        quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# sniffCSV ---------------------------------------------------------------------
#' Identify type of \pkg{R/qtl} data in CSV file.
#' 
#' This function identifies the type of data in the input CSV file. This can
#' be \code{'cross'}, \code{'geno'}, \code{'pheno'}, or \code{'map'}.
#' Alternatively, it can be \code{'unknown'} if the data cannot be identified.
#' 
#' Note that this function should be used with some discretion, as some
#' \pkg{R/qtl} input file types cannot be distinguished. For example, a
#' covariate data file containing a sample ID column will be misidentified
#' as containing phenotype data.
#' 
#' @param infile Input CSV file path.
#'
#' @return A string describing the type of data in the input CSV file.
#'
#' @export
#' @family CSV functions
#' @rdname sniffCSV
sniffCSV <- function(infile) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    
    data.class <- 'unknown'
    
    tryCatch({
        params <- getMetadataCSV(infile)
        data.class <- params$class
    }, error=function(e) {})
    
    return(data.class)
}

# writeCrossCSV ----------------------------------------------------------------
#' Write yeast \code{cross} to a CSV file.
#' 
#' This function writes a yeast \code{cross} to an \pkg{R/qtl} CSV file.
#' Phenotype, genotype, and sample IDs are taken from the \code{'info'}
#' attribute of the \code{cross}, if present.
#'  
#' @param cross An \pkg{R/qtl} \code{cross} object.
#' @param outfile Output CSV file path.
#' @param chr Vector of sequences for which genotype data should be included in 
#' the output file. If none are specified, genotype data are output for all 
#' sequences.
#' @param digits If specified, round genetic map positions and
#' numeric phenotype values to the specified number of digits.
#' @param include.mapunit Include map unit information in map positions.
#'  
#' @export
#' @family CSV functions
#' @family cross object functions
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
        crosstype <- getCrosstype(cross.info)
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
#' This function writes a yeast \code{geno} object to an \pkg{R/qtl} CSV file.
#' The \code{geno} object must have an attribute \code{'info'} of type
#' \code{\linkS4class{CrossInfo}}, from which genotype and sample IDs are taken.
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
#' @family CSV functions
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
#' This function writes a yeast \code{mapframe} object to any \pkg{R/qtl} CSV
#' file that can contain map data, whether it is a map table file, genotype
#' file, or a full cross file. If an output genotype or cross file already
#' exists, the mapframe is written to the map portion of the output file.
#' 
#' @param map An \pkg{R/qtl} \code{map} object.
#' @param outfile Output CSV file path.
#' @param include.mapunit Include map unit information in output file.
#'  
#' @export
#' @family CSV functions
#' @family map utility functions
#' @rdname writeMapCSV
writeMapCSV <- function(map, outfile, include.mapunit=TRUE) {
    
    stopifnot( 'map' %in% class(map) )
    stopifnot( isSingleString(outfile) )
    stopifnot( isBOOL(include.mapunit) )
    
    # Assume that we are not pushing map into output file.
    pushing <- FALSE
    
    # If output file exists and is a cross or geno CSV
    # file, we are pushing map into output file.
    if ( file.exists(outfile) ) {
        
        guess <- sniffCSV(outfile)
        
        if ( guess %in% c('cross', 'geno') ) {
            pushing <- TRUE
        } else if ( guess != 'map' ) {
            stop("cannot write map to existing file with ", guess, " data - '", outfile,"'")
        }
    }
    
    # If pushing map, push into cross or geno object and write result to file..
    if (pushing) {
        
        # Create output temp file.
        tmp <- tempfile()
        on.exit( file.remove(tmp) )
        
        if ( guess == 'cross' ) {
            
            cross <- readCrossCSV(outfile, require.mapunit=FALSE)
            cross <- qtl::replace.map(cross, map)
            writeCrossCSV(cross, tmp, include.mapunit=include.mapunit)
            
        } else { # guess == 'geno'
            
            geno <- readGenoCSV(outfile, require.mapunit=FALSE)
            geno <- pushMap(geno, map)
            writeGenoCSV(geno, tmp, include.mapunit=include.mapunit)
        }
        
        # Move temp file to final output file.
        # NB: file.copy is used here instead of file.rename because the latter
        # can sometimes fail when moving files between different file systems.
        file.copy(tmp, outfile, overwrite=TRUE)
        
    } else { # ..otherwise convert to mapframe and write to file.
        
        writeMapframeCSV(as.mapframe(map), outfile,
            include.mapunit=include.mapunit)
    }
    
    return( invisible() )
}

# writeMapframeCSV -------------------------------------------------------------
#' Write \code{mapframe} to a CSV file.
#' 
#' This function writes a yeast \code{mapframe} object to any \pkg{R/qtl} CSV
#' file that can contain map data, whether it is a map table file, genotype
#' file, or a full cross file. If an output genotype or cross file already
#' exists, the mapframe is written to the map portion of the output file.
#' 
#' @param x A \code{mapframe} object.
#' @param outfile Output CSV file path.
#' @param include.mapunit Include map unit information in output file.
#'  
#' @export
#' @family CSV functions
#' @family map utility functions
#' @importFrom utils write.csv
#' @rdname writeMapframeCSV
writeMapframeCSV <- function(x, outfile, include.mapunit=TRUE) {
    
    stopifnot( 'mapframe' %in% class(x) )
    stopifnot( isSingleString(outfile) )
    stopifnot( isBOOL(include.mapunit) )
    
    # Assume that we are not pushing mapframe into output file.
    pushing <- FALSE
    
    # If output file exists and is a cross or geno CSV
    # file, we are pushing mapframe into output file.
    if ( file.exists(outfile) ) {
        
        guess <- sniffCSV(outfile)
        
        if ( guess %in% c('cross', 'geno') ) {
            pushing <- TRUE
        } else if ( guess != 'map' ) {
            stop("cannot write map to existing file with ", guess, " data - '", outfile,"'")
        }
    }
    
    # If pushing mapframe, convert to map and push into output file..
    if (pushing) {
        
        writeMapCSV(as.map(x), outfile, include.mapunit=include.mapunit)
        
    } else { # ..otherwise write mapframe to file.
        
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
    }
    
    return( invisible() )
}

# writePhenoCSV ----------------------------------------------------------------
#' Write yeast pheno data to a CSV file.
#' 
#' This function writes a yeast \code{pheno} object to an \pkg{R/qtl} CSV file.
#' The \code{pheno} object must have an attribute \code{'info'} of type
#' \code{\linkS4class{CrossInfo}}, from which phenotype and sample IDs are taken.
#' 
#' @param pheno An \pkg{R/qtl} \code{cross} \code{pheno} object.
#' @param outfile Output CSV file path.
#' @param digits If specified, round numeric phenotype
#' values to the specified number of digits.
#' 
#' @export
#' @family CSV functions
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