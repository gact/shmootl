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

# readCrossCSV -----------------------------------------------------------------
#' Read yeast \code{cross} from a CSV file.
#' 
#' This function reads yeast cross data from an \pkg{R/qtl} CSV file and returns 
#' an \pkg{R/qtl} \code{cross} object with an attribute \code{'info'} of type 
#' \code{CrossInfo}. 
#' 
#' @param infile Input CSV file path.
#' @param missing.value Missing data value. This can be any single character 
#' that is not a possible phenotype or genotype value.
#' @param error.prob Genotyping error rate (ignored unless estimating genetic map).
#' @param map.function Genetic map function (ignored unless estimating genetic map).
#'  
#' @return An \pkg{R/qtl} \code{cross} object with an attribute \code{'info'} of 
#' type \code{CrossInfo}.
#' 
#' @export
#' @family csv utilities
#' @importFrom methods new
#' @rdname readCrossCSV
readCrossCSV <- function(infile, missing.value='-', error.prob=0.0001, 
    map.function=c('haldane', 'kosambi', 'c-f', 'morgan')) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    stopifnot( isSingleChar(missing.value) )
    error.prob <- as.numeric(error.prob)
    stopifnot( isSingleProbability(error.prob) )
    
    map.function <- match.arg(map.function)
    
    # Read cross input data as CSV file. Don't check names now, will   
    # check them soon. Don't allow headers, we want to see them as 
    # they are. Replace any whitespace/empty cells with NA values.
    cross.table <- read.csv(infile, header=FALSE, check.names=FALSE, quote='', 
        stringsAsFactors=FALSE, strip.white=TRUE, na.strings=missing.value)
    
    # Get logical vector indicating which columns are blank 
    # in the first row after the initial heading row.
    seq.is.blank <- cross.table[2, ] == ''
    
    # Set phenotype columns from those with blank sequence row.   
    pheno.cols <- which( seq.is.blank )
    
    # Set marker columns from those with nonempty sequence row.   
    geno.cols <- which( ! seq.is.blank )
    
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
    
    # Trim any empty rows from the bottom.
    while ( allNA( cross.table[last.data.row, ] ) ) {
        last.data.row <- last.data.row - 1
    }
    
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
    genotypes <- alleles <- geno.symbols[ ! geno.symbols %in% missing.value ]
    
    # Verify that there are exactly two genotypes.
    # TODO: handle cross with more than two genotypes.
    if ( length(genotypes) != 2 ) { 
        stop("unsupported number of genotypes - '", length(genotypes), "'")
    } 
    
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
        samples <- 1:length(dat.rows)
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
    cross.info <- setSamples(cross.info, samples)
    
    # Get normalised marker sequences, replace these in original table.
    locus.seqs <- getMarkerSeqs(cross.info)
    cross.table[2, geno.cols] <- locus.seqs
    
    # Set sorted sequences in CrossInfo object.
    sequences <- sortSeq( unique(locus.seqs) )
    cross.info <- setSequences(cross.info, sequences)
    
    # If map is present, check for map unit info.
    if (map.present) {
        
        # Set map table from cross map data.
        map.table <- data.frame( row.names=locus.ids, chr=locus.seqs,
            pos=as.character(cross.table[3, geno.cols]),
            stringsAsFactors=FALSE)
        
        # Validate genetic map unit information.
        tryCatch({
            validateGeneticMapUnit(map.table)
        }, error=function(e) {
            stop("cross input map positions must include genetic map units ",
                "(e.g. '47 cM')")
        })
        
        # Strip map unit information.
        map.table <- setPosColDataMapUnit(map.table, NULL)
        
        # Replace original map positions.
        cross.table[3, geno.cols] <- map.table$pos
    }
    
    # Create temp file for adjusted cross data.
    temp.file <- tempfile(fileext='csv')
    write.table(cross.table, file=temp.file, na='', sep=',', 
        quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    # Read adjusted cross data.
    cross <- qtl::read.cross('csv', '', temp.file, genotypes=genotypes, 
        alleles=alleles, na.strings=missing.value, error.prob=error.prob, 
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

# readMapCSV -------------------------------------------------------------------
#' Read \code{map} from a CSV file.
#' 
#' @param infile Input CSV file path.
#'     
#' @return An \pkg{R/qtl} \code{map} object.
#' 
#' @export
#' @family csv utilities
#' @rdname readMapCSV
readMapCSV <- function(infile) {
    return( as.map( readMapframeCSV(infile) ) )
}

# readMapframeCSV --------------------------------------------------------------
#' Read \code{mapframe} from a CSV file.
#' 
#' @param infile Input CSV file path.
#'     
#' @return A \code{mapframe} object.
#' 
#' @export
#' @family csv utilities
#' @rdname readMapframeCSV
readMapframeCSV <- function(infile) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    
    # Read mapframe from CSV file.
    map.table <- read.csv(infile, check.names=FALSE, quote='', strip.white=TRUE, 
        comment.char='', stringsAsFactors=FALSE, na.strings='')
    
    # Validate map unit information.
    tryCatch({
        validateMapUnit(x)
    }, error=function(e) {
        stop("input map positions must include map units ",
            "(e.g. '47 cM', '30 kb')")
    })
    
    # Set locus IDs from an input 'id' column, if present.
    if ( 'id' %in% colnames(x) ) {
        x <- setRownamesFromColumn(x, col.name='id')
    }
    
    return( as.mapframe(x) )
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
#' @param missing.value Missing data value. This can be any single character that 
#' is not a possible phenotype or genotype value.
#' @param include.mapunit Include map unit information in map positions.
#'  
#' @export
#' @family csv utilities
#' @rdname writeCrossCSV
writeCrossCSV <- function(cross, outfile, chr=NULL, digits=NULL, 
    missing.value='-', include.mapunit=TRUE) {
    
    stopifnot( isSingleString(outfile) )
    stopifnot( isSingleChar(missing.value) )
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
        compareCrossInfo(cross.info, cross)
        phenotypes <- getPhenotypes(cross.info)
        alleles <- getAlleles(cross.info)
        markers <- getSeqMarkers(cross.info, normSeq(chr), simplify=TRUE)
        marker.names <- getMarkerNames(cross.info, markers)
        sample.ids <- getSamples(cross.info)
    } else { # ..otherwise take directly from cross.
        phenotypes <- qtl::phenames(cross)[pheno.col]
        alleles <- pull.alleles(cross)
        markers <- marker.names <- qtl::markernames(cross, chr)
        sample.ids <- pull.ind(cross)
    }

    # Get phenotypes, map, genotypes.
    pheno.table <- qtl::pull.pheno(cross, pheno.col)
    map.table <- as.data.frame(qtl::pull.map(cross, chr), map.unit=map.unit)
    geno.table <- qtl::pull.geno(cross)[, marker.names]
    
    # If digits specified, round numeric phenotype values and map positions.
    if ( ! is.null(digits) ) {
        stopifnot( isSinglePositiveWholeNumber(digits) )
        mask <- sapply(pheno.table, is.numeric)
        pheno.table[, mask] <- round(pheno.table[, mask], digits=digits)
        map.table$pos <- round(map.table$pos, digits=digits)
    }
    
    # Replace encoded genotypes with actual genotype values.
    for ( i in 1:length(alleles) ) {
        geno.table[ geno.table == i ] <- alleles[i]
    }
    
    # Replace missing genotype data with missing value symbol.
    geno.table[ is.na(geno.table) ] <- missing.value
    
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
    write.table(output.table, file=outfile, na=missing.value, sep=',', 
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
#' @param missing.value Missing data value. This can be any single character 
#' that is not a possible genotype value.
#' @param include.mapunit Include map unit information in map positions.
#'  
#' @export
#' @family csv utilities
#' @rdname writeGenoCSV
writeGenoCSV <- function(geno, outfile, chr=NULL, digits=NULL, 
    missing.value='-', include.mapunit=TRUE) {
        
    stopifnot( isSingleString(outfile) )
    stopifnot( isSingleChar(missing.value) )
    
    # Convert geno data to a data frame for output.
    geno.table <- makeGenoTable(geno, chr=chr, digits=digits, 
        include.mapunit=include.mapunit)
    
    # Write cross geno data to CSV file.
    write.table(geno.table, file=outfile, na=missing.value, sep=',', 
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
    write.csv(x, file=outfile, quote=FALSE, row.names=FALSE)
    
    return( invisible() )
}

# End of csv.R #################################################################