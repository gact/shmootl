# Start of excel.R #############################################################

# listWorksheets ---------------------------------------------------------------
#' List available output worksheets.
#' 
#' @return Character vector containing the names of worksheets that can be set
#' as output to a spreadsheet file in \pkg{shmootl} spreadsheet functions. This
#' excludes summary worksheets such as \code{'README'} and \code{'Overview'},
#' as these are included in every output spreadsheet file.
#' 
#' @export
#' @family Excel functions
#' @rdname listWorksheets
listWorksheets <- function() {
    return( const$result.worksheets )
}

# makeWorksheetNames -----------------------------------------------------------
#' Make \pkg{shmootl} output worksheet names.
#' 
#' @param x A list of character vectors in which the names of the list are
#' analyses, and the elements of each character vector are results for the
#' given analysis.
#' 
#' @return Vector of \pkg{shmootl} output worksheet names.
#' 
#' @keywords internal
#' @rdname makeWorksheetNames
makeWorksheetNames <- function(x) {
    
    stopifnot( is.list(x) )
    
    worksheets <- character()
    
    if ( length(x) > 0 ) {
        
        stopifnot( length(x) > 0 )
        stopifnot( all( sapply(x, is.character) ) )
        stopifnot( all( lengths(x) > 0 ) )
        stopifnot( all( names(x) %in% names(const$supported.analyses) ) )
        
        worksheets <- unname( unlist( lapply( names(x), function(analysis)
            sapply( x[[analysis]], function(result) paste(analysis, result) ) ) ) )
        
        stopifnot( all( worksheets %in% const$result.worksheets ) )
    }
    
    return(worksheets)
}

# parseWorksheetNames ----------------------------------------------------------
#' Parse \pkg{shmootl} output worksheet names.
#' 
#' @param worksheets Vector of \pkg{shmootl} output worksheet names.
#' 
#' @return A list of character vectors in which the names of the list are
#' analyses, and the elements of each character vector are results for the
#' given analysis.
#' 
#' @keywords internal
#' @rdname parseWorksheetNames
parseWorksheetNames <- function(worksheets) {
    
    stopifnot( is.character(worksheets) )
    
    parsed <- list()
    
    if ( length(worksheets) > 0 ) {
        
        stopifnot( all( isValidID(worksheets) ) )
        stopifnot( all( worksheets %in% const$result.worksheets ) )
        
        m <- regexec(const$pattern$worksheet, worksheets)
        regmatch.list <- regmatches(worksheets, m)
        
        analyses <- sapply(regmatch.list, getElement, 2)
        results <- sapply(regmatch.list, getElement, 3)
        
        aset <- names(const$supported.analyses)[
            names(const$supported.analyses) %in% analyses ]
        
        for ( analysis in aset ) {
            parsed[[analysis]] <- results[ analyses == analysis ]
        }
    }
    
    return(parsed)
}

# writeDigestExcel -------------------------------------------------------------
#' Write an Excel digest of QTL scan results.
#' 
#' @param digest Path of output digest Excel file.
#' @inheritParams writeWorkbookExcel
#' 
#' @template author-thomas-walsh
#' @template author-yue-hu
#' 
#' @export
#' @family digest functions
#' @family Excel functions
#' @rdname writeDigestExcel
writeDigestExcel <- function(scanfiles, digest, phenotypes=NULL,
    analyses=NULL, worksheets=NULL, scanfile.pattern=NULL) {
    
    if ( is.null(worksheets) ) {
        worksheets <- const$default.worksheets[['digest']]
    }
    
    writeWorkbookExcel(scanfiles, digest, phenotypes=phenotypes,
        analyses=analyses, worksheets=worksheets,
        scanfile.pattern=scanfile.pattern)
    
    return( invisible() )
}

# writeReportExcel -------------------------------------------------------------
#' Write an Excel report of QTL scan results.
#' 
#' @param scanfile A QTL scan result HDF5 file.
#' @param report Path of output report Excel file.
#' @inheritParams writeWorkbookExcel
#' 
#' @template author-yue-hu
#' @template author-thomas-walsh
#' 
#' @export
#' @family report functions
#' @family Excel functions
#' @rdname writeReportExcel
writeReportExcel <- function(scanfile, report, phenotypes=NULL,
    analyses=NULL, worksheets=NULL, scanfile.pattern=NULL) {
    
    stopifnot( isSingleString(scanfile) )
    
    if ( is.null(worksheets) ) {
        worksheets <- const$default.worksheets[['report']]
    }
    
    writeWorkbookExcel(scanfile, report, phenotypes=phenotypes,
        analyses=analyses, worksheets=worksheets,
        scanfile.pattern=scanfile.pattern)
    
    return( invisible() )
}

# writeWorkbookExcel -----------------------------------------------------------
#' Write an Excel workbook of QTL analysis results.
#' 
#' @param scanfiles One or more QTL scan result HDF5 files.
#' @param workbook Path of output Excel file.
#' @param phenotypes Phenotypes for which results should be included in the
#' output file. If none are specified, results are output for all phenotypes.
#' @param analyses Analyses for which results should be included in the output
#' file. If none are specified, results are output for all available analyses.
#' @param worksheets Worksheets to include in the output file. If none are
#' specified, default worksheets are output. To view available worksheets,
#' call the function \code{\link{listWorksheets}}.
#' @param scanfile.pattern Pattern for extracting experiment info from scan file
#' names. This must be a valid Perl regex with named capture groups. Neither the
#' capture groups nor the pattern itself are required to match any given scan
#' file, but all capture groups must have a name, and that name cannot clash
#' with other names that might be used alongside the extracted information.
#' 
#' @importFrom utils installed.packages
#' @keywords internal
#' @rdname writeWorkbookExcel
writeWorkbookExcel <- function(scanfiles, workbook, phenotypes=NULL,
    analyses=NULL, worksheets=NULL, scanfile.pattern=NULL) {
    
    if ( ! 'xlsx' %in% rownames(utils::installed.packages()) ) {
        stop("cannot write Excel workbook without R package 'xlsx'")
    }
    
    stopifnot( is.character(scanfiles) )
    stopifnot( length(scanfiles) > 0 )
    stopif( anyDuplicated(scanfiles) )
    stopifnot( all( file.exists(scanfiles) ) )
    stopifnot( isSingleString(workbook) )
    workbook.format <- inferFormatFromFilename(workbook)
    stopifnot( workbook.format == 'Excel' )
    
    # If scanfile pattern specified, get experiment info from scan file names.
    if ( ! is.null(scanfile.pattern) ) {
        xinfo <- getXInfoFromFilenames(scanfiles, scanfile.pattern)
    } else {
        xinfo <- NULL
    }
    
    # Attach required package namespaces (if needed) ---------------------------
    
    req.pkgs <- c('rJava', 'xlsxjars', 'xlsx')
    names(req.pkgs) <- paste0('package:', req.pkgs)
    att.pkgs <- req.pkgs[ ! names(req.pkgs) %in% search() ]
    sapply(att.pkgs, attachNamespace)
    on.exit( sapply(names(att.pkgs), detach, character.only=TRUE), add=TRUE )
    
    # Check scan files for results of interest ---------------------------------
    
    # Set possible results to be sought in scan file.
    results.sought <- supported.results <- list(
        'Scanone' = c('QTL Intervals'),
        'Scantwo' = c('QTL Pairs')
    )
    
    # If analyses specified, filter results sought by given analyses.
    if ( ! is.null(analyses) ) {
        
        analyses <- unique( resolveAnalysisTitle(analyses) )
        
        unsupported <- analyses[ ! analyses %in% names(supported.results) ]
        if ( length(unsupported) > 0 ) {
            stop("Excel output not supported for analyses - '",
                toString(unsupported), "'")
        }
        
        results.sought <- results.sought[ names(results.sought) %in% analyses ]
    }
    
    # If worksheets specified, filter results sought for given worksheets.
    if ( ! is.null(worksheets) ) {
        
        unknown <- worksheets[ ! worksheets %in% const$result.worksheets ]
        if ( length(unknown) > 0 ) {
            stop("cannot output unknown Excel worksheets - '",
                toString(unknown), "'")
        }
        
        supported.worksheets <- makeWorksheetNames(supported.results)
        unsupported <- worksheets[ ! worksheets %in% supported.worksheets ]
        if ( length(unsupported) > 0 ) {
            stop("Excel output not supported for worksheets - '",
                toString(unsupported), "'")
        }
        
        sheets.sought <- makeWorksheetNames(results.sought)
        sheets.sought <- sheets.sought[ sheets.sought %in% worksheets ]
        results.sought <- parseWorksheetNames(sheets.sought)
    }
    
    if ( length(results.sought) == 0 ) {
        stop("cannot search for results with the given search parameters")
    }
    
    # Get result info.
    rinfo <- getResultInfoHDF5(scanfiles, phenotypes=phenotypes,
        analyses=names(results.sought))
    
    # Get results of interest.
    roi <- list()
    for ( scanfile in scanfiles ) {
        for ( phenotype in names(rinfo[[scanfile]]) ) {
            for ( analysis in names(rinfo[[scanfile]][[phenotype]]) ) {
                for ( result in rinfo[[scanfile]][[phenotype]][[analysis]] ) {
                    if ( analysis %in% names(results.sought) &&
                        result %in% results.sought[[analysis]] ) {
                        roi[[analysis]] <- union(roi[[analysis]], result)
                    }
                }
            }
        }
    }
    
    if ( length(roi) == 0 ) {
        stop("no relevant results found with the given search parameters")
    }
    
    # Set worksheet names ------------------------------------------------------
    
    sheet.names <- c(const$summary.worksheets, makeWorksheetNames(roi))
    
    # Init worksheet tables ----------------------------------------------------
    
    tables <- vector('list', length(sheet.names))
    names(tables) <- sheet.names
    
    # Setup 'README' worksheet -------------------------------------------------
    
    descriptions <- unlist( lapply(sheet.names, function(k)
        const$excel[[k]]$description) )
    
    tables[['README']] <- data.frame(Worksheet=sheet.names,
        Description=descriptions, stringsAsFactors=FALSE)
    
    # Setup 'Overview' worksheet -----------------------------------------------
    
    headings <- const$excel[['Overview']]$headings
    
    # Init exhaustive results overview (i.e. including all possible columns).
    tables[['Overview']] <- as.data.frame( matrix(nrow=0,
        ncol=length(headings), dimnames=list(NULL, headings) ),
        stringsAsFactors=FALSE)
    
    # Get exhaustive results overview of HDF5 scan files.
    for ( scanfile in scanfiles ) {
        
        input.oview <- readResultsOverviewHDF5(scanfile)
        
        output.oview <- as.data.frame( matrix(nrow=nrow(input.oview),
            ncol=length(headings), dimnames=list(NULL, headings) ),
            stringsAsFactors=FALSE)
        
        for ( k in colnames(input.oview) ) {
            output.oview[, k] <- as.character(input.oview[, k])
        }
        
        output.oview[, 'File'] <- scanfile
        
        tables[['Overview']] <- rbind(tables[['Overview']], output.oview)
    }
    
    # Remove unused results overview columns.
    tables[['Overview']] <- removeColsNA(tables[['Overview']])
    
    # Setup 'Scanone QTL Intervals' worksheet, if relevant ---------------------
    
    if ( 'Scanone QTL Intervals' %in% sheet.names ) {
        
        headings <- const$excel[['Scanone QTL Intervals']]$headings
        
        # Init table with all headings; empty columns will be deleted later.
        tab <- matrix(NA_character_, nrow=0, ncol=length(headings),
            dimnames=list(NULL, headings) )
        
        for ( scanfile in scanfiles ) {
            
            for ( phenotype in names(rinfo[[scanfile]]) ) {
                
                if ( 'Scanone' %in% names(rinfo[[scanfile]][[phenotype]]) &&
                    'QTL Intervals' %in% rinfo[[scanfile]][[phenotype]][['Scanone']] ) {
                    
                    qtl.intervals <- readResultHDF5(scanfile,
                        phenotype, 'Scanone', 'QTL Intervals')
                    num.qtls <- length(qtl.intervals)
                    
                    if ( num.qtls == 0 ) {
                        next
                    }
                    
                    attrs <- attributes(qtl.intervals)
                    
                    # Get any threshold attributes for this set of QTL intervals.
                    threshold <- if ( 'threshold' %in% names(attrs) ) { attrs$threshold } else { NA_real_ }
                    alpha <- if ( 'alpha' %in% names(attrs) ) { attrs$alpha } else { NA_real_ }
                    fdr <- if ( 'fdr' %in% names(attrs) ) { attrs$fdr } else { NA_real_ }
                    
                    # If possible, get interval type from associated attributes.
                    if ( 'prob' %in% names(attrs) ) {
                        interval.type <- paste0(attrs$prob * 100, '% Bayes Credible Interval')
                    } else if ( 'drop' %in% names(attrs) ) {
                        interval.type <- paste0(attrs$drop, '-LOD Support Interval')
                    } else {
                        interval.type <- NA_character_
                    }
                    
                    # Check if QTL intervals object contains physical positions.
                    physical.positions <- hasPhysicalPositions(qtl.intervals)
                    
                    # Get QTL interval features.
                    qtl.features <- readResultHDF5(scanfile, phenotype, 'Scanone', 'QTL Features')
                    
                    for ( i in seq_along(qtl.intervals) ) {
                        
                        # Get QTL interval.
                        qtl.interval <- qtl.intervals[[i]]
                        
                        # Get name of this QTL interval.
                        qtl.name <- names(qtl.intervals)[i]
                        
                        # Get reference sequence for this QTL interval.
                        interval.seq <- unique(qtl.interval[, 'chr'])
                        
                        # Get peak LOD value in this QTL interval.
                        peak.LOD <- qtl.interval[2, 'lod']
                        
                        # Get genetic map positions of this QTL interval.
                        start.cM <- qtl.interval[1, 'pos (cM)']
                        peak.cM <- qtl.interval[2, 'pos (cM)']
                        end.cM <- qtl.interval[3, 'pos (cM)']
                        
                        # Get physical map positions of QTL interval, if available.
                        if (physical.positions) {
                            start.bp <- qtl.interval[1, 'pos (bp)']
                            peak.bp <- qtl.interval[2, 'pos (bp)']
                            end.bp <- qtl.interval[3, 'pos (bp)']
                        } else {
                            start.bp <- peak.bp <- end.bp <- NA_integer_
                        }
                        
                        # Get QTL interval features, if available.
                        if ( ! is.null(qtl.features) && qtl.name %in% names(qtl.features) ) {
                            feature.ids <- paste(qtl.features[[qtl.name]]$ID, collapse=' ')
                        } else {
                            feature.ids <- NA_character_
                        }
                        
                        # Set row of info for this QTL interval.
                        row <- matrix( c(
                            scanfile,      # File
                            phenotype,     # Phenotype
                            qtl.name,      # QTL Name
                            interval.seq,  # Chromosome
                            peak.LOD,      # Peak LOD
                            threshold,     # LOD Threshold
                            alpha,         # alpha
                            fdr,           # FDR
                            interval.type, # Interval Type
                            start.cM,      # Start (cM)
                            peak.cM,       # Peak (cM)
                            end.cM,        # End (cM)
                            start.bp,      # Start (bp)
                            peak.bp,       # Peak (bp)
                            end.bp,        # End (bp)
                            feature.ids    # Scanone QTL Features
                        ), nrow=1, ncol=length(headings),
                        dimnames=list(NULL, colnames(tab)))
                        
                        # Add row to table.
                        tab <- rbind(tab, row)
                    }
                }
            }
        }
        
        # Remove empty columns.
        tab <- removeColsNA(tab)
        
        # Add experiment info, if available.
        if ( ! is.null(xinfo) ) {
            tab <- addXInfo(tab, xinfo)
        }
        
        # Set QTL intervals table.
        tables[['Scanone QTL Intervals']] <- data.frame(tab, check.names=FALSE,
            stringsAsFactors=FALSE)
    }
    
    # Setup 'Scantwo QTL Pairs’ worksheet, if relevant -------------------------
    
    if ( 'Scantwo QTL Pairs' %in% sheet.names ) {
        
        headings <- const$excel[['Scantwo QTL Pairs']]$headings
        
        # Init table with all headings; empty columns will be deleted later.
        tab <- matrix(NA_character_, nrow=0, ncol=length(headings),
            dimnames=list(NULL, headings) )
        
        for ( scanfile in scanfiles ) {
            
            for ( phenotype in names(rinfo[[scanfile]]) ) {
                
                if ( 'Scantwo' %in% names(rinfo[[scanfile]][[phenotype]]) &&
                    'QTL Pairs' %in% rinfo[[scanfile]][[phenotype]][['Scantwo']] ) {
                    
                    qtl.pairs <- readResultHDF5(scanfile,
                        phenotype, 'Scantwo', 'QTL Pairs')
                    num.pairs <- nrow(qtl.pairs)
                    
                    if ( num.pairs == 0 ) {
                        next
                    }
                    
                    qtl.pairs <- as.matrix(qtl.pairs, rownames.force=FALSE)
                    qtl.pairs <- insertColumn(qtl.pairs, col.index=1,
                        col.name='File', data=scanfile)
                    qtl.pairs <- insertColumn(qtl.pairs, col.index=2,
                        col.name='Phenotype', data=phenotype)
                    colnames(qtl.pairs) <- headings
                    
                    # Add row to table.
                    tab <- rbind(tab, qtl.pairs)
                }
            }
        }
        
        # Remove empty columns.
        tab <- removeColsNA(tab)
        
        # Add experiment info, if available.
        if ( ! is.null(xinfo) ) {
            tab <- addXInfo(tab, xinfo)
        }
        
        # Set QTL Pairs table.
        tables[['Scantwo QTL Pairs']] <- data.frame(tab, check.names=FALSE,
            stringsAsFactors=FALSE)
    }
    
    # Output Excel workbook ----------------------------------------------------
    
    # Create workbook, setting format from file extension.
    wb <- xlsx::createWorkbook( type=tools::file_ext(workbook) )
    
    # Setup table heading style.
    heading.style <- xlsx::CellStyle(wb) +
        xlsx::Alignment(horizontal='ALIGN_CENTER') +
        xlsx::Font(wb, isBold=TRUE)
    
    # Add each worksheet to workbook.
    for ( i in seq_along(tables) ) {
        ws <- xlsx::createSheet(wb, sheetName=sheet.names[i])
        stopifnot( is.data.frame(tables[[i]]) )
        xlsx::addDataFrame(tables[[i]], ws, col.names=TRUE,
            row.names=FALSE, colnamesStyle=heading.style, showNA=FALSE)
        xlsx::autoSizeColumn(ws, getColIndices(tables[[i]]))
    }
    
    # Save workbook to file.
    xlsx::saveWorkbook(wb, workbook)
    
    return( invisible() )
}

# End of excel.R ###############################################################