# Start of excel.R #############################################################

# getExcelDigestInfo -----------------------------------------------------------
#' Get info for digest Excel file.
#' 
#' @return A list containing descriptions
#' and headings for Excel digest worksheets.
#' 
#' @keywords internal
#' @rdname getExcelDigestInfo
getExcelDigestInfo <- function() {
    
    info <- list(
        
        `README` = list(
        
            description = 'This describes the contents of every worksheet in this workbook.',
        
            headings = c('Worksheet', 'Description')
        ),
            
        `Overview` = list(
            
            description = 'Results overview from across the set of scan files.',
            
            headings = c('File', 'Phenotype', getPkgAnalysisNames())
        ),
        
        `Scanone QTL Intervals` = list(
            
            description = paste(
                'Table of QTL intervals as obtained by a single- or multi-QTL scan.',
                'Genomic features within the QTL interval are included, if available.'
            ),
            
            headings = c('File', 'Phenotype', 'QTL Name', 'Chromosome',
                'Peak LOD', 'LOD Threshold', 'alpha', 'FDR',
                'Interval Type', 'Start (cM)', 'Peak (cM)', 'End (cM)',
                'Start (bp)', 'Peak (bp)', 'End (bp)', 'Scanone QTL Features')
        )
    )
    
    return(info)
}

# writeDigestExcel -------------------------------------------------------------
#' Write an Excel digest of QTL scan results.
#' 
#' @param scanfiles One or more QTL scan result HDF5 files.
#' @param digest Path of output digest Excel file.
#' @param scanfile.pattern Optional pattern for extracting experiment info from
#' scan file names. This must be a valid Perl regex with named capture groups.
#' Neither the capture groups nor the pattern itself are required to match any
#' given scan file, but all capture groups must have a name. Each such name is
#' used in the digest file, as a heading for a column that contains matches to
#' that capture group in the scan file names.
#' 
#' @template author-thomas-walsh
#' @template author-yue-hu
#' 
#' @export
#' @family digest functions
#' @family Excel functions
#' @importFrom utils installed.packages
#' @rdname writeDigestExcel
writeDigestExcel <- function(scanfiles, digest, scanfile.pattern=NULL) {
    
    if ( ! 'xlsx' %in% rownames(utils::installed.packages()) ) {
        stop("cannot write Excel digest without R package 'xlsx'")
    }
    
    stopifnot( is.character(scanfiles) )
    stopifnot( length(scanfiles) > 0 )
    stopif( anyDuplicated(scanfiles) )
    stopifnot( all( file.exists(scanfiles) ) )
    stopifnot( isSingleString(digest) )
    stopifnot( tools::file_ext(digest) %in% const$ext$excel )
    
    # Get Excel digest info.
    digest.info <- getExcelDigestInfo()
    
    # If scanfile pattern specified, get scanfile info from scan file names.
    scanfile.info <- NULL
    if ( ! is.null(scanfile.pattern) ) {
        
        stopifnot( isSingleString(scanfile.pattern) )
        
        # Parse scan file names by the given pattern.
        parsed <- parseFilenames(scanfiles, scanfile.pattern)
        
        # Set reserved headings that cannot be used as capture-group names.
        # NB: this should include headings for all
        # tables that will include scan file info.
        reserved.headings <- digest.info[['Scanone QTL Intervals']]$headings
        
        # Check for capture-group names clashing with reserved headings.
        clashing <- colnames(parsed)[ colnames(parsed) %in% reserved.headings ]
        if ( length(clashing) > 0 ) {
            stop("scan file capture-group names clash with headings - '",
                toString(clashing), "'")
        }
        
        # Remove empty columns.
        nonempty <- sapply( getColIndices(parsed),
            function(i) ! allNA(parsed[, i]) )
        parsed <- parsed[, nonempty]
        
        if ( ncol(parsed) > 0 ) {
            scanfile.info <- parsed
        }
    }
    
    # Attach required package namespaces (if needed) ---------------------------
    
    req.pkgs <- c('rJava', 'xlsxjars', 'xlsx')
    names(req.pkgs) <- paste0('package:', req.pkgs)
    att.pkgs <- req.pkgs[ ! names(req.pkgs) %in% search() ]
    sapply(att.pkgs, attachNamespace)
    on.exit( sapply(names(att.pkgs), detach, character.only=TRUE), add=TRUE )
    
    # Check scan files for results of interest ---------------------------------
    
    results.sought <- c('Scanone/QTL Intervals', 'Scanone/QTL Features')
    result.info <- list()
    roi <- character()
    
    for ( scanfile in scanfiles ) {
        
        if ( ! hasResultsOverviewHDF5(scanfile) ) {
            stop("cannot create digest - result overview not found in file '",
                scanfile, "'")
        }
        
        # Get result info for this HDF5 scan file.
        info <- list()
        for ( phenotype in getResultPhenotypesHDF5(scanfile) ) {
            analyses <- getResultAnalysesHDF5(scanfile, phenotype)
            info[[phenotype]] <- lapply( analyses, function(analysis)
                getResultNamesHDF5(scanfile, phenotype, analysis) )
            names(info[[phenotype]]) <- analyses
        }
        
        # Get results of interest.
        for ( phenotype in names(info) ) {
            for ( analysis in names(info[[phenotype]]) ) {
                for ( result in info[[phenotype]][[analysis]] ) {
                    h5name <- joinH5ObjectNameParts(c(analysis, result), relative=TRUE)
                    if ( h5name %in% results.sought ) {
                        roi <- union(roi, h5name)
                    }
                }
            }
        }
        
        result.info[[scanfile]] <- info
    }
    
    # Set worksheet names ------------------------------------------------------
    
    sheet.names <- c('README', 'Overview')
    
    if ( 'Scanone/QTL Intervals' %in% roi ) {
        sheet.names <- c(sheet.names, 'Scanone QTL Intervals')
    }
    
    # Init worksheet tables ----------------------------------------------------
    
    tables <- vector('list', length(sheet.names))
    names(tables) <- sheet.names
    
    # Setup 'README' worksheet -------------------------------------------------
    
    descriptions <- unlist( lapply(sheet.names, function(k)
        digest.info[[k]]$description) )
    
    tables[['README']] <- data.frame(Worksheet=sheet.names,
        Description=descriptions, stringsAsFactors=FALSE)
    
    # Setup 'Overview' worksheet -----------------------------------------------
    
    headings <- digest.info[['Overview']]$headings
    
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
    nonempty <- sapply( getColIndices(tables[['Overview']]),
        function(i) ! allNA(tables[['Overview']][, i]) )
    tables[['Overview']] <- tables[['Overview']][, nonempty]
    
    # Setup 'Scanone QTL Intervals' worksheet, if relevant ---------------------
    
    if ( 'Scanone QTL Intervals' %in% sheet.names ) {
        
        headings <- digest.info[['Scanone QTL Intervals']]$headings
        
        # Init table with all headings; empty columns will be deleted later.
        tab <- matrix(NA_character_, nrow=0, ncol=length(headings),
            dimnames=list(NULL, headings) )
        
        for ( scanfile in scanfiles ) {
            
            for ( phenotype in names(result.info[[scanfile]]) ) {
                
                if ( 'Scanone' %in% names(result.info[[scanfile]][[phenotype]]) &&
                    'QTL Intervals' %in% result.info[[scanfile]][[phenotype]][['Scanone']] ) {
                    
                    qtl.intervals <- readResultHDF5(scanfile,
                        phenotype, 'Scanone', 'QTL Intervals')
                    
                    if ( length(qtl.intervals) == 0 ) {
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
        nonempty <- sapply( getColIndices(tab),
            function(i) ! allNA(tab[, i]) )
        tab <- tab[, nonempty, drop=FALSE]
        
        # Add scanfile info, if available.
        if ( ! is.null(scanfile.info) ) {
            
            row.indices <- sapply( getRowIndices(tab), function(i)
                which( rownames(scanfile.info) == tab[i, 'File'] ) )
            
            scanfile.table <- scanfile.info[row.indices, ]
            rownames(scanfile.table) <- NULL
            
            tab <- cbind(tab[, 1, drop=FALSE], scanfile.table,
                tab[, -1, drop=FALSE])
        }
        
        # Set QTL intervals table.
        tables[['Scanone QTL Intervals']] <- data.frame(tab, check.names=FALSE,
            stringsAsFactors=FALSE)
    }
    
    # Output Excel workbook ----------------------------------------------------
    
    # Create workbook, setting format from file extension.
    workbook <- xlsx::createWorkbook( type=tools::file_ext(digest) )
    
    # Setup table heading style.
    heading.style <- xlsx::CellStyle(workbook) + xlsx::Alignment(horizontal='ALIGN_CENTER') +
        xlsx::Font(workbook, isBold=TRUE)
    
    # Add each worksheet to workbook.
    for ( i in seq_along(tables) ) {
        worksheet <- xlsx::createSheet(workbook, sheetName=sheet.names[i])
        stopifnot( is.data.frame(tables[[i]]) )
        xlsx::addDataFrame(tables[[i]], worksheet, col.names=TRUE,
            row.names=FALSE, colnamesStyle=heading.style, showNA=FALSE)
        xlsx::autoSizeColumn(worksheet, getColIndices(tables[[i]]))
    }
    
    # Save workbook to file.
    xlsx::saveWorkbook(workbook, digest)
    
    return( invisible() )
}

# End of excel.R ###############################################################