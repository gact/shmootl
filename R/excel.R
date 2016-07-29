# Start of excel.R #############################################################

# writeDigestExcel -------------------------------------------------------------
#' Write an Excel digest of QTL scan results.
#' 
#' @param scanfiles One or more QTL scan result HDF5 files.
#' @param digest Path of output digest Excel file.
#' 
#' @export
#' @importFrom tools file_ext
#' @importFrom xlsx addDataFrame
#' @importFrom xlsx Alignment
#' @importFrom xlsx autoSizeColumn
#' @importFrom xlsx CellStyle
#' @importFrom xlsx createSheet
#' @importFrom xlsx createWorkbook
#' @importFrom xlsx Font
#' @importFrom xlsx saveWorkbook
#' @include const.R
#' @include hdf5.R
#' @rdname writeDigestExcel
writeDigestExcel <- function(scanfiles, digest) {
    
    stopifnot( is.character(scanfiles) )
    stopifnot( length(scanfiles) > 0 )
    stopifnot( all( file.exists(scanfiles) ) )
    stopifnot( isSingleString(digest) )
    stopifnot( tools::file_ext(digest) %in% const$ext$excel )
    
    # Attach required package namespaces ---------------------------------------
    
    req.pkgs <- c('rJava', 'xlsxjars', 'xlsx')
    sapply(req.pkgs, attachNamespace)
    on.exit( sapply(paste0('package:', req.pkgs), detach, character.only=TRUE) )
    
    # Check scan files for results of interest ---------------------------------
    
    results.sought <- 'QTL Intervals'
    
    # Init per-scanfile info.
    result.info <- results.found <- pheno.info <- vector('list', length(scanfiles))
    names(result.info) <- names(results.found) <- names(pheno.info) <- scanfiles
    
    for ( scanfile in scanfiles ) {
        
        if ( ! hasObjectHDF5(scanfile, 'Results/Overview') ) {
            stop("cannot create digest - result overview not found in file '",
                scanfile, "'")
        }
        
        # Get phenotypes from scan result file.
        pheno.info[[scanfile]] <- getPhenotypesHDF5(scanfile)
        
        # Get all result names for each phenotype.
        result.info[[scanfile]] <- lapply(pheno.info[[scanfile]], function(phenotype)
            getResultNamesHDF5(scanfile, phenotype))
        names(result.info[[scanfile]]) <- pheno.info[[scanfile]]
        
        # Get names of results of interest for each phenotype.
        results.found[[scanfile]] <- lapply(result.info[[scanfile]], function(results)
            results[ results %in% results.sought ])
    }
    
    # Get set of results of interest.
    roi <- unique( unlist(results.found) )
    
    # Set worksheet names ------------------------------------------------------
    
    sheet.names <- c('README', 'Overview')
    
    if ( 'QTL Intervals' %in% roi ) {
        sheet.names <- c(sheet.names, 'QTL Intervals')
    }
    
    # Init worksheet tables ----------------------------------------------------
    
    tables <- vector('list', length(sheet.names))
    names(tables) <- sheet.names
    
    # Setup 'README' worksheet -------------------------------------------------
    
    descriptions <- unlist( lapply(sheet.names, function(k)
        const$excel$digest[[k]]$description) )
    
    tables[['README']] <- data.frame(Worksheet=sheet.names,
        Description=descriptions, stringsAsFactors=FALSE)
               
    # Setup 'Overview' worksheet -----------------------------------------------
    
    headings <- const$excel$digest[['Overview']]$headings
    
    tables[['Overview']] <- as.data.frame( matrix(nrow=0, ncol=length(headings),
        dimnames=list(NULL, headings) ) )
    
    for ( scanfile in scanfiles ) {
        overview <- readOverviewHDF5(scanfile)
        overview <- insertColumn(overview, 1, col.name='File', data=scanfile)
        tables[['Overview']] <- rbind(tables[['Overview']], overview)
    }
    
    # Setup 'QTL Intervals' worksheet, if relevant -----------------------------
    
    if ( 'QTL Intervals' %in% sheet.names ) {
        
        headings <- const$excel$digest[['QTL Intervals']]$headings
        
        # Init table with all headings; empty columns will be deleted later.
        tab <- matrix(NA_character_, nrow=0, ncol=length(headings),
            dimnames=list(NULL, headings) )
        
        for ( scanfile in scanfiles ) {
            
            phenotypes <- pheno.info[[scanfile]]
            
            for ( phenotype in phenotypes ) {
                
                if ( 'QTL Intervals' %in% results.found[[scanfile]][[phenotype]] ) {
                    
                    qtl.intervals <- readResultHDF5(scanfile, phenotype, 'QTL Intervals')
                    
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
                    
                    for ( qtl.interval in qtl.intervals ) {
                        
                        # Get reference sequence for this QTL interval.
                        interval.seq <- unique(qtl.interval[, 'chr'])
                        
                        # Set row of info for this QTL interval.
                        row <- matrix( c(
                            scanfile,                    # File
                            phenotype,                   # Phenotype
                            interval.seq,                # Chromosome
                            qtl.interval[2, 'lod'],      # Peak LOD
                            threshold,                   # LOD Threshold
                            alpha,                       # alpha
                            fdr,                         # FDR
                            interval.type,               # Interval Type
                            qtl.interval[1, 'pos (cM)'], # Start (cM)
                            qtl.interval[2, 'pos (cM)'], # Peak (cM)
                            qtl.interval[3, 'pos (cM)'], # End (cM)
                            NA_integer_,                 # Start (bp)
                            NA_integer_,                 # End (bp)
                            NA_character_                # Features
                        ), nrow=1, ncol=length(headings),
                        dimnames=list(NULL, colnames(tab)))
                        
                        # Add row to table.
                        tab <- rbind(tab, row)
                    }
                }
            }
        }
        
        # Set table from nonempty columns.
        nonempty <- sapply( getColIndices(tab), function(i) ! allNA(tab[, i]) )
        tables[['QTL Intervals']] <- data.frame(tab[, nonempty],
            check.names=FALSE, stringsAsFactors=FALSE)
    }
    
    # Output Excel workbook ----------------------------------------------------
    
    # Create workbook, setting format from file extension.
    workbook <- createWorkbook( type=tools::file_ext(digest) )
    
    # Setup table heading style.
    heading.style <- CellStyle(workbook) + Alignment(horizontal='ALIGN_CENTER') +
        Font(workbook, isBold=TRUE)
    
    # Add each worksheet to workbook.
    for ( i in seq_along(tables) ) {
        worksheet <- createSheet(workbook, sheetName=sheet.names[i])
        stopifnot( is.data.frame(tables[[i]]) )
        addDataFrame(tables[[i]], worksheet, col.names=TRUE,
            row.names=FALSE,  colnamesStyle=heading.style, showNA=FALSE)
        autoSizeColumn(worksheet, getColIndices(tables[[i]]))
    }
    
    # Save workbook to file.
    saveWorkbook(workbook, digest)
    
    return( invisible() )
}

# End of excel.R ###############################################################