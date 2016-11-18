# Start of run_scantwo.R #######################################################

# run_scantwo ------------------------------------------------------------------
#' Do 2-QTL scan.
#' 
#' Read cross data from the specified cross input file, run 2-dimensional
#' QTL analysis using \pkg{R/qtl} \code{scantwo} (Broman \emph{et al.} 
#' 2003), and write the results of that scan to the specified output file.
#' 
#' If the input cross contains enumerated genotypes, marker regression
#' is performed regardless of the value of the \code{method} parameter.
#' 
#' In typical usage, LOD threshold stringency can be set by the significance
#' level (\code{alpha}), with a default value of \code{0.05}. The significance
#' level is used to estimate LOD thresholds from \code{scantwo} permutations.
#' 
#' This can be bypassed by setting six fixed LOD \code{thresholds}, along with
#' a nominal significance level (\code{alpha}), in which case permutations are
#' skipped and the fixed LOD thresholds are applied directly for assessing
#' significance. When calling this function from the command line using
#' \code{Rscript}, the \code{thresholds} parameter must be specified as a YAML
#' string mapping \code{scantwo} LOD type names to threshold values (e.g.
#' \code{"full: 8, fv1: 7, int: 6, add: 5, av1: 4, one: 3"}). When called from
#' within the \code{R} environment, the \code{thresholds} parameter must be
#' specified as a mapping object (e.g.
#' \code{mapping( c(full=8, fv1=7, int=6, add=5, av1=4, one=3) )}). 
#' 
#' @param infile input cross file
#' @param h5file scan result file
#' @param chr sequences [default: all]
#' @param pheno phenotypes [default: all]
#' @param model phenotype model
#' @param method method of QTL analysis
#' @param n.perm number of permutations
#' @param n.cluster number of threads
#' @param alpha significance level for LOD thresholds
#' @param thresholds fixed LOD thresholds
#' @param step step size for genotype probabilities
#' @param map.function genetic map function
#' @param error.prob genotyping error rate
#' @param acovfile additive covariates file
#' @param icovfile interactive covariates file
#' 
#' @template author-yue-hu
#' @template author-thomas-walsh
#' @template ref-broman-2003
#' @template seealso-rqtl-manual
#' 
#' @concept shmootl:analysis
#' @export
#' @family pipeline functions
#' @rdname run_scantwo
run_scantwo <- function(infile=NA_character_, h5file=NA_character_,
    chr=character(), pheno=character(), model=c('normal', 'binary'),
    method=c('em', 'imp', 'hk', 'mr', 'mr-imp', 'mr-argmax'), n.perm=1000L,
    n.cluster=1L, alpha=NA_real_, thresholds=mapping(), step=0, 
    map.function=c('haldane', 'kosambi', 'c-f', 'morgan'), error.prob=0.0001,
    acovfile=NA_character_, icovfile=NA_character_) {
    
    # TODO: infer more precise single-QTL intervals with two-QTL scan results,
    # by using the conditional-interactive and conditional-additive models to
    # distinguish linked QTLs where possible.
    
    pipeline <- getPipelineName( as.character(match.call()[[1]]) )
    analysis <- getAnalysisTitle(pipeline)
    
    stopifnot( isSingleString(h5file) )
    stopifnot( isSingleNonNegativeNumber(step) )
    stopifnot( isSingleProbability(error.prob) )
    
    model <- match.arg(model)
    method <- match.arg(method)
    map.function <- match.arg(map.function)
    
    infile <- if ( ! identical(infile, NA_character_) ) { infile } else { NULL }
    chr <- if ( ! identical(chr, character()) ) { chr } else { NULL }
    pheno <- if ( ! identical(pheno, character()) ) { pheno } else { NULL }
    alpha <- if ( ! identical(alpha, NA_real_) ) { alpha } else { NULL }
    thresholds <- if ( ! identical(thresholds, mapping()) ) { thresholds } else { NULL }
    
    if ( ! is.null(alpha) ) {
        stopifnot( isSingleProbability(alpha) )
    } else if ( ! is.null(thresholds) ) {
        stop("cannot set fixed LOD thresholds without nominal significance level (alpha)")
    } else {
        alpha <- 0.05
    }
    
    # Attach required package namespaces (if needed) ---------------------------
    
    req.pkgs <- 'qtl'
    names(req.pkgs) <- paste0('package:', req.pkgs)
    att.pkgs <- req.pkgs[ ! names(req.pkgs) %in% search() ]
    sapply(att.pkgs, attachNamespace)
    on.exit( sapply(names(att.pkgs), detach, character.only=TRUE), add=TRUE )
    
    # --------------------------------------------------------------------------
    
    # Create temp output file, ensure will be removed.
    tmp <- tempfile()
    on.exit( file.remove(tmp), add=TRUE )
    
    # Read cross from input CSV file or HDF5 file as appropriate.
    # If HDF5 file does not exist, create it and add cross object.
    # If both input files exist, their cross objects must match.
    if ( file.exists(h5file) ) {
        updating.h5file <- TRUE
        stopifnot( hasCrossHDF5(h5file) )
        cross <- h5file.cross <- readCrossHDF5(h5file)
        if ( ! is.null(infile) ) {
            infile.cross <- readCrossCSV(infile, error.prob=error.prob,
                map.function=map.function)
            stopifnot( crossesEqual(h5file.cross, infile.cross) )
        }
    } else if ( ! is.null(infile) ) {
        updating.h5file <- FALSE
        cross <- readCrossCSV(infile, error.prob=error.prob,
            map.function=map.function)
        writeCrossHDF5(cross, tmp)
    } else {
        stop("no input CSV file or HDF5 file specified")
    }
    
    # Read additive covariate matrix file, if specified.
    if ( ! is.na(acovfile) ) {
        addcovar <- readCovarCSV(acovfile, cross=cross)
    } else {
        addcovar <- NULL
    }
    
    # Read interactive covariate matrix file, if specified.
    if ( ! is.na(icovfile) ) {
        intcovar <- readCovarCSV(icovfile, cross=cross)
    } else {
        intcovar <- NULL
    }
    
    # Get cross info.
    cross.info <- attr(cross, 'info')
    sequences <- getSequences(cross.info, chr)
    phenotypes <- getPhenotypes(cross.info, pheno)
    pheno.names <- getPhenotypeNames(cross.info, pheno)
    pheno.col <- getPhenoColIndices(cross, pheno.names)
    genotypes <- getGenotypes(cross.info)
    
    # If cross contains enumerated genotypes, do marker regression.
    if ( hasEnumGenotypes(cross) ) {
        
        if ( sum( qtl::nmissing(cross) ) > 0 ) {
            stop("cannot scan with enumerated genotypes - missing genotype data")
        }
        
        cat(" --- Using marker regression on enumerated genotypes\n")
        
        method <- 'mr'
    }
    
    # Calculate probabilities of true underlying genotypes given observed marker data.
    cross <- qtl::calc.genoprob(cross, step=step, error.prob=error.prob, 
        map.function=map.function)
    
    # Run two-dimensional QTL analysis using R/qtl.
    scantwo.result <- batchPhenoScantwo(cross, chr=sequences,
        pheno.col=pheno.col, model=model, method=method, n.cluster=n.cluster,
        addcovar=addcovar, intcovar=intcovar)
    
    # If no fixed thresholds specified, estimate from permutations..
    if( is.null(thresholds) ) {
        
        # Run permutation scans.
        scantwo.perms <- batchPermScantwo(cross, chr=sequences,
            pheno.col=pheno.col, model=model, method=method, n.perm=n.perm,
            n.cluster=n.cluster, addcovar=addcovar, intcovar=intcovar)
        
       # Get LOD thresholds from permutation results. 
       scantwo.thresholds <- summary(scantwo.perms, alpha=alpha)
       
    } else { # ..otherwise make scantwo LOD threshold object.
        
        scantwo.thresholds <- makeScantwoThresholdObject(thresholds,
            phenotypes, alpha)
        scantwo.perms <- NULL
    }
    
    scantwo.overview <- structure(rep('', length(phenotypes)), names=phenotypes)
    
    # If updating HDF5 file, init HDF5 object names of updated objects.
    if (updating.h5file) {
        updated.objects <- character()
    }
    
    # Output results of 2-QTL analysis for each phenotype.
    for ( i in seq_along(phenotypes) ) {
        
        # Output scantwo result for this phenotype.
        # NB: phenotype count check needed because qtl:::subset.scantwo
        # fails when subsetting from a single-phenotype scantwo object.
        if ( length(phenotypes) > 1 ) {
            pheno.result <- subset(scantwo.result, lodcolumn=i) # qtl:::subset.scantwo
            attr(pheno.result, 'phenotypes') <- attr(scantwo.result, 'phenotypes')[i]
        } else {
            pheno.result <- scantwo.result
        }
        writeResultHDF5(pheno.result, tmp, phenotypes[i], analysis, 'Result')
        
        # Output any permutation scan results for this phenotype.
        if ( ! is.null(scantwo.perms) ) {
            pheno.perms <- subset(scantwo.perms, lodcolumn=i) # qtl:::subset.scantwoperm
            writeResultHDF5(pheno.perms, tmp, phenotypes[i], analysis, 'Permutations')
        }
        
        # Output scantwo threshold for this phenotype.
        scantwo.threshold <- subset(scantwo.thresholds, lodcolumns=i)
        writeResultHDF5(scantwo.threshold, tmp, phenotypes[i], analysis, 'Thresholds')
        
        # Get scantwo threshold values for this phenotype.
        thresholds <- sapply(scantwo.threshold, function(x) unname(x[1, 1]))
        
        # Output scantwo-derived QTL model penalties.
        penalties <- data.frame( matrix( c(
            thresholds[6],                # main effect penalty
            thresholds[3],                # heavy interaction penalty
            thresholds[2] - thresholds[6] # light interaction penalty
        ), nrow=1, dimnames=list(NULL, c('main', 'heavy', 'light') ) ) )
        
        writeResultHDF5(penalties, tmp, phenotypes[i], analysis, 'Penalties')
        
        # Get significant 2-QTL pairs.
        qtl.pairs <- summary(pheno.result, thresholds=thresholds[1:5])
        num.pairs <- nrow(qtl.pairs)
        
        # Output any significant QTL pairs.
        if ( num.pairs > 0 ) {
            writeResultHDF5(qtl.pairs, tmp, phenotypes[i], analysis, 'QTL Pairs')
        }
        
        # Set scantwo results overview info for this phenotype.
        suffix <- if ( num.pairs == 1 ) { 'QTL Pair' } else { 'QTL Pairs' }
        scantwo.overview[i] <- paste(num.pairs, suffix)
        
        # If updating HDF5 file, add scantwo results to updated objects.
        if (updating.h5file) {
            h5name <- joinH5ObjectNameParts( c('Results', phenotypes[i], analysis) )
            updated.objects <- c(updated.objects, h5name)
        }
    } 
    
    # If updating existing HDF5 file..
    if (updating.h5file) {
      
        # ..output updated results overview..
        stopifnot( hasResultsOverviewHDF5(h5file) )
        overview <- readResultsOverviewHDF5(h5file)
        overview <- updateResultsOverview(overview, analysis, scantwo.overview)
        writeResultsOverviewHDF5(overview, tmp)
      
        # ..and add results overview to updated objects..
        h5name <- joinH5ObjectNameParts( c('Results', 'Overview') )
        updated.objects <- c(updated.objects, h5name)
      
    } else { # ..otherwise output new results overview.
        
        overview <- makeResultsOverview(phenotypes, analysis, scantwo.overview)
        writeResultsOverviewHDF5(overview, tmp)
    }
    
    # If updating HDF5 scan file, transfer all existing but
    # unchanged objects from HDF5 scan file to temp file.
    if (updating.h5file) {
        updated <- unique( unlist( lapply( updated.objects,
            function(h5name) getObjectNamesHDF5(tmp, h5name) ) ) )
        existing <- getObjectNamesHDF5(h5file)
        unchanged <- setdiff(existing, updated) # NB: guarantees unique
        copyObjectsHDF5(h5file, tmp, h5names=unchanged)
    }
    
    # Move temp file to final scan result file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, h5file, overwrite=TRUE)
      
    return( invisible() )
}

# End of run_scantwo.R #########################################################