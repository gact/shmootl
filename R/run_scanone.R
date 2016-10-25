# Start of run_scanone.R #######################################################

# run_scanone ------------------------------------------------------------------
#' Do single QTL scan.
#' 
#' Read cross data from the specified cross input file, run a single QTL
#' analysis using \pkg{R/qtl} \code{scanone} (Broman \emph{et al.} 2003),
#' and write the results of that scan to the specified output file.
#' 
#' If the input cross contains enumerated genotypes, marker regression is
#' performed regardless of the value of the \code{method} parameter.
#' 
#' In typical usage, LOD threshold stringency can be set through either the
#' significance level (\code{alpha}), or the false-discovery rate (\code{fdr}),
#' but not both. If neither is specified, a significance level \code{alpha}
#' of \code{0.05} is used. The given stringency is then used to estimate a
#' LOD threshold from \code{scanone} permutations.
#' 
#' This can be bypassed by setting a fixed LOD \code{threshold}, along with a
#' nominal stringency (\code{alpha} or \code{fdr}), in which case permutations
#' are skipped and the fixed LOD threshold is applied directly for assessing
#' significance.
#' 
#' LOD interval estimation can be controlled with the \code{'ci.function'}
#' parameter: set to \code{'lodint'} for LOD support intervals (adjusting
#' stringency with the \code{'drop'} parameter), or to \code{'bayesint'}
#' for Bayesian credible intervals (adjusting stringency with the \code{'prob'}
#' parameter). For more information on the QTL interval methods used, see
#' functions \code{'lodint'} and \code{'bayesint'} in the \pkg{R/qtl} manual,
#' as well as Section 4.5 of Broman and Sen (2009).
#' 
#' @param infile input cross file
#' @param h5file scan result file
#' @param chr sequences [default: all]
#' @param pheno phenotypes [default: all]
#' @param model phenotype model
#' @param method method of QTL analysis
#' @param n.perm number of permutations
#' @param n.cluster number of threads
#' @param alpha significance level for LOD threshold
#' @param fdr FDR for LOD threshold
#' @param threshold fixed LOD threshold
#' @param step step size for genotype probabilities
#' @param map.function genetic map function
#' @param error.prob genotyping error rate
#' @param ci.function QTL interval function
#' @param drop LOD support interval drop
#' @param prob Bayesian credible interval probability
#' @param acovfile additive covariates file
#' @param icovfile interactive covariates file
#' 
#' @template author-thomas-walsh
#' @template author-yue-hu
#' @template ref-broman-2003
#' @template ref-broman-2009
#' @template seealso-rqtl-manual
#' 
#' @concept shmootl:pipelines
#' @export
#' @family pipeline functions
#' @rdname run_scanone
run_scanone <- function(infile=NA_character_, h5file=NA_character_,
    chr=character(), pheno=character(), model=c('normal','binary','2part','np'),
    method=c('em','imp','hk','ehk','mr','mr-imp','mr-argmax'), n.perm=1000L,
    n.cluster=1L, alpha=NA_real_, fdr=NA_real_, threshold=NA_real_, step=0,
    map.function=c('haldane','kosambi','c-f','morgan'), error.prob=0.0001,
    ci.function=c('lodint', 'bayesint'), drop=1.5, prob=0.95,
    acovfile=NA_character_, icovfile=NA_character_) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    stopifnot( isSingleString(h5file) )
    stopifnot( isSingleNonNegativeNumber(step) )
    stopifnot( isSingleProbability(error.prob) )
    
    model <- match.arg(model)
    method <- match.arg(method)
    map.function <- match.arg(map.function)
    
    chr <- if ( ! identical(chr, character()) ) { chr } else { NULL }
    pheno <- if ( ! identical(pheno, character()) ) { pheno } else { NULL }
    alpha <- if ( ! identical(alpha, NA_real_) ) { alpha } else { NULL }
    fdr <- if ( ! identical(fdr, NA_real_) ) { fdr } else { NULL }
    threshold <- if ( ! identical(threshold, NA_real_) ) { threshold } else { NULL }
    
    if ( ! is.null(alpha) && ! is.null(fdr) ) {
        stop("cannot set both significance level (alpha) and FDR")
    } else if ( ! is.null(alpha) ) { # set specified alpha value
        stopifnot( isSingleProbability(alpha) )
        perm.type <- 'max'
    } else if ( ! is.null(fdr) ) { # set specified FDR
        stopifnot( isSingleFiniteNumber(fdr) )
        stopifnot( fdr > 0 & fdr < 1 )
        perm.type <- 'bins'
    } else if ( is.null(threshold) ) { # set default alpha value
        perm.type <- 'max'
        alpha <- 0.05
    } else {
        stop("cannot set fixed LOD threshold without nominal significance level (alpha) or FDR")
    }
    
    # Attach required package namespaces (if needed) ---------------------------
    
    req.pkgs <- 'qtl'
    names(req.pkgs) <- paste0('package:', req.pkgs)
    att.pkgs <- req.pkgs[ ! names(req.pkgs) %in% search() ]
    sapply(att.pkgs, attachNamespace)
    on.exit( sapply(names(att.pkgs), detach, character.only=TRUE), add=TRUE )
    
    # --------------------------------------------------------------------------
    
    # Read cross input file.
    cross <- readCrossCSV(infile, error.prob=error.prob,
        map.function=map.function)
    
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
    
    # Get cross map.
    cross.map <- qtl::pull.map(cross)
    
    # Calculate probabilities of true underlying genotypes given observed marker data.
    cross <- qtl::calc.genoprob(cross, step=step, error.prob=error.prob, 
        map.function=map.function)
    
    # Run single QTL analysis using R/qtl.
    scanone.result <- batchPhenoScanone(cross, chr=sequences, 
        pheno.col=pheno.col, model=model, method=method, n.cluster=n.cluster,
        addcovar=addcovar, intcovar=intcovar)
    
    # If no fixed threshold specified, estimate from permutations..
    if ( is.null(threshold) ) {
        
        # Run permutation scans.
        scanone.perms <- batchPermScanone(cross, chr=sequences,
            pheno.col=pheno.col, model=model, method=method, n.perm=n.perm,
            n.cluster=n.cluster, perm.type=perm.type, addcovar=addcovar,
            intcovar=intcovar)
        
        # Get LOD thresholds from permutation results.
        if ( ! is.null(alpha) ) {
            scanone.thresholds <- summary(scanone.perms, alpha=alpha) # qtl:::summary.scanoneperm
            multiple.lodcolumns <- TRUE
        } else { # fdr
            scanone.threshold <- summary(scanone.perms, scanone.result, fdr=fdr) # summary.scanonebins
            # NB: set FDR from scanonebins summary, as can differ from requested FDR
            fdr <- 0.01 * as.numeric(sub('%', '', rownames(scanone.threshold)[1]))
            multiple.lodcolumns <- FALSE
        }
        
    } else { # ..otherwise make LOD threshold object.
        
        scanone.threshold <- makeScanoneThresholdObject(threshold,
            alpha=alpha, fdr=fdr)
        multiple.lodcolumns <- FALSE
        scanone.perms <- NULL
    }
    
    # Create temp output file, ensure will be removed.
    tmp <- tempfile()
    on.exit( file.remove(tmp), add=TRUE )
    
    # Write map to temp output file.
    writeMapHDF5(cross.map, tmp, 'Genetic Map')
    
    comments <- character( length(phenotypes) )
    status <- logical( length(phenotypes) )
    
    # Output results of single QTL analysis for each phenotype.
    for ( i in seq_along(phenotypes) ) {
        
        # Output scanone result for this phenotype.
        pheno.result <- getLODProfile(scanone.result, lodcolumn=i)
        writeResultHDF5(pheno.result, tmp, phenotypes[i], 'Scanone')
        
        # Output any permutation scan results for this phenotype.
        if ( ! is.null(scanone.perms) ) {
            
            if ( ! is.null(alpha) ) {
                pheno.perms <- subset(scanone.perms, lodcolumn=i) # qtl:::subset.scanoneperm
            } else { # fdr
                pheno.perms <- scanone.perms[,, i]
            }
            writeResultHDF5(pheno.perms, tmp, phenotypes[i], 'Scanone Perms')
        }
        
        # Output scanone threshold for this phenotype.
        if (multiple.lodcolumns) {
            scanone.threshold <- subset(scanone.thresholds, lodcolumn=i)
        }
        writeResultHDF5(scanone.threshold, tmp, phenotypes[i], 'Scanone Threshold')
        
        # Get significant QTL intervals.
        qtl.intervals <- getQTLIntervals(pheno.result, ci.function=ci.function,
            drop=drop, prob=prob, threshold=scanone.threshold)
        
        # Output any significant QTL intervals.
        if ( length(qtl.intervals) > 0 ) {
            writeResultHDF5(qtl.intervals, tmp, phenotypes[i], 'Scanone QTL Intervals')
            comments[i] <- paste(length(qtl.intervals), 'QTLs')
            status[i] <- TRUE
        }
    }
    
    # Output results overview.
    overview <- data.frame(Phenotype=phenotypes, Status=status,
        Comments=comments, stringsAsFactors=FALSE)
    writeOverviewHDF5(overview, tmp)
    
    # Move temp file to final scan result file.
    # NB: file.copy is used here instead of file.rename because the latter
    # can sometimes fail when moving files between different file systems.
    file.copy(tmp, h5file, overwrite=TRUE)
    
    return( invisible() )
}

# End of run_scanone.R #########################################################