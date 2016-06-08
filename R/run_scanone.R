# Start of run_scanone.R #######################################################

# run_scanone ------------------------------------------------------------------
#' Do single QTL scan.
#' 
#' @description This script will read cross data from the specified input file, 
#' run a single QTL analysis using \pkg{R/qtl} \code{scanone}, and write the
#' results to the specified output file.
#' 
#' @param infile input cross file
#' @param outfile output result file
#' @param chr sequences [default: all]
#' @param pheno phenotypes [default: all]
#' @param model phenotype model
#' @param method method of QTL analysis
#' @param n.perm number of permutations
#' @param n.cluster number of threads
#' @param alpha significance level for LOD threshold
#' @param step step size for genotype probabilities
#' @param error.prob genotyping error rate 
#' @param map.function genetic map function
#' 
#' @export
#' @rdname run_scanone
run_scanone <- function(infile, outfile, chr=NA, pheno=NA, model=c('normal',
    'binary', '2part', 'np'), method=c('em', 'imp', 'hk', 'ehk', 'mr', 'mr-imp',
    'mr-argmax'), n.perm=1000L, n.cluster=1L, alpha=0.05, step=0,
    error.prob=0.0001, map.function=c('haldane', 'kosambi', 'c-f', 'morgan')) {
    
    stopifnot( isSingleString(infile) )
    stopifnot( file.exists(infile) )
    stopifnot( isSingleString(outfile) )
    stopifnot( isSingleProbability(alpha) )
    stopifnot( isSingleNonNegativeNumber(step) )
    stopifnot( isSingleProbability(error.prob) )
    
    model <- match.arg(model)
    method <- match.arg(method)
    map.function <- match.arg(map.function)
    
    chr <- if ( ! is.na(chr) ) { splitCSL(chr) } else { NULL }
    pheno <- if ( ! is.na(pheno) ) { splitCSL(pheno) } else { NULL }

    # Get parameter alpha as a percentage.
    pct.alpha <- paste0(alpha * 100, '%')
    
    # Read cross input file.
    cross <- readCrossCSV(infile, error.prob=error.prob,
        map.function=map.function)
    
    # Get cross info.
    cross.info <- attr(cross, 'info')
    sequences <- getSequences(cross.info, chr)
    phenotypes <- getPhenotypes(cross.info, pheno)
    pheno.names <- getPhenotypeNames(cross.info, pheno)
    pheno.col <- getPhenoColIndices(cross, pheno.names)
    
    # Get cross map.
    cross.map <- qtl::pull.map(cross)
    
    # Calculate probabilities of true underlying genotypes given observed marker data.
    cross <- qtl::calc.genoprob(cross, step=step, error.prob=error.prob, 
        map.function=map.function)
    
    # Run single QTL analysis using R/qtl.
    scanone.result <- batchPhenoScanone(cross, chr=sequences, 
        pheno.col=pheno.col, model=model, method=method, n.cluster=n.cluster)
    
    # Run permutation scans.
    scanone.perms <- batchPermScanone(cross, chr=sequences, pheno.col=pheno.col, 
        model=model, method=method, n.perm=n.perm, n.cluster=n.cluster)
    
    # Get LOD thresholds from permutation results. 
    thresholds <- qtl:::summary.scanoneperm(scanone.perms, alpha=alpha)
    
    if ( file.exists(outfile) ) {
        file.remove(outfile)
    }
    
    # Write map to output file.
    writeMapHDF5(cross.map, outfile)
    
    comments <- character( length(phenotypes) )
    status <- logical( length(phenotypes) )
    
    # Output results of single QTL analysis for each phenotype.
    for ( i in getIndices(phenotypes) ) {
        
        # Output scan result for this phenotype.
        pheno.result <- qtl:::subset.scanone(scanone.result, lodcolumn=i)
        writeResultHDF5(pheno.result, outfile, phenotypes[i])
        
        # Output permutation scan results for this phenotype.
        pheno.perms <- qtl:::subset.scanoneperm(scanone.perms, lodcolumn=i)
        writeResultHDF5(pheno.perms, outfile, phenotypes[i])
        
        # Get significant QTL intervals.
        qtl.intervals <- getQTLIntervals(scanone.result, lodcolumn=i, 
            threshold=thresholds[pct.alpha, i])
        
        # Output any significant QTL intervals.
        if ( ! is.null(qtl.intervals) ) {
            writeResultHDF5(qtl.intervals, outfile, phenotypes[i])
            comments[i] <- paste(length(qtl.intervals), 'QTLs')
            status[i] <- TRUE
        }
    }
    
    # Output results overview.
    overview <- data.frame(Phenotype=phenotypes, Status=status, 
        Comments=comments, stringsAsFactors=FALSE)
    writeOverviewHDF5(overview, outfile)
    
    return( invisible() )
}

# End of run_scanone.R #########################################################