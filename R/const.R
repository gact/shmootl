# Start of const.R #############################################################

# Package Constants ------------------------------------------------------------
#' Package constants for \pkg{shmootl}.
#'  
#' @docType package
#' @include util.R
#' @keywords internal
#' @name Package Constants
NULL

# Create environment for package constants -------------------------------------
# NB: only change these settings if you know exactly what you're doing.
const <- new.env()
with(const, {
    
    # Dataframe constructor keyword arguments
    dataframe.keywords <- c('row.names', 'check.rows', 'check.names', 
        'fix.empty.names', 'stringsAsFactors')
    
    # Special R object attributes.
    special.attributes <- list(
        default = c('class', 'comment', 'dim', 'dimnames', 'names', 
            'row.names', 'tsp'),
        DNAStringSet = c('pool', 'ranges', 'elementType', 
            'elementMetadata', 'metadata', 'class'),
        QualityScaledDNAStringSet = c('pool', 'ranges', 'elementType', 
            'elementMetadata', 'metadata', 'quality', 'class')
    )
    
    # Phenotypes not used in yeast QTL analysis.
    disallowed.phenotypes <- c('pgm', 'sex')
    
    # Reserved phenotypes.
    reserved.phenotypes <- c(disallowed.phenotypes, 'id')
    
    # Tetrad sample labels (for alphanumeric tetrad sample IDs).
    tetrad.sample.labels <- c('A', 'B', 'C', 'D')
    
    # Quality score info.
    qual = list(
        phred = list(
            range = c(0, 93),
            offset = 33
        )
    )
    
    # Item IDs -----------------------------------------------------------------
    reserved.characters <- c(
        "'", # single quote: quote character
        '"', # double quotes: quote character
        '`', # backtick: quote character
        ',', # comma: delimiter in CSV files
        '/', # forward slash: delimiter in HDF5 object names
        '\\' # backslash: common escape character
    )
    
    printable.ascii <- unlist( strsplit(rawToChar( as.raw(32:126) ), '') )
    
    id.charset <- printable.ascii[ ! printable.ascii %in% reserved.characters ]
    id.charset <- c( ']', id.charset[ id.charset != ']' ] )
    
    # Defaults -----------------------------------------------------------------
    
    default = list(
        genome = 'S288C_R64-1-1'
    )
    
    # Regex patterns -----------------------------------------------------------
    pattern <- list(
        
        # Pattern for simplification of chromosome names. This is based on the 
        # method used by SnpEff (Cingolani et al. 2012 [PMID:22728672]) to 
        # simplify chromosome names, as used in 'ChromosomeSimpleName.java'.
        # Prefixes of form 'c', 'chr', 'chromo', and 'chromosome' are removed. 
        # For more info, see: <https://github.com/pcingola/SnpEff>
        # [Accessed: 16 Feb 2016].
        chromosome = '^(?:c(?:hr(?:omo(?:some)?)?)?)?([^[:space:]]+)$',
        
        # Default marker ID.
        default.marker.id = '^c([[:digit:]]{2})[-:]([[:digit:]]{7})$',
        
        # Default QTL name.
        default.qtl.name = '^([[:alnum:]]+)@(-?[[:digit:]]+(?:[.][[:digit:]]+)?)$',
        
        # File name.
        file.with.ext = '^(.+)[.]([^.]+)$',
        
        # Genetic map units.
        gmap = c(cM='\\s*cM$'),
        
        # Haploid VCF genotype.
        vcf.haploid.geno = '^([.]|[[:digit:]]+)$',
        
        # HDF5 default element name.
        h5element = '^ELT([[:digit:]]+)$',
        
        # Valid ID.
        item.id = paste0('^[', paste(id.charset,  collapse=''), ']+$'),
        
        # Numeric name.
        numeric.name = 'X([.])?([[:digit:]]+(?:[.][[:digit:]]+)?)',
        
        # Package pipeline docs.
        pipe.docs = '^(run_([[:alpha:]][[:alpha:]._]*))[.]Rd$',
        
        # Package pipeline functions.
        pipe.func = '^run_([[:alpha:]][[:alpha:]._]*)$',
        
        # Mapframe position column heading.
        poscol = '(?:^|^.*[^[:alpha:]])pos(?:[^[:alpha:]].*$|$)',
        
        # Physical map units.
        pmap = c(Mb='\\s*Mbp?$', kb='\\s*kbp?$', bp='\\s*bp$'),
        
        # Pseudomarker ID.
        pseudomarker.id = '^c([[:alnum:]]+)[.]loc(-?[[:digit:]]+(?:[.][[:digit:]]+)?)$',
        
        # Tetrad sample ID.
        tetrad = c( 
            alphanumeric = paste0('.*?([[:digit:]]+)([', 
                                  paste(tetrad.sample.labels, collapse=''), '])'), 
            numeric = '^.+?([[:digit:]]+)$'
        )
    )
    
    # Chromosome/sequence info -------------------------------------------------
    
    # Load chromosome info.
    chrinfo <- loadChrInfo()
    
    # Load sequence info.
    seqinfo <- loadSeqInfo()
    
    # Get mapping of sequence aliases to standard names.
    alias2chrom <- unlist( lapply(1:nrow(chrinfo), function(i) {
        aliasmap <- character()
        if ( ! is.na(chrinfo$aliases[i]) ) {
            for ( alias in strsplit(chrinfo$aliases[i], paste0('[^', 
                paste(id.charset,  collapse=''), ']'))[[1]] ) {
                aliasmap[[ toupper(alias) ]] <- chrinfo$seqnames[i]
            }
        }
        return(aliasmap)
    }) )
    
    # Map info -----------------------------------------------------------------
    map.info <- list()
    
    # Genetic map info.
    known.gmap.units <- 'cM'
    map.info[['gmap']] <- data.frame(pattern=pattern$gmap, 
        factor=c(cM=1), stringsAsFactors=FALSE)
    
    # Physical map info.
    known.pmap.units <- c('Mb', 'kb', 'bp')
    map.info[['pmap']] <- data.frame(pattern=pattern$pmap, 
        factor=c(Mb=1000000, kb=1000, bp=1), stringsAsFactors=FALSE)
    
    # Map units.
    known.map.units <- c(known.gmap.units, known.pmap.units)
    basic.map.unit <- c(gmap='cM', pmap='bp')
    
    # Map types (named by map unit).
    known.map.types <- c(rep_len('gmap', length(known.gmap.units)), 
        rep_len('pmap', length(known.pmap.units)))
    names(known.map.types) <- known.map.units
    
    # Minimum number of sequences per map.
    min.spm = 1
    
    # Minimum number of loci per sequence.
    min.lps = 2
    
    # CrossInfo object sample info ---------------------------------------------
    
    sample.levels <- c('sample', 'strain', 'tetrad')
    sample.aspects <- data.frame(
        id=c('sample.id', 'sample.id', NA), 
        name=c('sample.name', 'sample.name', NA),
        index=c('sample.index', 'strain.index', 'tetrad.index'), 
        row.names=sample.levels, 
        stringsAsFactors=FALSE
    )
    sample.headings <- unique( na.omit( c(sample.aspects$id, 
        sample.aspects$name, sample.aspects$index) ) )
    marker.headings <- c('marker', 'seq')
    
    # Default QTL analysis result names ----------------------------------------
    default.result.names <- data.frame( matrix(c(
        #           class                    name
        'scanone',              'Scanone',
        'scantwo',              'Scantwo',
        'scanoneperm', 'Scanone Permutations',
        'scantwoperm', 'Scantwo Permutations',
        'qtlintervals',        'QTL Intervals'
    ), ncol=2, byrow=TRUE, dimnames=list(NULL, c('class', 'name'))), 
    stringsAsFactors=FALSE)
    
    # Scan function arguments --------------------------------------------------
    
    # Set documented scan function arguments.
    scan.args <- list(
        'funqtl::scanoneF' = c('cross', 'pheno.cols', 'n.perm'),
        'funqtl::scantwoF' = c('cross', 'pheno.cols', 'usec', 'n.perm'),
        'qtl::scanone' = c('cross', 'chr', 'pheno.col', 'model', 'method', 
            'addcovar', 'intcovar', 'weights', 'use', 'upper', 'ties.random', 
            'start', 'maxit', 'tol', 'n.perm', 'perm.Xsp', 'perm.strata', 
            'verbose', 'batchsize', 'n.cluster', 'ind.noqtl'),
        'qtl::scantwo' = c('cross', 'chr', 'pheno.col', 'model', 'method', 
            'addcovar', 'intcovar', 'weights', 'use', 'incl.markers', 
            'clean.output', 'clean.nmar', 'clean.distance', 'maxit', 'tol', 
            'verbose', 'n.perm', 'perm.Xsp', 'perm.strata', 'assumeCondIndep', 
            'batchsize', 'n.cluster')
    )
    
    # Scan result attributes.
    scan.attributes <- list(
        'scanone' = c('method', 'type', 'model'),
        'scantwo' = c('method', 'type', 'fullmap', 'phenotypes'),
        'scantwoperm' = c('method')
    )
    
    # R/argparser settings -----------------------------------------------------
    
    # Special arguments.
    special.params <- c('help', 'opts')
    
    # Parameter groups.
    param.groups <- c('positional', 'optional', 'flags')
})

# Lock constants ---------------------------------------------------------------
lockEnvironment(const, bindings=TRUE)

# End of const.R ###############################################################