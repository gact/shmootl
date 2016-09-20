#' @section Raw Alleles and Genotypes:
#' 
#' A raw allele is the original marker value. Single-nucleotide variants (i.e.
#' \code{'A'}, \code{'C'}, \code{'G'}, or \code{'T'}) are commonly used, and
#' are the only kind of marker supported by the \pkg{shmootl} package. These
#' must be converted into another kind of allele before doing QTL analysis.
#' 
#' For a given ploidy \code{N}, a raw genotype is the concatenation of the set
#' of \code{N} raw alleles. For example, given two \code{'A'} alleles at a
#' homozygous diploid marker, the genotype is \code{'AA'}.
