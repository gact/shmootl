% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_makegeno.R
\name{run_makegeno}
\alias{run_makegeno}
\title{Make genotype data from VCF.}
\usage{
run_makegeno(datafile = NA_character_, genfile = NA_character_,
  fdrfile = NA_character_, alleles = mapping(), digits = NA_integer_)
}
\arguments{
\item{datafile}{sample VCF file [required]}

\item{genfile}{output genotype CSV file [required]}

\item{fdrfile}{optional founder VCF file}

\item{alleles}{founder allele symbol mapping}

\item{digits}{numeric precision [default: unrounded]}
}
\description{
Read raw genotype data from the sample VCF file, encode the given genotype
data as founder or enumerated genotypes, and write the encoded genotypes
to an \pkg{R/qtl} genotype CSV file.
}
\details{
If no founder VCF file is specified, markers are assigned enumerated
genotypes, so that each raw allele is converted to a number in the order
in which it is first observed at a given marker (i.e. \code{'1'},
\code{'2'}, etc.).

If a founder VCF file is given, markers are assigned a genotype symbol
consisting of alleles that each correspond to a specific founder.

If the \code{alleles} parameter is specified, this must be a mapping of
founder sample IDs to allele symbols. If calling this function from within
the \code{R} environment, this must be specified as a mapping object (e.g.
\code{mapping( c(DBVPG6044 = 'W', Y12 = 'S') )}). When called from the
command line using \code{Rscript}, the \code{alleles} parameter must be
specified as a YAML string (or YAML file) mapping founders to allele
symbols (e.g. \code{"DBVPG6044: W, Y12: S"}). If the \code{alleles}
parameter is not specified, allele symbols are taken from the letters
of the alphabet (i.e. \code{'A'}, \code{'B'} etc.).
}
\seealso{
Other pipeline functions: \code{\link{run_annoqtl}},
  \code{\link{run_digest}}, \code{\link{run_estimap}},
  \code{\link{run_interptimes}},
  \code{\link{run_makecross}}, \code{\link{run_prep}},
  \code{\link{run_pullmap}}, \code{\link{run_pushmap}},
  \code{\link{run_recode}}, \code{\link{run_report}},
  \code{\link{run_scanone}}, \code{\link{run_scantwo}}
}
\concept{
shmootl:preparation
}
