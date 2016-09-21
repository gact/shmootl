# ShmooTL: QTL analysis utilities for yeast

This package contains pipelines and utilities for QTL analysis of yeast cross
data. Pipelines can be run from the command line or from within the R
environment. Utility functions are available from within R.

## Dependencies 

This package depends on the following:

- [abind](https://cran.r-project.org/web/packages/abind/index.html)
- [argparser](https://cran.r-project.org/web/packages/argparser/index.html)
- [Bioconductor](http://www.bioconductor.org/)
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
- [funqtl](https://github.com/ikwak2/funqtl)
- [qtl](http://www.rqtl.org/)
- [qtlcharts](http://kbroman.org/qtlcharts/)
- [rhdf5](http://bioconductor.org/packages/release/bioc/html/rhdf5.html)
- [rtracklayer](http://bioconductor.org/packages/devel/bioc/html/rtracklayer.html)
- [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
- [xlsx](https://cran.r-project.org/web/packages/xlsx/index.html)
- [yaml](https://cran.r-project.org/web/packages/yaml/index.html)

## Installation 

To install ShmooTL from GitHub, you will first need to install [devtools](https://github.com/hadley/devtools) using the R command:

```
install.packages('devtools')
```

The ShmooTL package and its dependencies can then be installed with the R commands:

```
library(devtools)
install_github('gact/shmootl')
```

## Usage 

Package pipelines can be run from the command line using Rscript as follows:

```
Rscript -e 'library(shmootl)' -e 'run()' <pipeline> [-h] [<args>]
```

...where `run` is the general ShmooTL pipeline-running function, `<pipeline>`
is the name of the pipeline to run, and `<args>` represents any arguments to be
passed to the given pipeline. To see the available options for a pipeline, input
the help flag (`-h`) after the name of the pipeline.

## Contact

For information or issues email tw164 (-a-) le.ac.uk

