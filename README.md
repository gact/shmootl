# ShmooTL: QTL analysis utilities for yeast

This package contains pipelines and utilities for QTL analysis of yeast cross
data. Pipelines can be run from the command line or from within the R
environment. Utility functions are available from within R.

To find more information on using **ShmooTL**, see the [Usage](#usage) and
[Documentation](#documentation) sections.

## Dependencies

This package depends on the following:

- [argparser](https://cran.r-project.org/web/packages/argparser/index.html)
- [Bioconductor](http://www.bioconductor.org/)
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
- [funqtl](https://github.com/ikwak2/funqtl)
- [qtl](http://www.rqtl.org/)
- [qtlcharts](http://kbroman.org/qtlcharts/)
- [rhdf5](http://bioconductor.org/packages/release/bioc/html/rhdf5.html)
- [rtracklayer](http://bioconductor.org/packages/devel/bioc/html/rtracklayer.html)
- [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
- [yaml](https://cran.r-project.org/web/packages/yaml/index.html)

In addition, the following suggested packages are each needed for a specific use:

- [abind](https://cran.r-project.org/web/packages/abind/index.html): estimating FDR in single-QTL analysis
- [viridis](https://cran.r-project.org/web/packages/viridis/index.html): plotting with the Viridis palette
- [xlsx](https://cran.r-project.org/web/packages/xlsx/index.html): writing Excel output

## Installation

Before installing **ShmooTL** from GitHub, you will first need to install its
dependency packages.

### Installing packages from R

Some required packages can be installed from within R as follows:

```
install.packages('argparser')
install.packages('qtl')
install.packages('yaml')
```

Suggested packages `abind`, `viridis` and `xlsx` can also be installed in this way:

```
install.packages('abind')
install.packages('viridis')
install.packages('xlsx')
```

### Installing Bioconductor packages

To install required Bioconductor packages (e.g. `VariantAnnotation`, `rhdf5`), see the [Bioconductor installation documentation page](https://www.bioconductor.org/install/).

### Installing packages from GitHub

Packages `funqtl` and `qtlcharts` can be installed directly from the GitHub
website, but you must first install the `devtools` package, which can be
installed from within R:

```
install.packages('devtools')
```

To install `funqtl` from GitHub, input:

```
library(devtools)
install_github('ikwak2/funqtl')
```

To install `qtlcharts` from GitHub, install the package `hmtlwidgets` as
follows:

```
install.packages('htmlwidgets')
```

Then install `qtlcharts` from GitHub:

```
library(devtools)
install_github('kbroman/qtlcharts')
```

### Installing ShmooTL from GitHub

With the `devtools` package and other dependencies installed,
**ShmooTL** can be installed from GitHub as follows:

```
library(devtools)
install_github('gact/shmootl')
```

If you want to access vignettes for the installed package,
set the `build_vignettes` option to `TRUE`, as follows:

```
library(devtools)
install_github('gact/shmootl', build_vignettes=TRUE)
```

**ShmooTL** vignettes can then be accessed as described in the
[Documentation](#documentation) section.

## Usage

Package pipelines can be run from the command line using Rscript as follows:

```
Rscript -e 'library(shmootl)' -e 'run()' <pipeline> [-h] [<args>]
```

...where `run` is the general **ShmooTL** pipeline-running function,
`<pipeline>` is the name of the pipeline to run, and `<args>` represents
any arguments to be passed to the given pipeline. To see the available
options for a pipeline, input the help flag (`-h`) after the name of the
pipeline.

## Documentation

To view the **ShmooTL** help pages within R, input the following command:

```
help(package=shmootl)
```

To view the help page for a specific **ShmooTL**
function, input a command such as the following:

```
help(readCrossCSV)
```

...where the first argument (in this case `readCrossCSV`) is the name
of the function for which help is being sought.

If package vignettes were built when
[installing **ShmooTL**](#installing-shmootl-from-github), the
**ShmooTL** user guide can be viewed with the following R command:

```
vignette('shmootl-guide')
```

Passing the string `'shmootl-devnotes'` to the same function displays the **ShmooTL** developer notes.

## Contact

For further information, contact the package maintainer listed in the package DESCRIPTION file.
Please submit any issue on the [ShmooTL GitHub repository](https://github.com/gact/shmootl).
