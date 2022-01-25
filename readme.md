# RExMap: Reference-based exact Mapping Pipeline

A 16S pipeline for a reference-based local alignment of 16S sequences to the reference database RExMapDB of known isolate microbial strains.

## Installation from source

### Linux

For Linux, we need to make sure you have development version of OpenSSL installed, so e.g. for Ubuntu run from terminal:
```
sudo apt-get -y install libssl-dev libcurl4-openssl-dev
```
Then to install the latest version of R, we need to add the repository, authentication key and install it:
```
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
sudo apt-get update
sudo apt-get -y install r-base
```

In R, install devtools package, then use `install_github` to install the package from source:

```R
install.packages('devtools')
devtools::install_github('benjjneb/dada2')
# Select 1 when prompted to install biocLite
devtools::install_github('sapientbio/rexmap-open')
```

This may take up to about 10 minutes to install on a fresh R installation.

### macOS

First, make sure you have the latest Xcode and [latest R](https://cran.r-project.org/bin/macosx/R-3.5.0.pkg) (> 3.4) installed. Then, install Xcode command line tools. To do this, open a macOS Terminal, then type
```
xcode-select --install
```

Then open either R console or RStudio and run this:
```R
# Install devtools first
options('install.packages.compile.from.source'='never')
install.packages('devtools')

# DADA2 pre-requisites
source('http://bioconductor.org/biocLite.R')
biocLite('GenomeInfoDbData', suppressUpdates=T)


# DADA2
devtools::install_github('benjjneb/dada2')

# Finally install RExMap
devtools::install_github('sapientbio/rexmap-open')
```


## Tutorials

* [Human gut microbiome data tutorial](tutorial.ipynb) (main tutorial with all explanations)
* [mock community tutorial](tutorial_mock.ipynb)

## Notes

By default, multithreading is enabled and supported only on Linux and macOS, through R package `parallel`.

If DADA2 package installation fails due to `latticeExtra` package requiring R version 3.6 and you are running 3.5, you can use this to install an older version of the `latticeExtra` package:

```R
devtools::install_version('latticeExtra', version='0.6-28', repos='http://cran.us.r-project.org')
```

RExMap is provided under the Attribution-NonCommercial-ShareAlike 4.0 International license.
