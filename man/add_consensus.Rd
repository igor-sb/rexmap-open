% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rexmap.R
\name{add_consensus}
\alias{add_consensus}
\title{Untrim trimmed sequences in a DADA2 object}
\usage{
add_consensus(
  dada_res,
  derep,
  fq_tri,
  fq_fil,
  truncLen,
  verbose = T,
  ncpu = rexmap_option("ncpu")
)
}
\arguments{
\item{dada_res}{DADA2 object with sequences to un-trim.}

\item{derep}{DADA2 derep object that was used to obtain dada_res.}

\item{fq_tri}{A character vector of FASTQ filenames pre-global-trimming
(but after PCR primer trim).}

\item{fq_fil}{A character vector of FASTQ filenames post-global-trimming and filtering.}

\item{truncLen}{Length used to trim.}

\item{verbose}{Boolean specifying whether to display progress bar.}

\item{ncpu}{Integer specifying number of CPU threads to use. This uses R package "parallel"
so works only on macOS and Linux.}
}
\description{
Untrim trimmed sequences in a DADA2 object
}
