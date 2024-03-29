% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastx_io.R
\name{sequences_to_fasta}
\alias{sequences_to_fasta}
\title{Saves sequences from RExMap sequence abundance table to FASTA file}
\usage{
sequences_to_fasta(abundance_table, fasta_out, remove_from_table = F)
}
\arguments{
\item{abundance_table}{Output from \code{\link{sequence_abundance}} function.}

\item{remove_from_table}{If TRUE, column with sequences (names sequences)
is removed from the data table after the FASTA file is written to disk.}
}
\description{
Saves sequences from RExMap sequence abundance table to FASTA file
}
