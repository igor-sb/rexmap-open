% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{read_files}
\alias{read_files}
\title{Read files from a folder satisfying a pattern.}
\usage{
read_files(path, pattern = "")
}
\arguments{
\item{path}{Full path to the folder to list files in.}

\item{pattern}{Pattern for pattern matching. If '' then
list all files.}
}
\description{
Shortcut for \code{sort(dir(path, pattern, full.names=T))}
}
\examples{
read_files('~/data/diabimmune/fastq_tutorial', 'R1') # Get forward reads
read_files('~/data/diabimmune/fastq_tutorial', 'R2') # Get reverse reads


}
