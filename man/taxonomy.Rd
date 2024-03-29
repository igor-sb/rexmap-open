% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{taxonomy}
\alias{taxonomy}
\title{Add NCBI taxonomy to an OSU table}
\usage{
taxonomy(
  osu_table,
  verbose = rexmap_option("verbose"),
  show_count = TRUE,
  ws = "_",
  split_char = ","
)
}
\arguments{
\item{table}{Data table containing columns 'osu_id', 'pctsim' and 'species'.
Typically used with the output table from then \code{\link{abundance}}
function.}
}
\description{
Add NCBI taxonomy to an OSU table
}
