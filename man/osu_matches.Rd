% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rexmap.R
\name{osu_matches}
\alias{osu_matches}
\title{Generate a table with a full list of all best-matched strains}
\usage{
osu_matches(blast_output, osu_offset = rexmap_option("osu_offset"))
}
\arguments{
\item{blast_output}{Blast output class. Output from \code{\link{blast}}
function.}

\item{osu_abundances}{Data table with OSU abundances; output from the
\code{\link{abundance}} function.}
}
\description{
Generate a table with a full list of all best-matched strains
}
