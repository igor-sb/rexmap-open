% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rexmap.R
\name{abundance}
\alias{abundance}
\title{OSU abundance table}
\usage{
abundance(
  abundance_table,
  blast_object,
  ncpu = rexmap_option("ncpu"),
  verbose = rexmap_option("verbose"),
  raw_strains = TRUE,
  pso_n = 1000,
  custom_sampleids = NULL,
  debug = FALSE
)
}
\arguments{
\item{abundance_table}{Sequence abundance table. Data table with 3 columns
(in this exact order): sample_id, qseqid, raw_count.}

\item{blast_object}{BLAST output class, output from \code{\link{blast}} function.}

\item{ncpu}{Integer specifying number of CPU threads to use. This uses R package "parallel"
(TRUE) or to simplify the output when multiple strains of the same species are in the
same OSU (FALSE).}

\item{verbose}{TRUE/FALSE: display status messages}

\item{raw_strains}{TRUE/FALSE: whether to report full strain information for each OSU.
If this is true "species" column will display a full list of best matching strains
for each OSU. Since this last may be very long, use \code{\link{print_strains}}
to make it more readable.}

\item{pso_n}{Integer specifying number of times to run Particle Swarm
Optimizer for an OSU abundance estimation algorithm. (Default: 1000)}

\item{custom_sampleids}{process only specified sampleids (separated by comma)
If NULL (default) process all sample ids.}
}
\description{
Combine BLAST object and a sequence abundance table into a final OSU
abundance table.
}
