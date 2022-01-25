# Read FASTA and FASTQ files regardless of number of characters in a row
# (a limitation of ShortRead::readFastq()).

#' Read file into a character vector
#'
#' @export
file_reader = function (file_name) {
  # read entire file at once into a list. used for reading our FASTA database
  # Basically same runtime as C++ code.
  size = file.info(file_name)$size
  strsplit(readChar(file_name, size, useBytes=T), "\n", fixed=T, useBytes=T)[[1]]
}

fasta_reader = function (file_name) {
  x = file_reader(file_name)
  list('meta'=sub('>', '', x[c(T,F)], fixed=T, useBytes=T), 'seqs'=x[c(F,T)])
}

fastq_reader = function (file_name) {
  x = file_reader(file_name)
  list('meta'=x[c(T,F,F,F)], 'seqs'=x[c(F,T,F,F)], 'qual'=x[c(F,F,F,T)])
}

sfastq_reader = function (file_name) {
  x = ShortRead::readFastq(file_name)
  xx = ShortRead::sread(x)
  list('meta'= as.character(ShortRead::id(x)),
       'seqs'=as.character(xx),
       'qual'=as.character(Biostrings::quality(Biostrings::quality(x))))
}

fastq_list_writer = function (fastq_list, output, ncpu=rexmap_option('ncpu')) {
  # Write FASTQ file at once. Input is a list where each element is a list
  # with names 'meta', 'seq' and 'qual'. I use this to save the result from (mc)mapply.
  list_paste = function (l) {
    paste(paste0('@', l$meta), l$seq, '+', paste0(l$qual, '\n'), sep='\n')
  }
  writeChar(paste(parallel::mclapply(fastq_list, list_paste, mc.cores=ncpu),
                  sep='\n'), output, eos=NULL)
}

fasta_writer = function (meta, seqs, output) {
  # Write FASTA sequences from meta data and sequences
  writeChar(paste(paste0('>', meta), seqs, sep='\n', collapse='\n'),
            output, eos=NULL)
}

#' Saves sequences from RExMap sequence abundance table to FASTA file
#'
#' @param abundance_table Output from \code{\link{sequence_abundance}} function.
#' @param remove_from_table If TRUE, column with sequences (names sequences)
#' is removed from the data table after the FASTA file is written to disk.
#'
#' @export
sequences_to_fasta = function (abundance_table, fasta_out,
                               remove_from_table=F) {
  if ('sequence' %in% names(abundance_table)) {
    with(unique(abundance_table[, .(qseqid, sequence)])[order(qseqid)],
      fasta_writer(
        qseqid,
        sequence,
        fasta_out
      )
    )
    if (remove_from_table) abundance_table[, sequence := NULL]
  } else {
    warning('No sequences column in the table. Nothing to do.')
  }
}
