# tempdir() to get temporary dir path

get_makeblastdb_path = function () {
  # Find makeblastdb location
  out = suppressWarnings(system2('which', 'makeblastdb', stdout=T))
  # If makeblastdb is not found return an error
  if (is.null(attr(out, 'status'))) return(out)
  if (attr(out, 'status') == 1) stop('Error: makeblastdb command not found.')
}

makeblastdb = function (fasta_in, db_out, verbose=TRUE,
                        makeblastdb_path=rexmap_option('path_makeblastdb')) {
  # Generate a BLAST database from fasta_in
  # db_out contains a full location and prefix for BLASTDB files
  # This function does not return anything, but can optionally
  # (if verbose == T) write the status to console.
  # makeblastdb_path = get_makeblastdb_path()
  if (verbose == TRUE) {
    output = ''
  } else {
    output = FALSE
  }
  system2(makeblastdb_path, c('-dbtype', 'nucl', '-in', fasta_in,
                              '-out', db_out), stdout=output)
}

cleanup_blastdb = function (db_out) {
  # Delete blastdb files
  db_files = paste(db_out, c('.nin', '.nhr', '.nsq'), sep='')
  for (f in db_files) {
    if (file.exists(f)) file.remove(f)
  }
}

# Normally formatted timestamp (to be used in file names)
ts = function () return(sub(' ', '_', gsub(':', '-', Sys.time())))

# Find minimal sequence length among the first N sequences from the fasta file
min_seq_len = function (fasta_files, n=50) {
   return(min(sapply(fasta_files, function (f) {
     return(nchar(fasta_reader(f)$seqs[1:n]))
   }), na.rm=T))
}

#' Collapse sequences that are exact matches, up to shifts and/or length
#'
#' Like \code{\link{dada2::collapseNoMismatch}} but a lot faster, multithreaded
#' and scales well to a very large number of unique sequences (100,000 and more).
#'
#' @importFrom igraph make_graph
#' @importFrom igraph groups
#' @importFrom igraph clusters
#' @export
collapse = function (ab_in, verbose=rexmap_option('verbose'),
                     ncpu=rexmap_option('ncpu'), temp_dir=tempdir(),
                     ws_scale=0.8, max_target_seqs=50, copy_table=TRUE) {
  # Provide ab_tab_nochim as an input argument, same as dada2::collapseNoMismatch
  if (verbose) cat('collapse:', fill=T)
  if (verbose) cat('* generating temporary files...')
  if (copy_table) {
    ab = copy(ab_in)
    if (verbose) cat('(copy OK) ')
  } else {
    ab = ab_in
  }
  out_files = ab_to_files(ab, verbose=verbose, temp_dir=temp_dir)
  fasta = out_files[1]
  db = out_files[2]
  if (verbose) cat('OK.\n')

  # Generate temp output file
  if (verbose) cat('* blast word size: ')
  blast_out = paste0(db, '_blast_output.txt')

  # Automatically generate a good BLAST word size
  ws = max(round(min_seq_len(fasta) * ws_scale), 7)
  if (verbose) cat(ws, '\n')

  # BLAST fasta_in vs db
  if (verbose) cat('* running blast...')
  blast_status = blastn(fasta, ref_db=db,
                        output=blast_out, max_target_seqs=max_target_seqs,
                        outfmt=paste0('6 ', rexmap_option('blast_coll_fmt')),
                        perc_identity=100, word_size=ws, ncpu=ncpu)
  cat('blast status: ', blast_status, fill=T)
  if (blast_status != 0) {
     stop('\nError in BLAST alignment step.')
  }
  if (verbose) cat('OK.\n')

  # Load blast output and remove the temp blast output file
  if (verbose) cat('* selecting ends-free alignments...')
  blast.dt = data.table::fread(blast_out)

  names(blast.dt) = strsplit(rexmap_option('blast_coll_fmt'), ' ', fixed=T)[[1]]
  blast2.dt = blast.dt[(qseqid != sseqid) & (pident==100)]
  # Select only ends-free alignments
  blast3.dt = blast2.dt[(qstart==1 & qend==qlen) | (sstart==1 & send==slen) |
                        (qstart==1 & send==slen) | (sstart==1 & qend==qlen)]
  if (nrow(blast3.dt) == 0) {
    # This means there aren't any sequences that need collapsing.
    # Return the input back then.
    if (verbose) cat('OK.\n')
    if (verbose) cat('* no sequences need collapsing.\n')
    if (verbose) cat('* cleaning up temporary files...')
    file.remove(blast_out)
    file.remove(fasta)
    cleanup_blastdb(db)
    if (verbose) cat('OK.\n')
    if (verbose) cat('* returning input.\n')
    return(ab_in)
  }

  # Find all connected clusters of sequence IDs
  g = make_graph(as.vector(t(blast3.dt[, .(qseqid, sseqid)])), directed=F)
  # Find all connected clusters
  cls = groups(clusters(g))
  filtered_columns = c()
  # i = 0
  for (cl in cls) {
    # i = i + 1
    # cat(i, '\n')
    if (length(cl) > 1) {
      column_sums = colSums(ab[, cl, drop=FALSE])
      # max_column is the sequence with max total number of reads
      # in the case of multiple just pick the first one in the list
      max_column = cl[which(column_sums==max(column_sums))][1]
      # add all other columns to that one
      for (id in cl[cl != max_column]) {
        ab[, max_column] = ab[, max_column, drop=F] + ab[, id, drop=F]
        filtered_columns = c(filtered_columns, id)
      }
    }
  }
  ab = ab[, setdiff(1:ncol(ab), sort(filtered_columns)), drop=F]

  ab_colsums = unname(colSums(ab))
  if (verbose) cat('OK.\n')

  # Cleanup the temp database
  if (verbose) cat('* cleaning up temporary files...')
  file.remove(blast_out)
  file.remove(fasta)
  cleanup_blastdb(db)
  if (verbose) cat('OK.\n')

  # Output abundance table sorted by total abundance (add up samples) in
  # descending order.
  return(ab[, order(-ab_colsums), drop=F])
}


ab_to_files = function (ab, verbose=T, temp_dir=tempdir()) {
   # Convert the abundance table matrix into a FASTA file, together
   # with the read counts. Internal function used by collapse().
   meta = colnames(ab)
   colnames(ab) = 1:ncol(ab)

   rand_id = sample(LETTERS, 10)
   # Generate a temp fasta file and blast database prefix
   fasta_out = file.path(
      temp_dir,
      paste(c('collapse_', rand_id, '.fasta'), collapse='')
   )
   db_prefix = file.path(dirname(fasta_out), paste(c('collapse_', rand_id), collapse=''))

   # Generate a temp abundance table.
   tab_out = file.path(
      temp_dir,
      paste(c('collapse_', rand_id, '.txt'), collapse='')
   )

   tab.dt = data.table(ab)
   tab.dt[, sample_id := rownames(ab)]
   tab_m.dt = data.table::melt(tab.dt, variable.name='seq_id', value.name='count',
                   id.vars='sample_id')

   # Write a temp fasta
   fasta_writer(colnames(ab), meta, fasta_out)

   # Generate a BLAST database from this fasta file
   makeblastdb(fasta_out, db_prefix, verbose=verbose)

   # Return file names
   return(c(fasta_out, db_prefix))
}

