#' Merge paired-end reads
#'
#' Merges reads using overlap version of the Needleman-Wunsch alignment algorithm and calculate
#' merged quality scores based on the posterior probabilities (cite: Edgar, Fyelberg)
#'
#' @param fq_fwd A character vector for input forward FASTQ files (can be gzipped).
#' @param fq_rev A character vector for input reverse FASTQ files (can be gzipped).
#' @param fq_mer A character vector for output merged FASTQ files.
#' @param min_sim Minimum similarity to accept an alignment. A floating point number between 0 and 1. Default: 0.75, corresponding to minimum 75\% similarity in the aligned overlap region.
#' @param min_aln_len Minimum alignment length. Ignore alignments below this alignment length.
#' @param match Score for character match in the alignment (default: 5.
#' @param mismatch Score for character mismatch in the alignment (default: -5)
#' @param gap_p Score for gap penalty, either gap opening on gap extension (default: -7).
#' @param rc_reverse TRUE/FALSE Reverse complement reverse reads? (default: TRUE)
#' @param ncpu Number of CPU threads to use for multithreading.
#' @param verbose TRUE/FALSE Display of status messages.
#' @param timing TRUE/FALSE Time merging.
#'
#' @export
merge_pairs = Vectorize(function (fq_fwd, fq_rev, fq_mer, min_sim=0.75, min_aln_len=50,
                              match=5L, mismatch=-5L, gap_p=-7L, rc_reverse=TRUE,
                              ncpu = rexmap_option('ncpu'),
                              verbose=FALSE, timing=FALSE
                              ) {
  # fq_fwd = vector of filenames (incl. paths) to forward reads
  # fq_rev = vector of filenames (incl. paths) to reverse reads
  # fq_mer = vector of filenames (incl. paths) to write merged reads
  # min_pct_sim = minimum similarity of the overlap aligned region, *100 to get pct (default: 0.75)
  # min_aln_len = minimum overlap aligned region length in nts (default: 50)
  # By default it overwrites existing fq_mer files.
  # Ends-free alignment parameters
  #  match = (int) score for matching ACGT (default: 5)
  #  mismatch = score for mismatching ACGT or alignment vs N (default: -2)
  #  gap_p = gap opening and extension penalty (default: -7)
  #  rc_reverse = (bool) reverse complement reverse reads before alignment? (default: TRUE)
  #  threads = number of cores to use for multi-threading
  #   (tries to detect total number of cores using parallel::detectCores())
  # verbose = (bool) print status messages? (default: FALSE)
  # timing = (bool) print execution time to stderr (default: FALSE)
  # path_precomputed_posterior = path to tab-delimited headerless tables with precomputed
  #  quality scores for posterior probabilities for read merging
  #  output from mergepairs_generate_posterior_probabilities.py
  if (timing) start_time = Sys.time()

    if (is.na(ncpu) | ncpu == 0) {
      # If we can't figure out number of threads, set it to 2 or just stop?
      ncpu = 2
      stop('Merge pairs: invalid number of CPU threads.')
  }
  m = function (...) {
    # Print message if verbose is enabled.
    if (verbose) {
      message(..., appendLF=F)
    }
  }

  # Does output folder exist?
  output_folders = unique(dirname(fq_mer))
  for (out_f in output_folders) {
    if (!dir.exists(out_f)) dir.create(out_f, recursive=T)
  }

  # Load reads using ShortRead::FastqStreamer
  m(paste0('Loading FASTQ reads: ', basename(fq_fwd), ', ', basename(fq_rev), ' ...'))
  f_fwd = ShortRead::FastqStreamer(fq_fwd)
  f_rev = ShortRead::FastqStreamer(fq_rev)
  # Go through each read entry, do alignment
  r_fwd = ShortRead::yield(f_fwd)
  if (length(r_fwd) == 0) break
  if (rc_reverse) r_rev = ShortRead::reverseComplement(ShortRead::yield(f_rev))
  else r_rev = ShortRead::yield(f_rev)

  # Process chunk
  read_fwd = as.character(ShortRead::sread(r_fwd))
  read_rev = as.character(ShortRead::sread(r_rev))
  qual_fwd = as.character(Biostrings::quality(Biostrings::quality(r_fwd)))
  qual_rev = as.character(Biostrings::quality(Biostrings::quality(r_rev)))
  ids = gsub('^([^ ]+) .*', '\\1', as.character(ShortRead::id(r_fwd)))
  m(' OK.\n')

  # Apply mergepairs
  m('Merging pairs...')
  merged_list = parallel::mcmapply(C_mergepairs, read_fwd, read_rev, qual_fwd, qual_rev,
                         match=match, mismatch=mismatch, gap_p=gap_p,
                         min_pct_sim=min_sim, min_aln_len=min_aln_len,
                         posterior_match_file=rexmap_option('mergepairs_matchqs'),
                         posterior_mismatch_file=rexmap_option('mergepairs_mismatchqs'),
                         mc.cores=ncpu
                       )
  m(' OK.\n')
  merged_aln_filter = as.logical(unname(merged_list[3, ])) &
    as.logical(unname(merged_list[4, ]))

  # Free memory
  close(f_fwd)
  close(f_rev)
  rm(read_fwd, read_rev, qual_fwd, qual_rev, f_fwd, f_rev, r_fwd, r_rev)
  m('Writing output files')
  # sequences
  final_seqs = unname(merged_list[1, merged_aln_filter])
  m('.')
  final_qual = unname(merged_list[2, merged_aln_filter])
  m('.')
  final_ids  = ids[merged_aln_filter]
  m('.')

  final_filter = !is.null(final_seqs) & !is.null(final_qual)
  if (!all(final_filter)) m('.')
  final_seqs = final_seqs[final_filter]
  final_qual = final_qual[final_filter]
  final_ids  = final_ids[final_filter]

  # any NULL sequences? (no idea why/how this happens sometime) seems like a
  # multithreading issue

  merged_sread = ShortRead::ShortReadQ(
    sread = Biostrings::DNAStringSet(final_seqs),
    quality = Biostrings::BStringSet(final_qual),
    id = Biostrings::BStringSet(final_ids)
  )

  # Generate statistics for filtered-out reads
  stats = list(
    'total'=length(merged_aln_filter),
    'low_pct_sim'=length(merged_aln_filter) - sum(as.logical(merged_list[3, ])),
    'low_aln_len'=length(merged_aln_filter) - sum(as.logical(merged_list[4, ]))
  )

  # Free memory and close file connections
  rm(merged_list)

  # If file exists, delete it, then write new one
  if (file.exists(fq_mer)) file_remove_result = file.remove(fq_mer)

  ShortRead::writeFastq(merged_sread, fq_mer, compress = F)
  m(' OK.\n')
  if (timing) {
    end_time = Sys.time()
    diff_time = difftime(end_time, start_time, units='secs')
    m('Finished in ', round(as.numeric(diff_time)/60), ' m ',
      round(as.numeric(diff_time)%%60, 1), ' s.\n')
  }
  return(stats)
}, c('fq_fwd', 'fq_rev', 'fq_mer'))

#' Convert mergestats table to a normal data.table
#'
#' @export
mergeout_to_table = function (mergestats) {
  mergestats.dt = data.table::data.table(
    matrix(
      unlist(t(mergestats)),
      ncol=nrow(mergestats),
      dimnames=list(sampleids_from_filenames(colnames(mergestats), '_'), rownames(mergestats))
    ),
    keep.rownames='sample_id'
  )
  return(mergestats.dt)
}
