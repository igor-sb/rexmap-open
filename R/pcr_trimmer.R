



#' PCR primer trimmer
#'
#' Trims PCR primers from FASTQ sequences.
#'
#' @param fq_in A character vector of input file names.
#' @param fq_out A character vector of output file names.
#' @param region Hypervariable region used. Used to automatically retrieve PCR
#' primer sequences for trimming. If this parameter is used (not NULL), then
#' \code{pr_fwd} and \code{pr_rev} are ignored.
#' @param pr_fwd Sequence for the \strong{forward} PCR primer (5' -> 3') in the extended
#' nucleotide code.
#' @param pr_rev Sequence for the \strong{reverse} PCR primer (5' -> 3') in the extended
#' nucleotide code. If you already have 3' -> 5' sequence, use \code{\link{reverse_complement}}
#' to obtain the 5' -> 3' sequence.
#' @param pr_fwd_maxoff Maximum allowed offset from the 5' sequence end for the
#' forward PCR primer. Primers aligned further into the read ignored.
#' @param pr_rev_maxoff Same but for reverse PCR primer, from the 3' end.
#' @param ncpu Number of CPU cores (threads) to use for multithreading. By default
#' all available CPU cores are used (\code{parallel::detectCores()}).
#' @param max_mismatch Maximum allowed number of mismatches for each PCR primer,
#' ignoring any ambiguous nucleotides (not A,C,G or T).
#' @param timing Display run-time at the end of trimming (default: FALSE).
#' @param overwrite Overwrite target files if they exist (default: TRUE)
#'
#' @export
remove_pcr_primers = Vectorize(function (fq_in, fq_out, region=NULL,
                              pr_fwd=NULL,
                              pr_rev=NULL,
                              pr_fwd_maxoff=10, pr_rev_maxoff=10,
                              return_noprimer=T,
                              ncpu=rexmap_option('ncpu'), max_mismatch=2,
                              timing=F, verbose=rexmap_option('verbose'),
                              overwrite=TRUE) {
  # fq_in = input fastq file
  # fq_out = output fastq file (without primers)
  # pr_fwd = forward primer 5'->3'
  # pr_rev = reverse primer 5'->3'
  #
  if (verbose) cat('FASTQ input: ', fq_in, fill=T)
  if (file.size(fq_in) == 0) {
    # Input file has zero size, so no sequences. This can happen if all the reads
    # are too noisy and have been filtered out in the merging step, if we were
    # unable to merge anything.
    file.create(fq_out)
    return(c('fwd_trim'=0, 'rev_trim'=0))
  }

  # Check input
  if (!all_exist(fq_in)) stop('PCR primer remover: some input files are missing.')
  if (is.null(region) & (is.null(pr_fwd) | is.null(pr_rev))) {
    stop('PCR primer remover: either region of pr_fwd and pr_rev need to be specified.')
  }

  if (!overwrite & file.exists(fq_out)) {
    # We are not overwriting files and output already exists
    return(c('fwd_trim'=NA, 'rev_trim'=NA))
  }

  # if region is specified, check that primers exist in the reference table
  if (!is.null(region)) {
    if (nrow(rexmap_option('blast_dbs')[Hypervariable_region==region]) == 0) {
      # Missing hypervariable region
      stop('PCR primer remover: hypervariable region \"', region, '\" not found.')
    }
    pr_fwd = rexmap_option('blast_dbs')[Hypervariable_region==region, Primer1_sequence_5to3]
    pr_rev = reverse_complement(
      rexmap_option('blast_dbs')[Hypervariable_region==region, Primer2_sequence_3to5]
    )
    if (verbose) cat('* PCR primer remover mode: region', region, paste0('(fwd: ',
      pr_fwd, ', rev: ', pr_rev, ')'), fill=T)
  } else {
    if (verbose) cat('* PCR primer remover mode: primers', paste0('(fwd: ',
      pr_fwd, ', rev: ', pr_rev, ')'), fill=T)

  }

  # Check if all output folders exist, and if not create them
  output_folders = unique(dirname(fq_out))
  for (out_f in output_folders) {
    if (!dir.exists(out_f)) {
      if (verbose) cat('* created output folder.', fill=T)
      dir.create(out_f, recursive=T)
    }
  }

  # Load Input FASTQ file
  if (timing) start_time = Sys.time()
  # Count extended DNA symbols
  # Ignore any extended DNA symbols. Shouldn't have too many of them anyway.
  pr_fwd = gsub('[^ACGT-]', 'N', pr_fwd)
  pr_rev = gsub('[^ACGT-]', 'N', pr_rev)
  if (verbose) cat('* loading file...')
  in_fq = sfastq_reader(fq_in)
  if (verbose) cat('OK.', fill=T)
  # max_mismatch = 2

  seq = in_fq[['seqs']][1]
  qual = in_fq[['qual']][1]

  fastq_trimmer = function (meta, seq, qual, return_noprimer=return_noprimer) {

    # Search for forward primer.
    aln = C_nwalign(pr_fwd, seq, match=1, mismatch=-1, indel=-1)

    aln_left_offset = regexpr('[^-]', aln[2])[1]
    pr_fwd_left = max(regexpr('[^-]', aln[1])[1],
                      aln_left_offset)
    # pr_fwd_left = regexpr('[^-]', aln[1])[1]
    pr_fwd_right = regexpr('[-]{1,}$', aln[1])[1]-1
    pr_fwd_n = lengths(regmatches(pr_fwd, gregexpr('N', pr_fwd)))
    pr_rev_n = lengths(regmatches(pr_rev, gregexpr('N', pr_rev)))
    # Alignment statistics: match, mismatch, gapopen, gapextend
    aln_stat = compare_alignment(
      str_sub(aln[1], start=pr_fwd_left, end=pr_fwd_right),
      str_sub(aln[2], start=pr_fwd_left, end=pr_fwd_right)
    )
    # Forward primer alignment is acceptable if it's near beginning (within first 5 nts)
    # and if it doesn't have more than 2 mismatches for non-N symbols, which includes indels.
    pr_fwd_found = FALSE
    pr_rev_found = FALSE
    if (pr_fwd_left < pr_fwd_maxoff & aln_stat[2]+aln_stat[3]+aln_stat[4]-pr_fwd_n <= max_mismatch) {
      pr_fwd_found = TRUE
    }

    # Search for reverse primer
    aln2 = C_nwalign(pr_rev, seq, match=1, mismatch=-1, indel=-1)
    pr_rev_left = regexpr('[^-]', aln2[1])[1]
    pr_rev_right = min(regexpr('[ACGTN][^ACGTN]*$', aln2[1])[1],
                       regexpr('[ACGTN][^ACGTN]*$', aln2[2])[1])
    aln2_stat = compare_alignment(
      str_sub(aln2[1], start=pr_rev_left, end=pr_rev_right),
      str_sub(aln2[2], start=pr_rev_left, end=pr_rev_right)
    )
    if (pr_rev_left > nchar(seq)-nchar(pr_rev)-pr_rev_maxoff & aln2_stat[2]+aln2_stat[3]+aln2_stat[4]-pr_rev_n <= max_mismatch) {
      # Reverse primer found
      pr_rev_found = TRUE
    }

    # Trim if needed
    start = 1L
    end = -1L
    if (pr_fwd_found) start = pr_fwd_right - aln_left_offset + 2
    if (pr_rev_found) end = pr_rev_left - 1

    return(list('meta' = meta,
                'seqs' = str_sub(seq, start=start, end=end),
                'qual' = str_sub(qual, start=start, end=end),
                'trim_fwd' = pr_fwd_found,
                'trim_rev' = pr_rev_found))
  }

  # Apply fastq_trimmer() to each sequence in this file
  if (verbose) cat('* removing primers...')
  ncpus = min(ncpu, length(seq))
  out_trimmed = unname(
    parallel::mcmapply(
      fastq_trimmer, in_fq[['meta']], in_fq[['seqs']], in_fq[['qual']],
      mc.cores=ncpus, SIMPLIFY=FALSE
    )
  )
  if (verbose) cat('OK.', fill=T)

  fwd_trimmed = sum(sapply(out_trimmed, function (x) x$trim_fwd))
  rev_trimmed = sum(sapply(out_trimmed, function (x) x$trim_rev))

  # Save results in a new file
  if (verbose) cat('* saving output...')
  fastq_list_writer(out_trimmed, fq_out, ncpu=ncpu)
  if (verbose) cat('OK.', fill=T)
  if (timing) {
    end_time = Sys.time()
    diff_time = difftime(end_time, start_time, units='secs')
    cat('Finished in ', round(as.numeric(diff_time)/60), ' m ',
      round(as.numeric(diff_time)%%60, 1), ' s.\n')
  }
  return(c('fwd_trim'=fwd_trimmed, 'rev_trim'=rev_trimmed))
}, vectorize.args=c('fq_in', 'fq_out'))






