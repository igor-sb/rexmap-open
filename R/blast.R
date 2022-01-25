# Run BLAST with denoised sequences from dada2 against a reference FASTA file
# Record all matches.

# We will use this to check how many exact matches, one off etc we have against
# sequences in the reference file.


blast_best_seq_matches = function (dt, id_col='dada_seqid', match_col='match_strains',
                                   exclude_match=FALSE) {
  # For each strain sequence, find the best match (by bitscore) for
  # fwd and rev primers.
  # exclude_match=T will not show the match_col in the output and will just show
  # statistics for the best matched hit.
  dt2 = dt[, {
    if (.N > 0) {
      best_match = .SD[score == max(score)]
      if (exclude_match) {
        unique(best_match[, setdiff(names(best_match), match_col), with=F])
      } else {
        best_match
      }
    }
  }, by=c(id_col)]
  return(dt2)
}

blast_out_to_best_cp = function (
   blast_output, region=NULL, ref_cp=NULL, aln_params=rexmap_option('aln_params'),
   blast_out_fmt=rexmap_option('blast_out_fmt'),
   verbose=T, ncpu=rexmap_option('ncpu'), variant_sep='-',
   show_alignment=F, skip_cp=F
  )
{
  # Select region all ref_cp from file
  if (!is.null(region)) { # If region is given, ignore ref_cp
    ref_cp = rexmap_option('blast_dbs')[Hypervariable_region==region, table]
    ref_cp = system.file('database', ref_cp, package=pname_l)
    if (length(ref_cp) == 0) stop('blast: invalid hyper-variable region.')
  } else { # Region is not given, just check if ref_cp file exist
    if (is.null(ref_cp)) stop('blast: reference copy number table argument missing.')
    if (!file.exists(ref_cp)) stop('blast: reference copy number table file missing.')
  }

  if (verbose) cat('* blast out: ')
  blast_out.dt = data.table::fread(blast_output)
  names(blast_out.dt) = strsplit(blast_out_fmt, ' ')[[1]]
  # Calculate BLAST scores.
  blast_out.dt[, c('score', 'match', 'mismatch', 'gapopen', 'gapextend') := data.table::transpose(
    unname(parallel::mcmapply(function (q,s) {
    aln = compare_alignment(q,s)
    as.integer(c(sum(aln_params*aln), aln[1], aln[2], aln[3], aln[4]))
  }, qseq, sseq, SIMPLIFY=F, mc.cores=ncpu)))]

  # For local alignments that are one-off at each end, reduce the score using
  # a mismatch penalty. (We deliberately use indel score that is very similar.)
  blast_out.dt[qstart >= 2 & sstart >= 2, score := score-aln_params[2]]
  blast_out.dt[qstart >= 2 & sstart >= 2, mismatch := mismatch+1L]
  blast_out.dt[qend <= qlen-1 & send < slen, score := score-aln_params[2]]
  blast_out.dt[qend <= qlen-1 & send < slen, mismatch := mismatch+1L]

  # Filter out alignments with large gaps at beginning or the end
  # blast_out.dt = blast_out.dt[qstart <= 2 & qend >= qlen-1]
  if (verbose) cat('OK. blast best: ')

  # For each dada2 sequence, keep only the best matches by bitscore
  blast_best.dt = blast_best_seq_matches(
    blast_out.dt,
    strsplit(blast_out_fmt, ' ')[[1]][1],
    strsplit(blast_out_fmt, ' ')[[1]][2]
  )
  blast_best.dt[, variant_id := gsub('^([0-9]+)-.*', '\\1', sseqid)]
  # Remove sequence alignment columns since we already calculated alignment score
  if (!show_alignment) {
     blast_best.dt[, c('qseq', 'sseq') := NULL]
  }
  # blast_best.dt[, c('qseq') := NULL]
  blast_best.dt[, c('sseqid') := NULL]
  # Calculate percentage similarity
  blast_best.dt[, pctsim := round(100*match/(match+mismatch+gapopen+gapextend), 2)]

  if (skip_cp) {
    return(list(blast_best.dt, data.table()))
  }

  if (verbose) cat('OK. copy number table: ')
  cp.dt = data.table::fread(ref_cp, colClasses=c('character', 'character', 'integer'))

  # Prepare copy number table columns
  cp.dt[, strain := gsub('_@rrn[0-9]+', '', strain_name)]

  # Generate strain index
  cp.dt[, strain_id := .GRP, by=strain]
  cp.dt[, strain_name := NULL]

  # How many unique 16S rrna variants of each strain?
  cp.dt[, rrn_uniq := length(unique(variant_id)), by=strain]

  # Sort by strain_id, then by
  cp.dt = cp.dt[order(strain_id, variant_id)]
  if (verbose) cat('OK.\n')
  # Merge with blast hits
  # First do perfect matches

  if (verbose) cat('merge: ')
  blast_best.dt = merge(blast_best.dt, cp.dt[, .(variant_id, strain, copy_number)],
                        by='variant_id', allow.cartesian=T)
  if (verbose) cat('OK. Fix overhang differences:')
  hang_diff_qseqids = blast_best.dt[pctsim==100, 1, by=qseqid][, unique(qseqid)]
  for (q in hang_diff_qseqids) {
    vids = blast_best.dt[pctsim==100 & qseqid==q, sort(unique(variant_id))]
    blast_best.dt[pctsim==100 & qseqid==q, variant_id := paste(vids, collapse=variant_sep)]
    cp.dt[variant_id %in% vids, variant_id := paste(vids, collapse=variant_sep)]
  }
  if (verbose) cat('.')
  blast_best.dt[pctsim==100, copy_number := sum(copy_number), by=.(variant_id, strain)]
  blast_best.dt = unique(blast_best.dt)

  if (verbose) cat('.')
  # Add up copy numbers for strain that differ in the undetected overhang
  cp.dt[, copy_number := sum(copy_number), by=.(variant_id, strain_id)]
  cp.dt = unique(cp.dt)
  if (verbose) cat('OK. ')


  # Recalculate rrn_uniq after re-numerating variant_ids (b/c the slightly shorter sequence now has the same
  # or smaller number of variants)
  cp.dt[, rrn_uniq := length(unique(variant_id)), by=strain_id]

  # Group sequences into OSU
  cp.dt[, spectrum := paste(copy_number, variant_id, sep=':', collapse=','), by=strain]

  # For strains which have only one unique 16S sequence variant, group them into the same
  # same OSU, regardless of the copy number.
  cp.dt[rrn_uniq==1, spectrum := paste0('X:', strsplit(spectrum[1], ':')[[1]][2]), by=strain]

  # OSU is a group of strains all with the same strain spectrum
  cp.dt[, osu_id := .GRP, by=spectrum]

  # Instead of X put mean copy number or 1
  cp.dt[rrn_uniq==1, spectrum := gsub('X', round(mean(copy_number)), spectrum), by=osu_id]
  # ref_cp.dt[rrn_uniq==1, spectrum := gsub('X', '1', spectrum), by=osu_id]

  # Count number of strains binned together
  cp.dt[, no_strains_in_osu := length(unique(strain)), by=osu_id]

  # Generate species names
  cp.dt[, species := gsub('^([^_]+)_([^_]+)_.*', '\\1_\\2', strain)]
  # Count unique number of species per OSU
  # cp.dt[, no_species_in_osu := length(unique(species)), by=osu_id]
  if (verbose) cat('OK.\n')
  return(list(blast_best.dt, cp.dt))
}


# OSU binning and abundance estimation
blast_cp_to_osu_dt = function (
  blast_best.dt, cp.dt, ab_tab_nochim_m.dt,
  ncpu = rexmap_option('ncpu'),
  verbose=T
) {
  pctsim_min = 100.00 # Min. pct sim for OSU binning
  if (verbose) cat('OSU table: ')
  # First generate a BLAST best results table with just alignment pctsim
  blast_reduced.dt = unique(blast_best.dt[, .(variant_id, qseqid, pctsim)])[order(qseqid)]
  blast_reduced.dt[, variant_id_new := variant_id[1], by=qseqid]
  # blast_reduced.dt = merge(blast_reduced.dt, ab_tab_nochim_m.dt[, .(qseqid, raw_count, sample_id)],
  #                          by='qseqid', all.x=T)
  blast_reduced.dt = blast_reduced.dt[, .SD[pctsim==max(pctsim)][1], by=variant_id]

  # Generate a separate reference table
  osu.dt = unique(cp.dt[, .(osu_id, spectrum, strain, no_strains_in_osu)])
  osu.dt[, no_variants := sapply(spectrum, function (s) length(strsplit(s, ',', fixed=T)[[1]]))]
  # osu_m.dt = osu.dt[, .(variant_id = sub('.*:([0-9]{5})$', '\\1', strsplit(spectrum, ',')[[1]]),
  #                       spectrum=spectrum, strain, no_variants), by=osu_id]

  # Now generate a smaller OSU table containing only osu id that contain
  # at least one variant ID from the data with > 0 raw count.
  varids_data = blast_reduced.dt[pctsim >= pctsim_min, unique(variant_id)]
  # varids_data = [pctsim > pctsim_min, unique(variant_id)]
  osu_data.dt = osu.dt[, {
    varids = strsplit(gsub('[0-9]+:', '', spectrum), ',')[[1]]
    overlap = intersect(varids_data, varids)
    if (length(overlap) > 0) .SD
  }, by=osu_id]

  if (verbose) cat('OK. Melt table: ')
  sample_ids = ab_tab_nochim_m.dt[, unique(sample_id)]

  # Optimized osu_data_m2.dt
  osu_data2.dt = copy(osu_data.dt)
  osu_data2.dt = osu_data.dt[,
   data.table::transpose(
      strsplit(strsplit(spectrum, ',', fixed=T)[[1]], ':', fixed=T)),
   by=.(osu_id, strain, no_strains_in_osu, no_variants)]
  names(osu_data2.dt)[5:6] = c('copy_number', 'variant_id')
  osu_data2.dt[, copy_number := as.integer(copy_number)]
  # Add BLAST alignment information to each variant_id
  osu_data3.dt = merge(
     osu_data2.dt,
     blast_reduced.dt[pctsim >= pctsim_min, .(variant_id, qseqid, pctsim)],
     by='variant_id', all.x=T)

  # Remove strain annotation here
  osu_data4.dt = copy(osu_data3.dt)
  osu_data4.dt[, strain := NULL]
  osu_data4.dt = unique(osu_data4.dt)

  # Fill in zeros for some OSUs where we have some, but not all variants
  sample_ids = ab_tab_nochim_m.dt[, unique(sample_id)]
  ab_tab_nochim_m_fill.dt = data.table::rbindlist(list(
     ab_tab_nochim_m.dt,
     data.table(sample_id=sample_ids, qseqid=as.integer(NA), raw_count=0L)
  ))

  # Add abundance information
  osu_data5.dt = merge(
     osu_data4.dt,
     ab_tab_nochim_m_fill.dt,
     by='qseqid', all.x=T, allow.cartesian=T
  ) # allow cartesian is due to multiple samples having same qseqids

  osu_data5.dt = osu_data5.dt[, .(osu_id, copy_number, variant_id, raw_count, sample_id, pctsim, no_strains_in_osu)]
  osu_data5.dt = osu_data5.dt[order(sample_id, osu_id, variant_id)]

  if (verbose) cat(' OK.')
  return(osu_data5.dt[, .(osu_id, copy_number, variant_id, raw_count, sample_id, pctsim, no_strains_in_osu)])
}

#' BLAST sequences against a reference database
#'
#' Aligns each sequence from the FASTA file, a character vector of DNA sequences,
#' or an abundance table, against a reference database.
#' Either \code{region} or both \code{ref_db} and \code{ref_cp} arguments must be specificed.
#' If \code{region}
#' is specified, then the \code{ref_db} and \code{ref_cp} arguments are ignored and a
#' reference database is chosen
#' from the pre-computed set for that hyper-variable region. For available databases
#' see columns DB and table in \code{rexmap_option('blast_dbs')}.
#' This function returns a blast-class object (named list; see ?blast-class for info).
#'
#' @param sequences (Required) Either a FASTA file, a character vector of DNA sequences, or
#' an abundance table (output from \code{\link{sequence_abundance}})
#' @param region Hyper-variable region. If this is not NULL, \code{ref_db} is ignored. Possible
#' values: 'V4', 'V3-V4'.
#' @param ref_db Full path to a custom database.
#' @param ref_cp Full path to a copy number table for a custom database.
#' @param max_target_seqs Maximum number of target sequences to return for each query.
#' Note that the temporary output file can be large if this is set > 1000.
#' @param word_size Word size used for BLAST alignment.
#'
#' @export
blast = function (sequences, blast_output=NULL, region=NULL, ref_db=NULL,
                  ref_cp=NULL, max_target_seqs=rexmap_option('blast_max_seqs'),
                  word_size=rexmap_option('blast_word_size'),
                  verbose=rexmap_option('verbose'),
                  show_args=F, output_error=F, show_alignment=F,
                  ncpu=rexmap_option('ncpu'), temp_dir=tempdir()) {

  # Pre-blastn sequence argument check
  # Sequences can be either FASTA file (ends with either .fa or .fasta),
  # an abundance data table (output from sequence abundance) or a character vector
  # that is a list of sequences.

  # If blast_output is not given, then we will generate it, so mark it for cleanup
  cleanup_blast_output = FALSE
  if (is.null(blast_output)) {
    cleanup_blast_output = TRUE
  }

  sequences_type = NULL
  if ('data.table' %in% class(sequences)) {
    if ('sequence' %in% names(sequences)) {
      # It's an abundance table. Create temporary FASTA file to run BLAST then
      # delete it afterwards.
      sequences_type = 'dt'
      rand_id = sample(LETTERS, 10)
      fasta_file = file.path(temp_dir, paste(c('blast_', rand_id, '.fasta'), collapse=''))
      sequences_to_fasta(sequences, fasta_file)
      if (verbose) cat('* blast input type: abundance table', fill=T)
    } else {
      stop('blast: input abundance table does not have \"sequence\" column.')
    }
  } else {
    # First check that the class is character
    if (class(sequences) == 'character') {
      if (length(sequences) > 1) {
        # Looks like a set of DNA sequences. Check letters.
        sequences_type = 'DNA'
        rand_id = sample(LETTERS, 10)
        fasta_file = file.path(temp_dir, paste(c('blast_', rand_id, '.fasta'), collapse=''))
        sequences_to_fasta(data.table(qseqid=1:length(sequences),
                                      sequence=sequences), fasta_file)
        if (verbose) cat('* blast input type: character vector', fill=T)
      } else {
        if (grepl('[ACGTN]+', sequences)) {
          # Looks like a single sequence
          sequences_type = 'DNA'
          rand_id = sample(LETTERS, 10)
          fasta_file = file.path(temp_dir, paste(c('blast_', rand_id, '.fasta'), collapse=''))
          sequences_to_fasta(data.table(qseqid=1:length(sequences),
                                        sequence=sequences), fasta_file)
          if (verbose) cat('* blast input type: character vector', fill=T)
        } else {
          # Looks like a FASTA file. Just check that it exists.
          if (!grepl('\\.f[n]?a|\\.fasta', sequences)) stop('blast: file needs to be in FASTA format. If it is, use normal file extensions: .fa and .fasta.')
          if (!file.exists(sequences)) stop('blast: sequence FASTA file does not exits.')
          fasta_file = sequences
          sequences_type = 'fasta'
          if (verbose) cat('* blast input type: FASTA file', fill=T)
        }
      }
    } else {
      # The class is neither character or data.table / data.frame.
      stop('blast: invalid sequences format.')
    }
  }

  # Run Blastn
  # Run blast_out_to_best_cp
  # Return himap object
  if (is.null(blast_output)) {
    # Generate temporary blast output file
    rand_id = sample(LETTERS, 10)
    blast_output = file.path(temp_dir, paste(c('blast_output_', rand_id, '.txt'),
                                              collapse=''))
  }
  if (!is.null(region)) {
    blast_status = blastn(fasta_file, blast_output, region=region,
                          max_target_seqs=max_target_seqs,
                          word_size=word_size, output_err=output_error,
                          ncpu=ncpu, verbose=verbose)
  } else {
    blast_status = blastn(fasta_file, blast_output, ref_db=ref_db,
                          max_target_seqs=max_target_seqs,
                          word_size=word_size, output_err=output_error,
                          ncpu=ncpu, verbose=verbose)
  }
  if (blast_status != 0) stop('blast: error running blastn.')
  # Load BLAST results (can take a while if there are lots of sequences and max_target_seqs
  # is large.
  if (!file.exists(blast_output)) stop('blast: ', blast_output, ' file does not exist.')
  if (!is.null(region)) {
    blast_cp = blast_out_to_best_cp(blast_output, region=region,
                                    show_alignment=show_alignment, verbose=verbose)
  } else {
    blast_cp = blast_out_to_best_cp(blast_output, ref_cp=ref_cp,
                                    show_alignment=show_alignment)
  }

  names(blast_cp) = c('alignments', 'cp')
  blast_cp$parameters = list(
    'max_target_seqs'=max_target_seqs, 'word_size'=word_size, 'alignment_parameters'=
    paste(c('match', 'mismatch', 'gap_open', 'gap_extend'), rexmap_option('aln_params'),
          collapse=', ', sep=': '))

  if (class(sequences)[1] == 'character') {
    return(as(blast_cp, 'blast'))
  }

  blast_cp$sequences = unique(blast_cp$cp[, .(osu_id, variant_id, spectrum)])
  blast_cp$sequences = merge(
    blast_cp$sequences,
    unique(blast_cp$alignments[pctsim==100, .(variant_id, qseqid)])[, .(qseqid=min(qseqid)), by=variant_id],
    by='variant_id', all.x=T
  )
  # The last part min(qseqid) selects one of the qseqids with the SAME variatn
  # this can only happen if two different sequences align EXACTLY to the same
  # variant, which only happens if they are reverse complement of each
  # other.
  blast_cp$sequences = data.table::rbindlist(list(
    blast_cp$sequences,
    unique(blast_cp$alignments[
      pctsim<100,
      .(variant_id=NA, osu_id=qseqid+1000000L, spectrum=NA, qseqid)])
  ))
  blast_cp$sequences[, copy_number := 1L]
  blast_cp$sequences[osu_id < 1000000L, copy_number := as.integer(
    gsub('^([0-9]+):.*', '\\1', strsplit(spectrum[1], ',', fixed=T)[[1]])
  ), by=osu_id]

  blast_cp$sequences = merge(
    blast_cp$sequences,
    unique(sequences[, .(qseqid, sequence)]),
    by='qseqid', all.x=T)


  # If we made a temp fasta file, remove it
  if (sequences_type %in% c('dt', 'DNA')) {
    file.remove(fasta_file)
  }
  # If we made temp blast file, remove it as well
  # if (dirname(blast_output) == temp_dir) file.remove(blast_output)
  if (cleanup_blast_output) {
    if (file.exists(blast_output)) {
      file.remove(blast_output)
    }
  }
  return(as(blast_cp, 'blast'))
}


blastn = function (
  seqs_fa, output, region=NULL, ref_db=NULL,
  blast_path=rexmap_option('path_blastn'),
  match=rexmap_option('aln_params')[1], mismatch=rexmap_option('aln_params')[2],
  gapopen=-rexmap_option('aln_params')[3], gapextend=-rexmap_option('aln_params')[4],
  word_size=rexmap_option('blast_word_size'), ncpu=rexmap_option('ncpu'),
  max_target_seqs=rexmap_option('blast_max_seqs'),
  perc_identity=75, outfmt=paste0('6 ', rexmap_option('blast_out_fmt')),
  dust='20 64 1',
  output_err=F, verbose=rexmap_option('verbose')
  ) {

  dbs = rexmap_option('blast_dbs')
  if (!is.null(region)) {
    # Region is specificed, ignore ref_db
    if (!(region %in% dbs[, Hypervariable_region])) {
      stop('blast: wrong hypervariable region specified.')
    } else {
      if (dbs[Hypervariable_region==region, DB] == '') {
        stop('blast: no database listed for this hypervariable region.')
      } else {
        ref_db = system.file('database',
                             paste0(dbs[Hypervariable_region==region, DB], '.nhr'),
                             package=pname_l)
        if (ref_db == '') {
          stop('blast: missing database.')
        } else {
          ref_db = sub('.nhr$', '', ref_db)
        }
      }
    }
  } else {
    # ref_db is specified so just use that and check if the files are there.
    missing_files = file.exists(paste0(ref_db, c('.nin', '.nhr', '.nsq')))
    if (any(!missing_files)) {
      stop('blast: missing database file(s).')
    }
  }
  blast_args = c(
      '-dust', shQuote(dust), '-word_size', word_size,
      '-reward', match, '-penalty', mismatch, '-gapopen', gapopen,
      '-gapextend', gapextend, '-outfmt', shQuote(outfmt),
      '-query', seqs_fa, '-db', ref_db, '-num_threads', ncpu,
      '-perc_identity', perc_identity
  )
  if (!is.null(max_target_seqs)) {
    blast_args = c(blast_args, '-max_target_seqs', max_target_seqs)
  }

  if (output_err != F) cat(paste(blast_args, sep=' '), fill=T)
  if (verbose) {
    cat('blastn: ', paste0(blast_args, collapse=' '), '\n')
  }
  x = system2(blast_path, args=blast_args, stdout=output, stderr=output_err)
  return(x)
}
