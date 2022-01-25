#' Pool multiple RExMap outputs
#'
#' Combine RExMap outputs from multiple datasets into a single table.
#'
#' @param result_list A \code{list()} with each element being a list of two
#' data tables.
#'
#' 1. OSU abundance data table (output from \code{\link{abundance}}) with
#' columns:
#'
#' \code{sample_id | osu_id | osu_count | pctsim | species}\cr
#'
#' Additionally, this data table can contain an extra \code{dataset} column
#' which contains a name of that dataset.
#'
#' 2. OSU sequences data table (output from \code{\link{osu_sequences}}) with
#' columns:
#'
#' \code{osu_id | species | qseqid | copy_number | pctsim | sequence}\cr
#'
#' @param dataset_names A character vector containing dataset names for each
#' element from the \code{result_list} list. If not given (NULL), each abundance
#' table (first element from each list) should have a 'dataset' column.
#'
#' @param verbose TRUE/FALSE: Print progress during dataset pooling.
#'
#' @export
pool_rexmap_results = function (result_list, dataset_names = NULL,
                               verbose=T) {
   # Add dataset name columns
   if (!is.null(dataset_names)){
      if (length(result_list) != length(dataset_names)) {
         stop('Error: The number of input data tables (', length(result_list),
              ') doesn\'t match the number of dataset names (',
              length(dataset_names), ').')
      } else {
         if (length(unique(dataset_names)) < 2) {
            stop('Error: Dataset names must be unique.')
         }
         for (i in 1:length(dataset_names)) {
            result_list[[i]][[1]]$dataset = dataset_names[i]
         }
      }
   }

   #                    v---internal function--v  v--list---v
   pooled_list = Reduce(pool_two_himap_outputs_2, result_list)

   return(pooled_list)
}

fremove = function (filename) {
   if (file.exists(filename)) {
      file.remove(filename)
   }
}

vcat = function (..., verbose=rexmap_option('verbose')) {
   if (verbose) {
      cat(...)
   }
}


strain_count = function (x, sep=',') {
   x_split = strsplit(x, sep, fixed=T)[[1]]
   x_split_expanded = lapply(x_split, function (x_i) {
      # Extract the number
      regex_count = sub('.*_\\[([0-9]+)\\]$', '\\1', x_i)
      if (regex_count == x_i) {
         # There is no count for this species/strain so just return as is
         return(x_i)
      } else {
         x_i_name = sub('^(.*)_\\[[0-9]+\\]$', '\\1', x_i)
         x_i_new = rep(x_i_name, as.integer(regex_count))
         x_i_new = paste(
            x_i_new,
            1:length(x_i_new),
            sep='_'
         )
         return(x_i_new)
      }
   })
   return(length(unlist(x_split_expanded)))
}

species_expand = function (x, sep=',', coll=T) {
   # Internal function for expanding collapsed strains from print_strains
   # Escherichia_coli_[3] to Escherichia_coli,Escherichia_coli,Escherichia_coli
   x_split = strsplit(x, sep, fixed=T)[[1]]
   x_split_expanded = lapply(x_split, function (x_i) {
      # Extract the number
      regex_count = sub('.*_\\[([0-9]+)\\]$', '\\1', x_i)
      if (regex_count == x_i) {
         # There is no count for this species/strain so just return as is
         return(x_i)
      } else {
         x_i_name = sub('^(.*)_\\[[0-9]+\\]$', '\\1', x_i)
         x_i_new = rep(x_i_name, as.integer(regex_count))
         x_i_new = paste(
            x_i_new,
            1:length(x_i_new),
            sep='_'
         )
         if (coll==TRUE) {
            return(paste(x_i_new, collapse=sep))
         } else {
            return(x_i_new)
         }

      }
   })
   if (coll==TRUE) {
      return(paste(x_split_expanded, collapse=sep))
   } else {
      return(unname(x_split_expanded))
   }
}

vector_split = function (x, sep=',') {
   strsplit(x, ',')[[1]]
}

species_relabel = function (x, sep=',') {
   print_strains(vector_split(species_expand(x, sep=sep)))
}

id_from_spectrum = function (v, sep1=',', sep2=':') {
   v_split = vector_split(v, sep=sep1)
   v_uniq_ids = sort(unique(sapply(
      v_split, function (x) strsplit(x, sep2)[[1]][1], USE.NAMES=F)))
   return(paste(v_uniq_ids, collapse=sep1))
}

spectrum_match = function (v1, v2, sep1=',', sep2=':',
                           mismatch=as.character(NA)) {
   # ASymmetric spectrum match: check if v1 is contained in any combination
   # of spectra from v2.
   # e.g. v1 = '1:5,140:1,209:1',
   # v2 = '1:5,1:2,1:1,140:1,209:1'
   # First number before colon is the identifier (seqid), second is copy
   # number.
   # we are checking to see if all numbers from v1 (1,140,209) occur exactly
   # once in v2, AND if there is a combination of these numbers and copy
   # numbers that match across v1 and v2.
   v1_split = vector_split(v1, sep=sep1)
   v2_split = vector_split(v2, sep=sep1)

   v1_uniq_ids = sort(unique(sapply(v1_split, function (x) strsplit(x, sep2)[[1]][1], USE.NAMES=F)))
   v2_uniq_ids = sort(unique(sapply(v2_split, function (x) strsplit(x, sep2)[[1]][1], USE.NAMES=F)))

   if (!identical(v1_uniq_ids, v2_uniq_ids)) {
      # one of the spectra has a unique ID that the other doesn't
      return(mismatch)
   } else {
      # OK, identical spectra, no check for matches from v1.
      if (all(sapply(v1_split, function (v) v %in% v2_split))) {
         # spectra match by copy number
         return(v1)
      } else {
         return(mismatch)
      }
   }
}


collapse_2 = function (ab_in, verbose=rexmap_option('verbose'),
                     ncpu=rexmap_option('ncpu'), temp_dir=tempdir(),
                     ws_scale=0.85, max_target_seqs=NULL,
                     blast_verbose=T) {
   # Provide ab_tab_nochim as an input argument, same as dada2::collapseNoMismatch
   if (verbose) cat('collapse:', fill=T)
   if (verbose) cat('* generating temporary files...')
   ab = copy(ab_in)
   out_files = ab_to_files(ab, verbose=blast_verbose, temp_dir=temp_dir)
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
   if (!is.null(blast_verbose)) {
      blast_verb = as.logical(blast_verbose)
   } else { # Blast verbose optionis not set (=NULL)
      blast_verb = verbose
   }
   blast_status = blastn(fasta, ref_db=db,
                         output=blast_out, max_target_seqs=max_target_seqs,
                         outfmt=paste0('6 ', rexmap_option('blast_coll_fmt')),
                         perc_identity=100, word_size=ws, ncpu=ncpu,
                         verbose=blast_verb)
   vcat('blast status: ', blast_status)
   if (blast_status != 0) {
      stop('\nError in BLAST alignment step.')
   }
   # if (verbose) cat('OK.\n')

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
   i = 0
   for (cl in cls) {
      i = i + 1
      if (length(cl) > 1) {
         # Always collapse to the shorter sequence
         # column_sums = colSums(ab[, cl, drop=FALSE])
         column_nchars = nchar(colnames(ab)[cl])
         # max_column is the sequence with max total number of reads
         # in the case of multiple just pick the first one in the list
         min_column = cl[which(column_nchars==min(column_nchars))][1]
         # cat('Pooling:', i, 'cl:', paste(cl, collapse=';'), 'min_col:',
         #     min_column, '\n')
         for (id in cl[cl != min_column]) {
            ab[1, min_column] = paste(ab[1, min_column],
                                     ab[1, id], sep=',')
            filtered_columns = c(filtered_columns, id)
         }
      }
   }
   ab = ab[, setdiff(1:ncol(ab), sort(filtered_columns)), drop=F]

   # ab_colsums = unname(colSums(ab))
   if (verbose) cat('OK.\n')

   # Cleanup the temp database
   if (verbose) cat('* cleaning up temporary files...')
   file.remove(blast_out)
   file.remove(fasta)
   cleanup_blastdb(db)
   if (verbose) cat('OK.\n')

   # Output abundance table sorted by total abundance (add up samples) in
   # descending order.
   return(ab)
}

collapse_dt = function (dt, seq_col='sequence', id_col='qseqid',
                        verbose=TRUE, ncpu=12) {

   # Generate fake abundance table so we can just use collapse
   # ab_tab = matrix(1, nrow=1, ncol=nrow(dt))
   # ab_tab = matrix(1, nrow=1, ncol=nrow(dt))
   # rownames(ab_tab) = 'qseqid'
   # colnames(ab_tab) = dt[, get(seq_col)]
   # ab_tab_coll2 = collapse(ab_in=ab_tab, verbose=verbose, ncpu=ncpu)
   # LOOKS GOOD. Numbers check out. TRY ON THE ENTIRE SET
   # We will need to pool some 99.6% OSUs to 100% OSUs here

   # DEBUG:
   # dt = unique(osu_seq_new.dt[, .(qseqid, sequence)])


   ab_tab = matrix(dt[, as.character(get(id_col))], nrow=1, ncol=nrow(dt))
   rownames(ab_tab) = 'qseqid'
   colnames(ab_tab) = dt[, get(seq_col)]
   ab_tab_coll = collapse_2(ab_in=ab_tab, verbose=verbose, ncpu=ncpu,
                            blast_verbose=verbose)
   # Now unpack the data
   dt_coll = data.table(qseqid_list=ab_tab_coll[1, ],
                        sequence=colnames(ab_tab_coll))
   dt_coll2 = dt_coll[, .(qseqid = as.integer(strsplit(qseqid_list, ',')[[1]])),
                      by=sequence]
   return(dt_coll2)
}



lu = function (x) length(unique(x))


pool_two_himap_outputs_2 = function (osu_1, osu_2,
                                   osu_offset=rexmap_option('osu_offset'),
                                   verbose=TRUE, temp_dir=tempdir(),
                                   blast_verbose=F, return_mapping=F) {


   # SET _1 = SUBJECT
   # SET _2 = QUERY
   pool_blast = function (qseq_1.dt, qseq_2.dt,
                          temp_dir=tempdir(), delete_temp=TRUE,
                          verbose=verbose, blast_verbose=F) {

      # Generate temp FASTA filenames
      rand_id = sample(LETTERS, 10)
      fasta_out_1 = file.path(
         temp_dir, paste(c('pool_', rand_id, '_1.fasta'), collapse='')
      )
      fasta_out_2 = file.path(
         temp_dir, paste(c('pool_', rand_id, '_2.fasta'), collapse='')
      )
      db_prefix_1 = file.path(dirname(fasta_out_1),
                              paste(c('pool_', rand_id, '_1'), collapse=''))

      # Generate temp BLAST DB
      fasta_writer(qseq_1.dt[, qseqid], qseq_1.dt[, sequence], fasta_out_1)
      fasta_writer(qseq_2.dt[, qseqid], qseq_2.dt[, sequence], fasta_out_2)

      # Generate a BLAST database from this fasta file
      makeblastdb(fasta_out_1, db_prefix_1, verbose=blast_verbose)

      # Automatically generate a good BLAST word size for set _1
      # Find the dataset that has shorter average sequence length
      nch_1 = qseq_1.dt[, mean(nchar(sequence))]
      nch_2 = qseq_2.dt[, mean(nchar(sequence))]
      shorter_set = 0
      if (nch_1 <= nch_2) {
         ws = max(round(min_seq_len(fasta_out_1) * 0.5), 7)
         shorter_set = 1
      } else {
         ws = max(round(min_seq_len(fasta_out_2) * 0.5), 7)
         shorter_set = 2
      }

      # Run BLAST and generate output table
      blast_out = file.path(
         temp_dir, paste(c('pool_', rand_id, '_blast_output.txt'), collapse='')
      )
      # QUERY IS SET _2 (qseqid = qseqid_2)
      # SUBJECT IS SET _1 (sseqid = qseqid_1)
      blast_status = blastn(fasta_out_2, ref_db=db_prefix_1,
                            output=blast_out, max_target_seqs=1000,
                            match=1, mismatch=-2, gapopen=1, gapextend=1,
                            outfmt=paste0('6 ', rexmap_option('blast_coll_fmt')),
                            perc_identity=100, word_size=ws, verbose=blast_verbose)
      vcat('blast status: ', blast_status, fill=T)

      # Cleanup temporary filenames
      if (delete_temp) {
         fremove(fasta_out_1)
         fremove(fasta_out_2)
         fremove(paste0(db_prefix_1, '.nhr'))
         fremove(paste0(db_prefix_1, '.nin'))
         fremove(paste0(db_prefix_1, '.nsq'))
      }

      # Load and parse BLAST results
      if (blast_status != 0) {
         stop('\nError in BLAST alignment step.')
      }
      vcat('OK.\n')

      # Load blast output and remove the temp blast output file
      vcat('* selecting ends-free alignments')
      blast.dt = data.table::fread(blast_out)
      vcat('.')
      names(blast.dt) = strsplit(rexmap_option('blast_coll_fmt'), ' ', fixed=T)[[1]]

      # Remove last temporary file (blast output)
      if (delete_temp) {
         fremove(blast_out)
      }


      # First pool out qseqids that do not show up in the alignment from both sets
      # unaligned_qseqids_1 = setdiff(qseq_1.dt[, qseqid], blast.dt[, qseqid])
      # unaligned_qseqids_2 = setdiff(qseq_2.dt[, qseqid], blast.dt[, sseqid])

      # For the pairs of qseqids that do match, keep only the longest alignments
      #### blast2.dt = blast.dt[, .SD[length==max(length)], by=.(qseqid)]
      blast2.dt = copy(blast.dt)
      vcat('.')
      # Filter out partial alignments !!!
      #
      #   Keep:
      #
      #   ------------------------------------------------------ 1 (q)
      #   ------------------------------------------------------ 2 (s)
      #
      #
      #      -----------------------   qstart=1, qend=qlen (fully aligned 1st)
      #     -------------------------
      #
      #     -------------------------
      #        --------------------     sstart=1, send=slen (fully aligned 2nd)
      #
      #       --------------------------- qstart=1, send=slen
      #    -------------------------
      #
      #   -----------------------------
      #      -----------------------------  qend=qlen, sstart=1
      #
      #  ------------------------------
      #  ---------------------------    sstart=1, send=slen
      #
      # Select only ends-free alignments
      blast3.dt = blast2.dt[
         (qstart==1 & qend==qlen) | (sstart==1 & send==slen) |
            (qstart==1 & send==slen) | (sstart==1 & qend==qlen)
         ]
      vcat('.')
      if (shorter_set == 2) {
         # Set _2 is shorter (qseqids in this table)
         blast3.dt[, qseqid_new := .GRP, by=qseqid]
      } else {
         # Set _1 is shorter (sseqids) so collapse to that one
         blast3.dt[, qseqid_new := .GRP, by=sseqid]
      }

      vcat(' OK.\n')


      # if (shorter_set == 2) {
      #    # _1 = sseqid (in debug this is SKOREA, longer)
      #    # _2 = qseqid (in debug example this is TwinsUK, shorter)
      #    # If longer sequence maps to multiple shorter, then we have a situation
      #    # like this:
      #    #                          -------------  (shorter)
      #    #  --------------------------------       (longer)
      #    # In these cases, let's not combine the reads, because then everything
      #    # else will have to be shortened to:
      #    #                          --------
      #    blast4.dt = blast3.dt[, if (length(unique(qseqid)) == 1) .SD, by=sseqid]
      #    # Now sort out the OSUs
      # } else {
      #    # Shorter_seq == 1
      #    # this is just like prev if block, but swap qseqid <--> sseqid
      #    blast4.dt = blast3.dt[, if (length(unique(sseqid)) == 1) .SD, by=qseqid]
      # }
      blast4.dt = copy(blast3.dt)
      return(list(blast4.dt, shorter_set, blast_status))
   }

   vcat('\n')
   vcat('------------------  RExMap POOL TWO DATASETS ----------------- \n')

   osu_ab_1.dt = osu_1[[1]]
   osu_seq_1.dt = na.omit(osu_1[[2]])
   osu_ab_2.dt = osu_2[[1]]
   osu_seq_2.dt = na.omit(osu_2[[2]])


   osu_seq_1.dt[, sequence := as.character(sequence)]
   osu_seq_2.dt[, sequence := as.character(sequence)]

   osu_ab_1_colorder = names(osu_ab_1.dt)
   osu_ab_2_colorder = names(osu_ab_2.dt)
   osu_seq_1_colorder = names(osu_seq_1.dt)
   osu_seq_2_colorder = names(osu_seq_2.dt)

   if (any(osu_ab_1_colorder != osu_ab_2_colorder)) {
      stop('Error: Abundance tables have mismatched column names/order.')
   }

   if (any(osu_seq_1_colorder != osu_seq_2_colorder)) {
      stop('Error: Sequence tables have mismatched column names/order.')
   }

   # Extract simple stats
   vcat('- DataSet 1: ', paste(osu_ab_1.dt[, unique(dataset)], collapse=','),
        paste0('(', osu_ab_1.dt[, lu(dataset)], ')'),
        '\n')
   vcat('             ')
   num_samples_1 = osu_ab_1.dt[, lu(sample_id)]
   vcat(' [# samples: ', num_samples_1, ']', sep='')
   num_osus_1 = osu_seq_1.dt[, lu(osu_id)]
   vcat(' [# OSUs: ', num_osus_1, ']', sep='')
   num_seqs_1 = osu_seq_1.dt[, lu(qseqid)]
   vcat(' [# seqs: ', num_seqs_1, ']', sep='')
   avg_read_len_1 = osu_seq_1.dt[, mean(nchar(sequence), na.rm=T)]
   vcat(' [avg seq len: ', avg_read_len_1, ']\n', sep='')

   # vcat('\n')
   vcat('- DataSet 2: ', paste(osu_ab_2.dt[, unique(dataset)], collapse=','),
        paste0('(', osu_ab_2.dt[, lu(dataset)], ')'),
        '\n')
   vcat('             ')
   num_samples_2 = osu_ab_2.dt[, lu(sample_id)]
   vcat(' [# samples: ', num_samples_2, ']', sep='')
   num_osus_2 = osu_seq_2.dt[, lu(osu_id)]
   vcat(' [# OSUs: ', num_osus_2, ']', sep='')
   num_seqs_2 = osu_seq_2.dt[, lu(qseqid)]
   vcat(' [# seqs: ', num_seqs_2, ']', sep='')
   avg_read_len_2 = osu_seq_2.dt[, mean(nchar(sequence), na.rm=T)]
   vcat(' [avg seq len: ', avg_read_len_2, ']\n', sep='')



   # To update the counts of all species, we will have to expand things like
   # Escherichia_coli_[3] to Escherichia_coli,Escherichia_coli,Escherichia_coli

   # ----------------- Extract sequences from Sets _1 and _2 -------------------
   vcat('* extracting sequences and sequence ids...')
   qseq_1.dt = na.omit(unique(osu_seq_1.dt[, .(qseqid, sequence)])[order(qseqid)])
   qseq_2.dt = na.omit(unique(osu_seq_2.dt[, .(qseqid, sequence)])[order(qseqid)])
   vcat(' OK.\n')

   vcat('* running blast alignment...')
   pool_blast_result = pool_blast(qseq_1.dt, qseq_2.dt, temp_dir=temp_dir,
                                  verbose=verbose, delete_temp=F)
   blast4.dt = pool_blast_result[[1]]
   shorter_set = pool_blast_result[[2]]
   blast_status = pool_blast_result[[3]]
   if (blast_status == 0) {
      # vcat(' OK.\n')
   } else {
      stop('\nError in blast alignment.\n')
   }
   rm(pool_blast_result)


   # Cleanup tables and prepare for merging with osu_seq tables
   names(blast4.dt)[1:2] = c('qseqid_2', 'qseqid_1')
   names(blast4.dt)[3:4] = c('q2len', 'q1len')
   names(blast4.dt)[6:9] = c('q2start', 'q2end', 'q1start', 'q1end')

   vcat('* generating mapping tables..')
   osu_seq_map_1.dt = unique(osu_seq_1.dt[
      , .(osu_id_1=osu_id, qseqid_1=qseqid, copy_number_1=copy_number,
          sequence_1=sequence, species_1=species, pctsim_1=pctsim)])
   osu_seq_map_1.dt[, osu_id_1_seqs_n := lu(qseqid_1), by=osu_id_1]
   vcat('.')
   osu_seq_map_2.dt = unique(osu_seq_2.dt[
      , .(osu_id_2=osu_id, qseqid_2=qseqid, copy_number_2=copy_number,
          sequence_2=sequence, species_2=species, pctsim_2=pctsim)])
   osu_seq_map_2.dt[, osu_id_2_seqs_n := lu(qseqid_2), by=osu_id_2]
   vcat(' OK.\n')

   #--------------------------- First pool 100% hits----------------------------
   # in blast4.dt, qseqid is qseqid from Set _2 and sseqid is qseqid from Set _1
   vcat('* matching 100% hits')
   osu_seq_maps_12.dt = merge(
      blast4.dt, osu_seq_map_1.dt[osu_id_1 < osu_offset], by='qseqid_1',
      allow.cartesian=T
   )
   osu_seq_maps_12.dt = merge(
      osu_seq_maps_12.dt, osu_seq_map_2.dt[osu_id_2 < osu_offset],
      by='qseqid_2', allow.cartesian=T
   )
   osu_seq_maps_12.dt[nchar(sequence_1) >= nchar(sequence_2), longer_seq := 1,
                      by=.(sequence_1, sequence_2)]
   osu_seq_maps_12.dt[nchar(sequence_1) < nchar(sequence_2), longer_seq := 2,
                      by=.(sequence_1, sequence_2)]
   vcat('.')

   # osu_seq_maps_12.dt = na.omit(osu_seq_maps_12.dt)
   osu_seq_maps_12_2.dt = copy(osu_seq_maps_12.dt)

   # Iterate over OSUs from shorter sequence and use that as a final reference
   # This works great 99.99% of the cases. In rare case V4 extends past
   # V3-V4 region by ~1 nt:
   #  (short seq)             -----------------------------------
   #  (long seq) -----------------------------------------------
   # We will take short seq as the final sequence, INCLUDING the extra nt.
   # So we don't recalculate pctsim and species.

   # Sometimes when the length difference is small, longer sequence might
   # need to be taken because it has less information, i.e.
   #
   # TACGGA...TGA
   # CACGGA...TGA
   #  ACGGA...TGACCGT
   #
   # In these cases we need to take, from shorter sequences
   #  - min(pctsim)
   #  - substring the sequence
   #  - max(osu_id)
   #  - species relabelled from the pooled species vector

   osu_id_new_offset = 0
   osu_seq_maps_12_2.dt[
      longer_seq==2,
      c('osu_id_new', 'sequence_new', 'pctsim_new', 'species_new',
        'copy_number_new') := list(
         .GRP + osu_id_new_offset,
         sequence_1,
         pctsim_1,
         species_1,
         copy_number_1
      ),
      by=osu_id_1
   ]
   vcat('.')
   if (nrow(na.omit(osu_seq_maps_12_2.dt)) == 0) {
      osu_id_new_offset = 0
   } else {
      osu_id_new_offset = osu_seq_maps_12_2.dt[, max(osu_id_new, na.rm=T)]
   }

   osu_seq_maps_12_2.dt[
      longer_seq==1,
      c('osu_id_new', 'sequence_new', 'pctsim_new', 'species_new',
        'copy_number_new') := list(
           .GRP + osu_id_new_offset,
           sequence_2,
           pctsim_2,
           species_2,
           copy_number_2
        ),
      by=osu_id_2
      ]

   vcat('.')
   osu_seq_maps_12_2.dt[
      ,
      osu_1_spectrum := unique(.SD[, .(qseqid_new, copy_number_new)])[
         , paste(sort(paste(qseqid_new, copy_number_new, sep=':')), collapse=',')],
      by=osu_id_1
      ]

   vcat('.')
   osu_seq_maps_12_2.dt[
      ,
      osu_2_spectrum := unique(.SD[, .(qseqid_new, copy_number_new)])[
         , paste(sort(paste(qseqid_new, copy_number_new, sep=':')), collapse=',')],
      by=osu_id_2
      ]

   vcat('.')
   if (osu_seq_maps_12_2.dt[longer_seq==2, .N] > 0) {
      osu_seq_maps_12_2.dt[
         longer_seq==2,
         spectrum_matched := spectrum_match(.BY[[1]], .BY[[2]]),
         by=.(osu_1_spectrum, osu_2_spectrum)
         ]
   }
   if (osu_seq_maps_12_2.dt[longer_seq==1, .N] > 0) {
      osu_seq_maps_12_2.dt[
         longer_seq==1,
         spectrum_matched := spectrum_match(.BY[[2]], .BY[[1]]),
         by=.(osu_1_spectrum, osu_2_spectrum)
         ]
   }
   vcat('.')
   osu_seq_maps_12_2.dt = osu_seq_maps_12_2.dt[!is.na(spectrum_matched)]

   # If we have multiple spectra matched and OSUs from SHORTER set that map
   # to single OSUs in the LONGER SET (which only happens due to optimization)
   # take the OSU match that has most matched species!
   # We had E. coli strains filling in the counts in my tests.

   osu_seq_maps_12_2.dt[, species_1_count := strain_count(.BY[[1]]), by=species_1]
   osu_seq_maps_12_2.dt[, species_2_count := strain_count(.BY[[1]]), by=species_2]

   osu_seq_maps_12_2a.dt = osu_seq_maps_12_2.dt[
      longer_seq==2,
      if (lu(osu_id_1) > 1) {
         .SD[species_1_count == max(species_1_count)][1] } else { .SD },
      by=osu_id_2
      ]

   osu_seq_maps_12_2b.dt = osu_seq_maps_12_2.dt[
      longer_seq==1,
      if (lu(osu_id_2) > 1) { .SD[species_2_count == max(species_2_count)][1] } else { .SD },
      by=osu_id_1
      ]


   osu_seq_maps_12_3.dt = merge(osu_seq_maps_12_2a.dt, osu_seq_maps_12_2b.dt,
                                by=names(osu_seq_maps_12_2a.dt), all=T)

   vcat(' OK.\n')

   #--------------------------------- Pool  < 100% hits ! ----------------------

   vcat('* matching non 100% hits')
   osu_seq_maps_12_non100.dt = merge(
      blast4.dt, osu_seq_map_1.dt[osu_id_1 >= osu_offset], by='qseqid_1'
   )
   osu_seq_maps_12_non100.dt = merge(
      osu_seq_maps_12_non100.dt, osu_seq_map_2.dt[osu_id_2 >= osu_offset],
      by='qseqid_2'
   )

   osu_seq_maps_12_non100.dt[nchar(sequence_1) >= nchar(sequence_2),
                             longer_seq := 1,
                             by=.(sequence_1, sequence_2)]
   osu_seq_maps_12_non100.dt[nchar(sequence_1) < nchar(sequence_2),
                             longer_seq := 2,
                             by=.(sequence_1, sequence_2)]

   vcat('.')
   # Again iterate over shorter seq and pool the rest
   osu_id_new_offset = osu_offset
   osu_seq_maps_12_non100.dt[
      longer_seq==2,
      c('osu_id_new', 'sequence_new', 'pctsim_new', 'species_new',
        'copy_number_new') := list(
           .GRP + osu_id_new_offset,
           sequence_1,
           pctsim_1,
           species_1,
           copy_number_1
        ),
      by=osu_id_1
   ]

   osu_seq_maps_12_non100_tofix.dt = osu_seq_maps_12_non100.dt[
      longer_seq==2,
      if (lu(osu_id_1) > 1) .(
         osu_id_1 = osu_id_1,
         osu_id_new = osu_id_new,
         osu_id_new_fix = min(osu_id_new),
         pctsim_new_fix = min(pctsim_1),
         species_new_fix = species_relabel(paste(species_1, collapse=',')),
         sequence_new_fix = substring(sequence_1[1], q1start, q1end)
      ), by=osu_id_2]

   if (nrow(osu_seq_maps_12_non100_tofix.dt) > 0) {
      osu_seq_maps_12_non100.dt = merge(
         osu_seq_maps_12_non100.dt, osu_seq_maps_12_non100_tofix.dt, all.x=T,
         by=c('osu_id_2', 'osu_id_new', 'osu_id_1')
      )
      osu_seq_maps_12_non100.dt[
         !is.na(osu_id_new_fix)
         , c('osu_id_new', 'species_new', 'pctsim_new', 'sequence_new') := list(
            as.integer(osu_id_new_fix), species_new_fix, pctsim_new_fix,
            sequence_new_fix
         )
      ]
      osu_seq_maps_12_non100.dt[
         , c('osu_id_new_fix', 'species_new_fix', 'pctsim_new_fix',
             'sequence_new_fix') := NULL
      ]
   }

   vcat('.')
   if (nrow(na.omit(osu_seq_maps_12_non100.dt)) == 0) {
      osu_id_new_offset = osu_offset
   } else {
      osu_id_new_offset = osu_seq_maps_12_non100.dt[, max(osu_id_new, na.rm=T)]
   }

   osu_seq_maps_12_non100.dt[
      longer_seq==1,
      c('osu_id_new', 'sequence_new', 'pctsim_new', 'species_new',
        'copy_number_new') := list(
           .GRP + osu_id_new_offset,
           sequence_2,
           pctsim_2,
           species_2,
           copy_number_2
        ),
      by=osu_id_2
   ]

   osu_seq_maps_12_non100_tofix2.dt = osu_seq_maps_12_non100.dt[
      longer_seq==1,
      if (lu(osu_id_2) > 1) .(
         osu_id_2 = osu_id_2,
         osu_id_new = osu_id_new,
         osu_id_new_fix = min(osu_id_new),
         pctsim_new_fix = min(pctsim_2),
         species_new_fix = species_relabel(paste(species_2, collapse=',')),
         sequence_new_fix = substring(sequence_2[1], q2start, q2end)
      ), by=osu_id_1]

   if (nrow(osu_seq_maps_12_non100_tofix2.dt) > 0) {
      osu_seq_maps_12_non100.dt = merge(
         osu_seq_maps_12_non100.dt, osu_seq_maps_12_non100_tofix2.dt, all.x=T,
         by=c('osu_id_1', 'osu_id_new', 'osu_id_2')
      )
      osu_seq_maps_12_non100.dt[
         !is.na(osu_id_new_fix)
         , c('osu_id_new', 'species_new', 'pctsim_new', 'sequence_new') := list(
            as.integer(osu_id_new_fix), species_new_fix, pctsim_new_fix,
            sequence_new_fix
         )
         ]
      osu_seq_maps_12_non100.dt[
         , c('osu_id_new_fix', 'species_new_fix', 'pctsim_new_fix',
             'sequence_new_fix') := NULL
         ]
   }

   vcat(' OK.\n')

   vcat('* renumerating unmatched OSUs  [Set 1 100% hits')
   #-------- Pool 100 % hit unmatched from Set _1 ------------------------------
   osu_id_new_offset = osu_seq_maps_12_3.dt[, max(osu_id_new, na.rm=T)]
   qseqid_new_offset = max(
      osu_seq_maps_12_3.dt[, max(qseqid_new, na.rm=T)],
      osu_seq_maps_12_non100.dt[, max(qseqid_new, na.rm=T)]
   )
   osu_seq_1_nonmatched.dt = osu_seq_1.dt[osu_id < osu_offset &
      !(osu_id %in% osu_seq_maps_12_3.dt[, unique(osu_id_1)])
   ]
   osu_seq_1_nonmatched.dt[, osu_id_new := .GRP + osu_id_new_offset, by=osu_id]
   osu_seq_1_nonmatched.dt[, qseqid_new := .GRP + qseqid_new_offset, by=qseqid]
   # Check if any of these qseqids were already assigned a new one before in the
   # matched table.
   osu_seq_1_nonmatched.dt = merge(
      osu_seq_1_nonmatched.dt,
      unique(osu_seq_maps_12_3.dt[, .(qseqid_prev=qseqid_new, sequence=sequence_new)]),
      by='sequence', all.x=T
   )
   osu_seq_1_nonmatched.dt[!is.na(qseqid_prev), qseqid_new := qseqid_prev]
   osu_seq_1_nonmatched.dt[, qseqid_prev := NULL]
   vcat(' OK] [Set 2 100% hits')

   #------------- Pool 100 % hits unmatched from Set _2 ------------------------
   osu_id_new_offset = osu_seq_1_nonmatched.dt[, max(osu_id_new, na.rm=T)]
   qseqid_new_offset = osu_seq_1_nonmatched.dt[, max(qseqid_new, na.rm=T)]

   osu_seq_2_nonmatched.dt = osu_seq_2.dt[osu_id < osu_offset &
      !(osu_id %in% osu_seq_maps_12_3.dt[, unique(osu_id_2)])
   ]
   osu_seq_2_nonmatched.dt[, osu_id_new := .GRP + osu_id_new_offset, by=osu_id]
   osu_seq_2_nonmatched.dt[, qseqid_new := .GRP + qseqid_new_offset, by=qseqid]
   vcat(' OK] [Set 1 < 100% hits')

   #------------- < 100% unmatched from set _1 ---------------------------------
   # if (nrow(na.omit(osu_seq_maps_12_non100.dt)) > 0) {
   #
   # }
   osu_id_new_offset = osu_seq_maps_12_non100.dt[, max(osu_id_new, na.rm=T)]
   qseqid_new_offset = osu_seq_2_nonmatched.dt[, max(qseqid_new, na.rm=T)]
   osu_seq_1_nonmatched_non100.dt = osu_seq_1.dt[
      osu_id >= osu_offset &
         !(osu_id %in% osu_seq_maps_12_non100.dt[, unique(osu_id_1)])
      ]
   osu_seq_1_nonmatched_non100.dt[, osu_id_new := .GRP + osu_id_new_offset,
                                  by=osu_id]
   osu_seq_1_nonmatched_non100.dt[, qseqid_new := .GRP + qseqid_new_offset,
                                  by=qseqid]
   vcat(' OK] [Set 2 < 100% hits')

   #------------- < 100% unmatched from set _2 ---------------------------------
   osu_id_new_offset = max(
      c(
         # osu_seq_1_nonmatched_non100.dt[, osu_id_new],
         osu_seq_1_nonmatched.dt[, osu_id_new],
         osu_offset
      ), na.rm=T)
   # qseqid_new_offset = osu_seq_1_nonmatched_non100.dt[
   #    , max(qseqid_new, na.rm=T)]
   qseqid_new_offset = max(
      c(
         osu_seq_1_nonmatched_non100.dt[, qseqid_new],
         osu_seq_1_nonmatched.dt[, qseqid_new]
      )
   )
   osu_seq_2_nonmatched_non100.dt = osu_seq_2.dt[
      osu_id >= osu_offset &
         !(osu_id %in% osu_seq_maps_12_non100.dt[, unique(osu_id_2)])
      ]
   osu_seq_2_nonmatched_non100.dt[, osu_id_new := .GRP + osu_id_new_offset,
                                  by=osu_id]
   osu_seq_2_nonmatched_non100.dt[, qseqid_new := .GRP + qseqid_new_offset,
                                  by=qseqid]
   vcat(' OK]\n')

   #---------------- Pool everything together into a new data table ------------
   vcat('* merging mapping tables..')

   osu_seq_remap_1.dt = data.table::rbindlist(list(
      # 100% that match something from _2
      unique(
         osu_seq_maps_12_3.dt[
            , .(osu_id_1, osu_id_new, species_new, pctsim_new, sequence_new,
                qseqid_new, copy_number_new)
            ]
      ),
      # 100% hits that do not match anything in _2
      unique(
         osu_seq_1_nonmatched.dt[
            !is.na(pctsim)
            , .(osu_id_1=osu_id, osu_id_new, species_new=species,
                pctsim_new=pctsim, sequence_new=sequence,
                qseqid_new, copy_number_new=copy_number)]
      ),
      # non 100% htis that match in _2
      unique(
         osu_seq_maps_12_non100.dt[
            , .(osu_id_1, osu_id_new, species_new, pctsim_new, sequence_new,
                qseqid_new, copy_number_new)
            ]
      ),
      # non 100% hits that don't match anything in _2
      unique(
         osu_seq_1_nonmatched_non100.dt[
            , .(osu_id_1=osu_id, osu_id_new, species_new=species,
                pctsim_new=pctsim, sequence_new=sequence,
                qseqid_new, copy_number_new=copy_number)]
      )
   ))
   vcat('.')

   osu_seq_remap_2.dt = data.table::rbindlist(list(
      # 100% that match something from _2
      unique(
         osu_seq_maps_12_3.dt[
            , .(osu_id_2, osu_id_new, species_new, pctsim_new, sequence_new,
                qseqid_new, copy_number_new)
            ]
      ),
      # 100% hits that do not match anything in _2
      unique(
         osu_seq_2_nonmatched.dt[
            !is.na(pctsim)
            , .(osu_id_2=osu_id, osu_id_new, species_new=species,
                pctsim_new=pctsim, sequence_new=sequence,
                qseqid_new, copy_number_new=copy_number)]
      ),
      # non 100% htis that match in _2
      unique(
         osu_seq_maps_12_non100.dt[
            , .(osu_id_2, osu_id_new, species_new, pctsim_new, sequence_new,
                qseqid_new, copy_number_new)
            ]
      ),
      # non 100% hits that don't match anything in _2
      unique(
         osu_seq_2_nonmatched_non100.dt[
            , .(osu_id_2=osu_id, osu_id_new, species_new=species,
                pctsim_new=pctsim, sequence_new=sequence,
                qseqid_new, copy_number_new=copy_number)]
      )
   ))
   vcat(' OK.\n')

   # Generate a new osu_seq table for _1
   vcat('* generating new sequence table for set _1 ...')
   osu_seq_1_new.dt = merge(
      osu_seq_1.dt, osu_seq_remap_1.dt, by.x='osu_id', by.y='osu_id_1'
   )
   osu_seq_1_new.dt[, c('sequence', 'species', 'qseqid', 'pctsim',
                        'copy_number', 'osu_id') := NULL]
   osu_seq_1_new.dt = unique(osu_seq_1_new.dt)
   vcat(' OK.\n')

   # Generate a new osu_seq table for _2
   vcat('* generating new sequence table for set _2 ...')
   osu_seq_2_new.dt = merge(
      osu_seq_2.dt, osu_seq_remap_2.dt, by.x='osu_id', by.y='osu_id_2'
   )
   osu_seq_2_new.dt[, c('sequence', 'species', 'qseqid', 'pctsim',
                        'copy_number', 'osu_id') := NULL]
   osu_seq_2_new.dt = unique(osu_seq_2_new.dt)
   vcat(' OK.\n')

   # Now use these mapping files to add new osu_id, species and pctsim
   # to osu_ab and osu_seq tables.
   # Then remove old OSU_IDs
   # For old OSU_AB add up osu_counts, then use unique to finalize the new
   # table.

   # Generate a new osu_ab table for set _1
   vcat('* generating new abundance tables for set _1')
   osu_ab_1_new.dt = merge(
      osu_ab_1.dt,
      unique(osu_seq_remap_1.dt[
         , .(osu_id=osu_id_1, osu_id_new, species_new, pctsim_new)]
      ),
      by='osu_id'
   )
   vcat('.')
   osu_ab_1_new.dt[, c('osu_id', 'pctsim', 'species') := NULL]
   vcat('.')
   names(osu_ab_1_new.dt)[names(osu_ab_1_new.dt) == 'osu_id_new'] = 'osu_id'
   names(osu_ab_1_new.dt)[names(osu_ab_1_new.dt) == 'species_new'] = 'species'
   names(osu_ab_1_new.dt)[names(osu_ab_1_new.dt) == 'pctsim_new'] = 'pctsim'
   setcolorder(osu_ab_1_new.dt, osu_ab_1_colorder)

   names(osu_seq_1_new.dt)[names(osu_seq_1_new.dt) == 'osu_id_new'] = 'osu_id'
   names(osu_seq_1_new.dt)[names(osu_seq_1_new.dt) == 'species_new'] = 'species'
   names(osu_seq_1_new.dt)[names(osu_seq_1_new.dt) == 'pctsim_new'] = 'pctsim'
   names(osu_seq_1_new.dt)[names(osu_seq_1_new.dt) == 'qseqid_new'] = 'qseqid'
   names(osu_seq_1_new.dt)[names(osu_seq_1_new.dt) == 'sequence_new'] = 'sequence'
   names(osu_seq_1_new.dt)[names(osu_seq_1_new.dt) == 'copy_number_new'] =
      'copy_number'
   setcolorder(osu_seq_1_new.dt, osu_seq_1_colorder)

   # Pool counts
   osu_ab_1_new.dt[, osu_count := sum(osu_count), by=.(osu_id, sample_id, dataset)]
   vcat('.')
   osu_ab_1_new.dt = unique(osu_ab_1_new.dt)
   vcat(' OK.\n')

   vcat('* generating new abundance tables for set _2')
   osu_ab_2_new.dt = merge(
      osu_ab_2.dt,
      unique(osu_seq_remap_2.dt[
         , .(osu_id=osu_id_2, osu_id_new, species_new, pctsim_new)]
      ),
      by='osu_id'
   )
   vcat('.')
   osu_ab_2_new.dt[, c('osu_id', 'pctsim', 'species') := NULL]
   vcat('.')
   names(osu_ab_2_new.dt)[names(osu_ab_2_new.dt) == 'osu_id_new'] = 'osu_id'
   names(osu_ab_2_new.dt)[names(osu_ab_2_new.dt) == 'species_new'] = 'species'
   names(osu_ab_2_new.dt)[names(osu_ab_2_new.dt) == 'pctsim_new'] = 'pctsim'
   setcolorder(osu_ab_2_new.dt, osu_ab_2_colorder)

   names(osu_seq_2_new.dt)[names(osu_seq_2_new.dt) == 'osu_id_new'] = 'osu_id'
   names(osu_seq_2_new.dt)[names(osu_seq_2_new.dt) == 'species_new'] = 'species'
   names(osu_seq_2_new.dt)[names(osu_seq_2_new.dt) == 'pctsim_new'] = 'pctsim'
   names(osu_seq_2_new.dt)[names(osu_seq_2_new.dt) == 'qseqid_new'] = 'qseqid'
   names(osu_seq_2_new.dt)[names(osu_seq_2_new.dt) == 'sequence_new'] = 'sequence'
   names(osu_seq_2_new.dt)[names(osu_seq_2_new.dt) == 'copy_number_new'] =
      'copy_number'
   setcolorder(osu_seq_2_new.dt, osu_seq_2_colorder)

   # Pool counts
   osu_ab_2_new.dt[, osu_count := sum(osu_count), by=.(osu_id, sample_id, dataset)]
   vcat('.')
   osu_ab_2_new.dt = unique(osu_ab_2_new.dt)
   vcat(' OK.\n')

   vcat('* final merging full spectrum matches...')
   osu_ab_new.dt = data.table::rbindlist(list(osu_ab_1_new.dt, osu_ab_2_new.dt))
   vcat('.')
   osu_seq_new.dt = unique(data.table::rbindlist(list(osu_seq_1_new.dt, osu_seq_2_new.dt)))
   vcat(' OK.\n')
   osu_seq_new.dt[, qseqid := .GRP, by=sequence]


   #------------------------ Final collapse ------------------------------------
   # As a final assignment just add up all these strain names, and pick lowest
   # OSU as a final OSU
   # OSU_SEQ CORRECTIOn
   # Final collapse, some OSUs that have partial matches have not yet been updated
   # so this takes care of those cases.
   vcat('* final collapse sequence table')
   qseq_seq.dt = collapse_dt(unique(osu_seq_new.dt[!is.na(sequence), .(qseqid, sequence)]),
                             verbose=blast_verbose)
   osu_seq_new.dt = merge(osu_seq_new.dt,
                          qseq_seq.dt[, .(qseqid, seq_new=sequence)],
                          by='qseqid', all.x=T)
   osu_seq_new.dt[sequence != seq_new, sequence := seq_new]
   osu_seq_new.dt[, seq_new := NULL]
   setcolorder(osu_seq_new.dt, osu_seq_1_colorder)

   # Renumerate qseqids again after this collapse
   osu_seq_new.dt[, qseqid := .GRP, by=sequence]

   # Now pool OSUs with same qseqid spectrum
   osu_seq_new.dt[, spectrum := paste(qseqid, copy_number, sep=':', collapse=','),
                  by=osu_id]

   # Check for OSUs that have only 1 qseqid in their spectrum, but different
   # copy numbers.
   # By def. !(spectrum %like% ',') we will get OSUs with only these qseqids
   # osu_seq_new_singles.dt = osu_seq_new.dt[!(spectrum %like% ',') & pctsim==100,
   #                if (lu(osu_id) > 1) .SD, by=qseqid]
   #
   # osu_seq_new_singles.dt[, ]

   vcat('')

   # osu_seq_new.dt[, c('osu_id_new', 'species_new', 'pctsim_new'):= list(
   #    min(osu_id),
   #    .SD[pctsim==max(pctsim), species_relabel(species)],
   #    max(pctsim)
   # ), by=spectrum]

   # Unchanged OSUs
   osu_seq_new_unch.dt = osu_seq_new.dt[, if (lu(osu_id) == 1) .SD, by=spectrum]
   osu_seq_new_ch.dt = osu_seq_new.dt[, if (lu(osu_id) > 1) .SD, by=spectrum]

   osu_seq_new_ch.dt[, c('osu_id_new', 'species_new', 'pctsim_new') := list(
      .SD[pctsim==max(pctsim), min(osu_id)],
      .SD[pctsim==max(pctsim), species_relabel(paste(species, collapse=','))],
      max(pctsim)
   ), by=spectrum]
   vcat('.')
   # Generate OSU remapping table
   osu_osu_remap.dt = unique(osu_seq_new_ch.dt[
      , .(osu_id, osu_id_new, species_new, pctsim_new)])

   osu_seq_new_ch.dt[, c('osu_id', 'species', 'pctsim') := list(
      osu_id_new, species_new, pctsim_new
   )]

   osu_seq_new_ch.dt[, c('osu_id_new', 'species_new', 'pctsim_new') := NULL]

   osu_seq_new2.dt = data.table::rbindlist(list(osu_seq_new_unch.dt, osu_seq_new_ch.dt))
   osu_seq_new2.dt[, spectrum := NULL]
   setcolorder(osu_seq_new2.dt, osu_seq_1_colorder)
   vcat('. OK.\n')

   # OSU_AB CORRECTION
   vcat('* final collapse abundance table...')
   osu_ab_new.dt = merge(osu_ab_new.dt, osu_osu_remap.dt, by='osu_id', all.x=T)
   osu_ab_new.dt[
      !is.na(osu_id_new),
      c('osu_id', 'osu_count', 'pctsim', 'species') := list(
         osu_id_new,
         sum(osu_count),
         pctsim_new,
         species_new
      ), by=.(osu_id_new, sample_id, dataset)]
   osu_ab_new.dt[, c('osu_id_new', 'species_new', 'pctsim_new') := NULL]
   osu_ab_new.dt = unique(osu_ab_new.dt)
   setcolorder(osu_ab_new.dt, osu_ab_1_colorder)
   osu_ab_new.dt[, pctsim := min(pctsim), by=.(sample_id, dataset, osu_id)]
   osu_ab_new.dt = unique(osu_ab_new.dt)
   vcat(' OK.\n')

   # if (nrow(osu_seq_new.dt[, if (lu(qseqid) > 1) .SD, by=sequence]) > 0) {
   #    vcat('WARNING: some sequences do not map to single qseqids!')
   # }

   vcat('Done.')

   vcat('\n')

   vcat('- Summary\n')
   vcat('* Input: DataSet 1: ', num_osus_1, ' OSUs. DataSet 2: ', num_osus_2,
        ' OSUs.\n', sep='')
   vcat('         ', num_osus_1+num_osus_2, ' total OSUs.\n', sep='')
   vcat('* 100% OSUs pooled: ', osu_seq_maps_12_3.dt[, lu(osu_id_1)],
        ' Set1 OSUs & ', osu_seq_maps_12_3.dt[, lu(osu_id_2)],
        ' Set2 OSUs => ', osu_seq_maps_12_3.dt[, lu(osu_id_new)], ' OSUs.\n',
        sep='')
   vcat('* < 100% OSUs pooled: ', osu_seq_maps_12_non100.dt[, lu(osu_id_1)],
        ' Set1 OSUs & ', osu_seq_maps_12_non100.dt[, lu(osu_id_2)],
        ' Set2 OSUs => ', osu_seq_maps_12_non100.dt[, lu(osu_id_new)], ' OSUs.\n',
        sep='')
   vcat('* Output: ', osu_seq_new2.dt[, lu(osu_id)], ' final OSUs.\n', sep='')


   #------------------------- Return output ------------------------------------
   if (!return_mapping) {
      # Default output when called from the main pooling function
      return(list(
         osu_ab_new.dt,
         unique(osu_seq_new2.dt[order(osu_id)])
      ))

   } else {
      # Output for OSU ID1 vs OSU ID2 from two datasets
      return(list(
         osu_seq_maps_12_3.dt,
         osu_seq_maps_12_non100.dt
      ))
   }
}


#' Match OSUs across datasets
#'
#'
#' @param osu_1 A list of two data tables.
#'
#' 1. OSU abundance data table (output from \code{\link{abundance}}) with
#' columns:
#'
#' \code{sample_id | osu_id | osu_count | pctsim | species}\cr
#'
#' Additionally, this data table can contain an extra \code{dataset} column
#' which contains a name of that dataset.
#'
#' 2. OSU sequences data table (output from \code{\link{osu_sequences}}) with
#' columns:
#'
#' \code{osu_id | species | qseqid | copy_number | pctsim | sequence}\cr
#'
#' @param osu_2 Same as \code{osu_1} just for a second dataset.
#'
#' @export
osu_dataset_match = function (
   osu_1, osu_2, dataset_names = c('dataset1', 'dataset2'),
   osu_offset=rexmap_option('osu_offset'), verbose=T,
   temp_dir=tempdir(), blast_verbose=T) {

   osu_1[[1]]$dataset = dataset_names[1]
   osu_2[[1]]$dataset = dataset_names[2]
   # Run
   mappings.list = pool_two_himap_outputs_2(
      osu_1=osu_1, osu_2=osu_2, osu_offset=osu_offset, verbose=verbose,
      temp_dir=temp_dir, blast_verbose=blast_verbose, return_mapping=T
   )
   map.dt = data.table::rbindlist(list(
      unique(
         mappings.list[[1]][, .(osu_id_1, species_1, pctsim_1, osu_id_2, species_2, pctsim_2)]
      ),
      unique(
         mappings.list[[2]][, .(osu_id_1, species_1, pctsim_1, osu_id_2, species_2, pctsim_2)]
      )
   ))
   map.dt = map.dt[order(osu_id_1)]
   return(map.dt)
}

pool_two_himap_outputs = function (osu_1, osu_2,
                                   osu_offset=rexmap_option('osu_offset'),
                                   verbose=T, temp_dir=tempdir()) {

   # SET _1 = SUBJECT
   # SET _2 = QUERY
   pool_blast = function (qseq_1.dt, qseq_2.dt, temp_dir=tempdir(), verbose=T) {

      # Generate temp FASTA filenames
      rand_id = sample(LETTERS, 10)
      fasta_out_1 = file.path(
         temp_dir, paste(c('pool_', rand_id, '_1.fasta'), collapse='')
      )
      fasta_out_2 = file.path(
         temp_dir, paste(c('pool_', rand_id, '_2.fasta'), collapse='')
      )
      db_prefix_1 = file.path(dirname(fasta_out_1),
                              paste(c('pool_', rand_id, '_1'), collapse=''))

      # Generate temp BLAST DB
      fasta_writer(qseq_1.dt[, qseqid], qseq_1.dt[, sequence], fasta_out_1)
      fasta_writer(qseq_2.dt[, qseqid], qseq_2.dt[, sequence], fasta_out_2)

      # Generate a BLAST database from this fasta file
      makeblastdb(fasta_out_1, db_prefix_1, verbose=verbose)

      # Automatically generate a good BLAST word size for set _1
      # Find the dataset that has shorter average sequence length
      nch_1 = qseq_1.dt[, mean(nchar(sequence))]
      nch_2 = qseq_2.dt[, mean(nchar(sequence))]
      shorter_set = 0
      if (nch_1 <= nch_2) {
         ws = max(round(min_seq_len(fasta_out_1) * 0.5), 7)
         shorter_set = 1
      } else {
         ws = max(round(min_seq_len(fasta_out_2) * 0.5), 7)
         shorter_set = 2
      }

      # Run BLAST and generate output table
      blast_out = file.path(
         temp_dir, paste(c('pool_', rand_id, '_blast_output.txt'), collapse='')
      )
      # QUERY IS SET _2 (qseqid = qseqid_2)
      # SUBJECT IS SET _1 (sseqid = qseqid_1)
      blast_status = blastn(fasta_out_2, ref_db=db_prefix_1,
                            output=blast_out, max_target_seqs=NULL,
                            match=1, mismatch=-2, gapopen=1, gapextend=1,
                            outfmt=paste0('6 ', rexmap_option('blast_coll_fmt')),
                            perc_identity=100, word_size=ws, verbose=T)
      vcat('blast status: ', blast_status, fill=T)

      # Cleanup temporary filenames
      fremove(fasta_out_1)
      fremove(fasta_out_2)
      fremove(paste0(db_prefix_1, '.nhr'))
      fremove(paste0(db_prefix_1, '.nin'))
      fremove(paste0(db_prefix_1, '.nseq'))

      # Load and parse BLAST results
      if (blast_status != 0) {
         stop('\nError in BLAST alignment step.')
      }
      vcat('OK.\n')

      # Load blast output and remove the temp blast output file
      vcat('* selecting ends-free alignments')
      blast.dt = data.table::fread(blast_out)
      vcat('.')
      names(blast.dt) = strsplit(rexmap_option('blast_coll_fmt'), ' ', fixed=T)[[1]]

      # Remove last temporary file (blast output)
      fremove(blast_out)

      # First pool out qseqids that do not show up in the alignment from both sets
      # unaligned_qseqids_1 = setdiff(qseq_1.dt[, qseqid], blast.dt[, qseqid])
      # unaligned_qseqids_2 = setdiff(qseq_2.dt[, qseqid], blast.dt[, sseqid])

      # For the pairs of qseqids that do match, keep only the longest alignments
      #### blast2.dt = blast.dt[, .SD[length==max(length)], by=.(qseqid)]
      blast2.dt = copy(blast.dt)
      vcat('.')
      # Filter out partial alignments !!!
      #
      #   Keep:
      #
      #   ------------------------------------------------------ 1 (q)
      #   ------------------------------------------------------ 2 (s)
      #
      #
      #      -----------------------   qstart=1, qend=qlen (fully aligned 1st)
      #     -------------------------
      #
      #     -------------------------
      #        --------------------     sstart=1, send=slen (fully aligned 2nd)
      #
      #       --------------------------- qstart=1, send=slen
      #    -------------------------
      #
      #   -----------------------------
      #      -----------------------------  qend=qlen, sstart=1
      #
      #  ------------------------------
      #  ---------------------------    sstart=1, send=slen
      #
      # Select only ends-free alignments
      blast3.dt = blast2.dt[
         (qstart==1 & qend==qlen) | (sstart==1 & send==slen) |
            (qstart==1 & send==slen) | (sstart==1 & qend==qlen)
         ]
      vcat('.')
      if (shorter_set == 2) {
         # Set _2 is shorter (qseqids in this table)
         blast3.dt[, qseqid_new := .GRP, by=qseqid]
      } else {
         # Set _1 is shorter (sseqids) so collapse to that one
         blast3.dt[, qseqid_new := .GRP, by=sseqid]
      }

      vcat(' OK.\n')


      # if (shorter_set == 2) {
      #    # _1 = sseqid (in debug this is SKOREA, longer)
      #    # _2 = qseqid (in debug example this is TwinsUK, shorter)
      #    # If longer sequence maps to multiple shorter, then we have a situation
      #    # like this:
      #    #                          -------------  (shorter)
      #    #  --------------------------------       (longer)
      #    # In these cases, let's not combine the reads, because then everything
      #    # else will have to be shortened to:
      #    #                          --------
      #    blast4.dt = blast3.dt[, if (length(unique(qseqid)) == 1) .SD, by=sseqid]
      #    # Now sort out the OSUs
      # } else {
      #    # Shorter_seq == 1
      #    # this is just like prev if block, but swap qseqid <--> sseqid
      #    blast4.dt = blast3.dt[, if (length(unique(sseqid)) == 1) .SD, by=qseqid]
      # }
      blast4.dt = copy(blast3.dt)
      return(list(blast4.dt, shorter_set))
   }


   # osu_1 = list(osu_ab_1.dt, osu_seq_1.dt)
   vcat('\n')
   vcat('POOLING: \n')
   vcat('- Dataset 1: ', paste(unique(osu_1[[1]]$dataset), collapse=','), '\n')
   vcat('- Dataset 2: ', paste(unique(osu_2[[1]]$dataset), collapse=','), '\n')
   vcat('\n')

   osu_ab_1.dt = osu_1[[1]]
   osu_seq_1.dt = osu_1[[2]]
   osu_ab_2.dt = osu_2[[1]]
   osu_seq_2.dt = osu_2[[2]]


   osu_seq_1.dt[, sequence := as.character(sequence)]
   # osu_seq_1.dt[, nqs := lu(qseqid), by=osu_id]
   osu_seq_2.dt[, sequence := as.character(sequence)]
   # osu_seq_2.dt[, nqs := lu(qseqid), by=osu_id]

   osu_ab_1_colorder = names(osu_ab_1.dt)

   # To update the counts of all species, we will have to expand things like
   # Escherichia_coli_[3] to Escherichia_coli,Escherichia_coli,Escherichia_coli

   species_expand = function (x, sep=',', coll=T) {
      # Internal function for expanding collapsed strains from print_strains
      # Escherichia_coli_[3] to Escherichia_coli,Escherichia_coli,Escherichia_coli
      x_split = strsplit(x, sep, fixed=T)[[1]]
      x_split_expanded = lapply(x_split, function (x_i) {
         # Extract the number
         regex_count = sub('.*_\\[([0-9]+)\\]$', '\\1', x_i)
         if (regex_count == x_i) {
            # There is no count for this species/strain so just return as is
            return(x_i)
         } else {
            x_i_name = sub('^(.*)_\\[[0-9]+\\]$', '\\1', x_i)
            x_i_new = rep(x_i_name, as.integer(regex_count))
            x_i_new = paste(
               x_i_new,
               1:length(x_i_new),
               sep='_'
            )
            if (coll==TRUE) {
               return(paste(x_i_new, collapse=sep))
            } else {
               return(x_i_new)
            }

         }
      })
      if (coll==TRUE) {
         return(paste(x_split_expanded, collapse=sep))
      } else {
         return(unname(x_split_expanded))
      }
   }

   vector_split = function (x, sep=',') {
      strsplit(x, ',')[[1]]
   }

   species_relabel = function (x, sep=',') {
      print_strains(vector_split(species_expand(x, sep=sep)))
   }


   # Extract sequences
   qseq_1.dt = na.omit(unique(osu_seq_1.dt[, .(qseqid, sequence)])[order(qseqid)])
   qseq_2.dt = na.omit(unique(osu_seq_2.dt[, .(qseqid, sequence)])[order(qseqid)])

   pool_blast_result = pool_blast(qseq_1.dt, qseq_2.dt, temp_dir=temp_dir,
                                  verbose=verbose)
   blast4.dt = pool_blast_result[[1]]
   shorter_set = pool_blast_result[[2]]
   rm(pool_blast_result)

   # New code
   names(blast4.dt)[1:2] = c('qseqid_2', 'qseqid_1')
   names(blast4.dt)[3:4] = c('q2len', 'q1len')
   names(blast4.dt)[6:9] = c('q2start', 'q2end', 'q1start', 'q1end')

   osu_seq_map_1.dt = unique(osu_seq_1.dt[
      , .(osu_id_1=osu_id, qseqid_1=qseqid, copy_number_1=copy_number,
          sequence_1=sequence)])
   osu_seq_map_1.dt[, osu_id_1_seqs_n := lu(qseqid_1), by=osu_id_1]

   osu_seq_map_2.dt = unique(osu_seq_2.dt[
      , .(osu_id_2=osu_id, qseqid_2=qseqid, copy_number_2=copy_number,
          sequence_2=sequence)])
   osu_seq_map_2.dt[, osu_id_2_seqs_n := lu(qseqid_2), by=osu_id_2]

   # in blast4.dt, qseqid is qseqid from Set _2 and sseqid is qseqid from Set _1
   osu_seq_maps_12.dt = merge(
      blast4.dt, osu_seq_map_1.dt[osu_id_1 < osu_offset], by='qseqid_1'
   )
   osu_seq_maps_12.dt = merge(
      osu_seq_maps_12.dt, osu_seq_map_2.dt[osu_id_2 < osu_offset], by='qseqid_2'
   )
   osu_seq_maps_12.dt[nchar(sequence_1) >= nchar(sequence_2), longer_seq := 1]
   osu_seq_maps_12.dt[nchar(sequence_1) < nchar(sequence_2), longer_seq := 2]

   # osu_seq_maps_12.dt = na.omit(osu_seq_maps_12.dt)

   osu_seq_maps_12_2.dt = osu_seq_maps_12.dt[
      , if (lu(qseqid_2) == osu_id_2_seqs_n[1] &
            lu(qseqid_1) == osu_id_1_seqs_n[1]) .SD, by=.(osu_id_1, osu_id_2)]

   # Generate BLAST DB from the first set, FASTA File from the second set then run
   # BLAST.


   vcat('* preparing blast tables..')
   osu_seq_map_1.dt = unique(osu_seq_1.dt[
      , .(osu_id, qseqid, copy_number, sequence)])

   osu_seq_map_2.dt = unique(osu_seq_2.dt[
      , .(osu_id, qseqid, copy_number, sequence)])

   # in blast4.dt, qseqid is qseqid from Set _2 and sseqid is qseqid from Set _1


   # Merge Set 1 with BLAST results, find shared OSUs
   osu_seq_map_blast_1.dt = merge(
      osu_seq_map_1.dt,
      blast4.dt[, .(sseqid, qseqid2=qseqid, qseqid_new, qlen, slen,
                    qstart, qend, sstart, send, length)],
      # blast4.dt,
      by.x='qseqid', by.y='sseqid', all.x=T, allow.cartesian=T
   )
   # osu_seq_map_blast_1.dt[, sequence := substring(sequence, sstart[1], send[1]),
   #                        by=qseqid]
   osu_seq_map_blast_1.dt = osu_seq_map_blast_1.dt[
      , .SD[qlen==max(qlen)], by=.(qseqid, osu_id)]
   osu_seq_map_blast_1.dt[
      , sequence := substring(sequence, max(sstart), min(send)),
                          by=qseqid]
   ## osu_seq_map_blast_1.dt = osu_seq_map_blast_1.dt[, if(all(!is.na(sequence))) .SD, by=osu_id]
   osu_seq_map_blast_1a.dt = copy(na.omit(osu_seq_map_blast_1.dt))

   vcat('.')
   osu_seq_map_blast_2.dt = merge(
      osu_seq_map_2.dt,
      blast4.dt[, .(qseqid, qseqid1=sseqid, qseqid_new, slen, qlen,
                    qstart, qend, sstart, send, length)],
      by.x='qseqid', by.y='qseqid', all.x=T, allow.cartesian=T
   )
   osu_seq_map_blast_2.dt = osu_seq_map_blast_2.dt[
      , .SD[slen==max(slen)], by=.(qseqid, osu_id)]

   osu_seq_map_blast_2.dt[
      , sequence := substring(sequence,  max(qstart), min(qend)),
      by=qseqid]

   osu_seq_map_blast_2a.dt = copy(na.omit(osu_seq_map_blast_2.dt))
   vcat('OK.\n')

   vcat('* collapsing OSUs in set _1 ')
   seqs_with_diff_qseqids = osu_seq_map_blast_1a.dt[
      , if (lu(qseqid) > 1) .SD, by=sequence][, unique(sequence)]

   #-------------------- COLLAPSE SEQ FOR SET _1 -=============================
   seqs_with_diff_qseqids = osu_seq_map_blast_2a.dt[
      , if (lu(qseqid) > 1) .SD, by=sequence][, unique(sequence)]


   # osu_seq_map_blast_2.dt[, osu_n_qseqids := lu(qseqid), by=osu_id]
   # osu_seq_map_blast_2.dt[, osu_n_qseqid_news := lu(qseqid_new), by=osu_id]
   # osu_seq_map_blast_2.dt[, copy_number := sum(copy_number), by=.(osu_id, qseqid_new)]
   # # osu_seq_map_blast_2.dt[osu_n_qseqid_news==1, osu_id , by=sequence]
   #
   # qseqid_new_collapse = osu_seq_map_blast_2.dt[
   #    osu_n_qseqid_news==1 & osu_n_qseqids != osu_n_qseqid_news,
   #    unique(qseqid_new)
   #    ]
   #
   # for (qn in qseqid_new_collapse) {
   #    qn_collapse_osus = osu_seq_map_blast_2.dt[
   #       osu_n_qseqid_news==1 & qseqid_new==qn,
   #       sort(unique(osu_id))
   #       ]
   #    qn_collapse_osus_best = qn_collapse_osus[qn_collapse_osus < osu_offset]
   #    if (length(qn_collapse_osus_best) == 0) {
   #       qn_collapse_osus_best = qn_collapse_osus
   #    }
   #    qn_new_osu = min(qn_collapse_osus_best)
   #    # Update the sequence map
   #    osu_seq_map_blast_2.dt[
   #       osu_n_qseqid_news==1 & qseqid_new==qn,
   #       osu_id := qn_new_osu
   #       ]
   #    # Update the sequence descriptions for that OSU
   #    qn_new_pctsim = osu_seq_2.dt[osu_id %in% qn_collapse_osus_best,
   #                                 max(pctsim, na.rm=T)]
   #    old_coll_species = osu_seq_2.dt[osu_id %in% qn_collapse_osus_best,
   #                                    paste(species, collapse=',')]
   #    old_species = strsplit(old_coll_species, ',')[[1]]
   #    qn_new_species = print_strains(old_species)
   #    osu_seq_2.dt[
   #       osu_id==qn_new_osu,
   #       c('species', 'pctsim') := list(qn_new_species, qn_new_pctsim)]
   # }

   for (seq in seqs_with_diff_qseqids) {
      # vcat('.')
      min_id = osu_seq_map_blast_1a.dt[sequence==seq, which(osu_id==min(osu_id))[1]]
      # min_osuid = osu_seq_map_blast_1a.dt[sequence==seq, min(osu_id)]
      all_osuids = osu_seq_map_blast_1a.dt[sequence==seq, osu_id]
      min_cp = osu_seq_map_blast_1a.dt[sequence==seq, .SD[min_id, copy_number]]
      # min_cp = max(
      #    osu_seq_map_blast_1a.dt[sequence==seq, .SD[min_id, copy_number]],
      #    osu_seq_map_blast_2a.dt[
      #       sequence==seq &
      #          qseqid==osu_seq_map_blast_1a.dt[sequence==seq, .SD[min_id, qseqid2]],
      #       copy_number
      #       ]
      # )

      # min_pctsim = osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, pctsim]]
      # min_species = osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, species]]
      min_osuid = min(all_osuids)
      osu_seq_map_blast_1.dt[
         sequence==seq, c('osu_id', 'copy_number', 'qseqid_new') := list(
            .SD[min_id, osu_id],
            .SD[min_id, copy_number],
            # min_cp,
            .SD[min_id, qseqid_new]
         )
         ]
      # Update abundance table by pooling OSU counts
      # If we already have single OSU with multiple sequences (spectra)
      # then we collapse to the first sequence
      # Update sequence reference table
      # osu_seq_1.dt[
      #    osu_id %in% all_osuids,
      #    c('osu_id', 'copy_number') := list(min_osuid, min_cp)
      # ]
   }
   osu_seq_map_blast_1.dt = unique(osu_seq_map_blast_1.dt)
   vcat(' OK.\n')

   osu_spectra_1.dt = osu_seq_map_blast_1.dt[
      , .(
         spectrum = paste(
            sort(unique(
               paste(sequence, copy_number, sep=':')
            )),
            collapse=','
          ),
          len_1=slen[1]),
      by=osu_id
   ]

   #------------------ Collapse sequences for Set _2 ---------------------------
   vcat('* collapsing OSUs in set _2 ')

   seqs_with_diff_qseqids = osu_seq_map_blast_2a.dt[
      , if (lu(qseqid) > 1) .SD, by=sequence][, unique(sequence)]

#
#    osu_seq_map_blast_2.dt[, osu_n_qseqids := lu(qseqid), by=osu_id]
#    osu_seq_map_blast_2.dt[, osu_n_qseqid_news := lu(qseqid_new), by=osu_id]
#    osu_seq_map_blast_2.dt[, copy_number := sum(copy_number), by=.(osu_id, qseqid_new)]
#    # osu_seq_map_blast_2.dt[osu_n_qseqid_news==1, osu_id , by=sequence]
#
#    qseqid_new_collapse = osu_seq_map_blast_2.dt[
#       osu_n_qseqid_news==1 & osu_n_qseqids != osu_n_qseqid_news,
#       unique(qseqid_new)
#    ]
#
#    for (qn in qseqid_new_collapse) {
#       qn_collapse_osus = osu_seq_map_blast_2.dt[
#          osu_n_qseqid_news==1 & qseqid_new==qn,
#          sort(unique(osu_id))
#       ]
#       qn_collapse_osus_best = qn_collapse_osus[qn_collapse_osus < osu_offset]
#       if (length(qn_collapse_osus_best) == 0) {
#          qn_collapse_osus_best = qn_collapse_osus
#       }
#       qn_new_osu = min(qn_collapse_osus_best)
#       # Update the sequence map
#       osu_seq_map_blast_2.dt[
#          osu_n_qseqid_news==1 & qseqid_new==qn,
#          osu_id := qn_new_osu
#       ]
#       # Update the sequence descriptions for that OSU
#       qn_new_pctsim = osu_seq_2.dt[osu_id %in% qn_collapse_osus_best,
#                                    max(pctsim, na.rm=T)]
#       old_coll_species = osu_seq_2.dt[osu_id %in% qn_collapse_osus_best,
#                                       paste(species, collapse=',')]
#       old_species = strsplit(old_coll_species, ',')[[1]]
#       qn_new_species = print_strains(old_species)
#       osu_seq_2.dt[
#          osu_id==qn_new_osu,
#          c('species', 'pctsim') := list(qn_new_species, qn_new_pctsim)]
#    }

   # osu_seq_map_blast_2.dt[osu_n_qseqid_news]

      for (seq in seqs_with_diff_qseqids) {
         # vcat('.', which(seqs_with_diff_qseqids == seq), '.')
         min_id = osu_seq_map_blast_2.dt[sequence==seq, which(osu_id==min(osu_id))[1]]
         all_osuids = osu_seq_map_blast_2.dt[sequence==seq, unique(osu_id)]
         min_cp = osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, copy_number]]
         # Always choose the larger copy number between the two matched sequences
         # min_cp = max(
         #    osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, copy_number]],
         #    osu_seq_map_blast_1.dt[
         #       sequence==seq &
         #       qseqid==osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, qseqid1]],
         #       copy_number
         #    ]
         # )
         # min_pctsim = osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, pctsim]]
         # min_species = osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, species]]
         min_osuid = min(all_osuids)
         osu_seq_map_blast_2.dt[
            sequence==seq, c('osu_id', 'copy_number', 'qseqid_new') := list(
               .SD[min_id, osu_id],
               .SD[min_id, copy_number],
               # min_cp,
               .SD[min_id, qseqid_new]
            )
         ]
         # Update abundance table by pooling OSU counts
         # If we already have single OSU with multiple sequences (spectra)
         # then we collapse to the first sequence
         # osu_ab_2.dt[
         #    osu_id %in% all_osuids,
         #    c('osu_id', 'osu_count', 'pctsim', 'species') := list(
         #       min_osuid, sum(.SD[, osu_count[1], by=osu_id]),
         #       .SD[osu_id==min_osuid, pctsim],
         #       .SD[osu_id==min_osuid, species]), by=.(sample_id, dataset)]
         # Update sequence reference table
         # osu_seq_2.dt[
         #    osu_id %in% all_osuids,
         #    c('osu_id', 'copy_number') := list(min_osuid, min_cp)
         # ]
      }
   osu_seq_map_blast_2.dt = unique(osu_seq_map_blast_2.dt)
   vcat(' OK.\n')

   ## osu_seq_map_blast_2.dt = osu_seq_map_blast_2.dt[, if(all(!is.na(sequence))) .SD, by=osu_id]
   osu_spectra_2.dt = osu_seq_map_blast_2.dt[
      , .(
         spectrum = paste(
            sort(unique(
               paste(sequence, copy_number, sep=':')
            )),
            collapse=','
         ),
         len_2=qlen[1]
      ),by=osu_id
   ]

   osu_matching_spectra.dt = merge(
      osu_spectra_1.dt[, .(osu_id_1=osu_id, spectrum, len_1)],
      osu_spectra_2.dt[, .(osu_id_2=osu_id, spectrum, len_2)],
      by='spectrum'
   )[order(osu_id_1)]

   # Renumerate matching OSUs
   # First exact matches
   if (shorter_set == 1) {
      # osu_matching_spectra.dt[osu_id_1 < osu_offset,
      #                         osu_id_new := .GRP, by=spectrum]
      # osu_matching_spectra.dt[osu_id_1 >= osu_offset,
      #                         osu_id_new := .GRP + as.integer(osu_offset),
      #                         by=spectrum]
      osu_matching_spectra.dt[
         , osu_id_new := ifelse(
            osu_id_1 < osu_offset,
            .GRP,
            .GRP + as.integer(osu_offset)
         )
         , by=spectrum]
   } else {  # shorter_set == 2
      # osu_matching_spectra.dt[osu_id_2 < osu_offset,
      #                         osu_id_new := .GRP, by=spectrum]
      # osu_matching_spectra.dt[osu_id_2 >= osu_offset,
      #                         osu_id_new := .GRP + as.integer(osu_offset),
      #                         by=spectrum]
      osu_matching_spectra.dt[
         , osu_id_new := ifelse(
            osu_id_2 < osu_offset,
            .GRP,
            .GRP + as.integer(osu_offset)
         )
         , by=spectrum]
   }

   # After sequence trimming sometimes < 100% match can collapse to 100% match
   # so the mapping will go to the higher pctsim OSU
   osu_matching_spectra.dt[, osu_id_new := min(osu_id_new), by=osu_id_1]
   osu_matching_spectra.dt[, osu_id_new := min(osu_id_new), by=osu_id_2]

   # Renumerate OSUs from _1 that are not matching
   vcat('* renumerating OSUs _1 ...')
   osu_renum_1.dt = unique(osu_seq_1.dt[, .(osu_id_1 = osu_id)])[
      !(osu_id_1 %in% osu_matching_spectra.dt[, osu_id_1])]
   osu_100_last = osu_matching_spectra.dt[osu_id_new < as.integer(osu_offset), max(osu_id_new)]
   osu_non100_last = osu_matching_spectra.dt[osu_id_new >= as.integer(osu_offset), max(osu_id_new)]
   osu_renum_1.dt[osu_id_1 < osu_offset,
                  osu_id_new := .GRP + osu_100_last, by=osu_id_1]
   osu_renum_1.dt[osu_id_1 >= osu_offset,
                  osu_id_new := .GRP + osu_non100_last, by=osu_id_1]
   osu_renum_all_1.dt = data.table::rbindlist(list(
      osu_renum_1.dt, unique(osu_matching_spectra.dt[, .(osu_id_1, osu_id_new)])
   ))[order(osu_id_new)]
   vcat(' OK.\n')
   # Renumerate OSUs from _2 that are not matching
   vcat('* renumerating OSUs _2 ...')
   osu_renum_2.dt = unique(osu_seq_2.dt[, .(osu_id_2 = osu_id)])[
      !(osu_id_2 %in% osu_matching_spectra.dt[, osu_id_2])]
   osu_100_last = osu_renum_1.dt[osu_id_new < osu_offset, max(osu_id_new)]
   osu_non100_last = osu_renum_1.dt[osu_id_new >= osu_offset, max(osu_id_new)]
   osu_renum_2.dt[osu_id_2 < as.integer(osu_offset),
                  osu_id_new := .GRP + osu_100_last, by=osu_id_2]
   osu_renum_2.dt[osu_id_2 >= as.integer(osu_offset),
                  osu_id_new := .GRP + osu_non100_last, by=osu_id_2]
   osu_renum_all_2.dt = data.table::rbindlist(list(
      osu_renum_2.dt, unique(osu_matching_spectra.dt[, .(osu_id_2, osu_id_new)])
   ))[order(osu_id_new)]
   vcat(' OK.\n')


   # Add new OSU IDs to abundance tables
   vcat('* adding new OSUs .')
   osu_ab_1_new.dt = merge(osu_ab_1.dt, osu_renum_all_1.dt, by.x='osu_id',
                           by.y='osu_id_1', all.x=T, all.y=T)

   osu_seq_1_new.dt = merge(osu_seq_1.dt, osu_renum_all_1.dt, by.x='osu_id',
                           by.y='osu_id_1', all.x=T, all.y=T)

   # Update non-matching OSUs with trimmed sequences, e.g.
   # we have a matching spectra 5:s1, 1:s2, 1:s3
   # but _2 contains 1:s1 spectra that is not matched in _1
   osu_seq_1_new.dt = merge(
      osu_seq_1_new.dt,
      unique(osu_seq_map_blast_1.dt[, .(qseqid, seq_new=sequence)]),
      by='qseqid', all.x=T
   )
   osu_seq_1_new.dt[!is.na(seq_new) & sequence != seq_new,
                    sequence := seq_new]
   osu_seq_1_new.dt[, seq_new := NULL]
   setcolorder(osu_seq_1_new.dt, c('osu_id', 'species', 'qseqid', 'copy_number',
                                   'pctsim', 'sequence', 'osu_id_new'))


   # If _2 set is shorter then update the sequences here in _1 with the shorter
   # ones from _2, for matching OSUs.

   if (osu_ab_1_new.dt[is.na(osu_id), .N] > 0 |
       osu_ab_1_new.dt[is.na(osu_id_new), .N] > 0) {
      stop('Error: Missing or extra OSU IDs in _1 after renumeration.')
   }
   osu_ab_1_new.dt[, osu_id := osu_id_new]
   osu_ab_1_new.dt[, osu_id_new := NULL]
   osu_ab_1_new.dt[, osu_count_new := sum(osu_count), by=.(osu_id, sample_id)]
   # Pool osu_counts for OSUs that have been collapsed because they exactly
   # matched shorter sequence from the other dataset.
   osu_ab_1_new.dt[, osu_count := osu_count_new]
   osu_ab_1_new.dt[, c('osu_count_new', 'pctsim', 'species') := NULL]
   osu_ab_1_new.dt = unique(osu_ab_1_new.dt)

   vcat('.')

   # DO THE SAME WITH OSU_SEQs
   osu_ab_2_new.dt = merge(osu_ab_2.dt, osu_renum_all_2.dt, by.x='osu_id',
                           by.y='osu_id_2', all.x=T, all.y=T)
   osu_seq_2_new.dt = merge(osu_seq_2.dt, osu_renum_all_2.dt, by.x='osu_id',
                            by.y='osu_id_2', all.x=T, all.y=T)

   # Update non-matching OSUs with trimmed sequences, e.g.
   # we have a matching spectra 5:s1, 1:s2, 1:s3
   # but _2 contains 1:s1 spectra that is not matched in _1
   osu_seq_2_new.dt = merge(
      osu_seq_2_new.dt,
      unique(osu_seq_map_blast_2.dt[, .(qseqid, seq_new=sequence)]),
      by='qseqid', all.x=T
   )
   osu_seq_2_new.dt[!is.na(seq_new) & sequence != seq_new,
                    sequence := seq_new]
   osu_seq_2_new.dt[, seq_new := NULL]
   setcolorder(osu_seq_2_new.dt, c('osu_id', 'species', 'qseqid', 'copy_number',
                                   'pctsim', 'sequence', 'osu_id_new'))

   if (osu_ab_2_new.dt[is.na(osu_id), .N] > 0 |
       osu_ab_2_new.dt[is.na(osu_id_new), .N] > 0) {
      stop('Error: Missing or extra OSU IDs in _2 after renumeration.')
   }
   osu_ab_2_new.dt[, osu_id := osu_id_new]
   osu_ab_2_new.dt[, osu_id_new := NULL]
   osu_ab_2_new.dt[, osu_count_new := sum(osu_count), by=.(osu_id, sample_id)]
   # Pool osu_counts for OSUs that have been collapsed because they exactly
   # matched shorter sequence from the other dataset.
   osu_ab_2_new.dt[, osu_count := osu_count_new]
   osu_ab_2_new.dt[, c('osu_count_new', 'pctsim', 'species') := NULL]
   osu_ab_2_new.dt = unique(osu_ab_2_new.dt)

   vcat('.')

   # Now sort out OSU_SEQ tables
   # From the matching table
   # If _1 is longer than _2 then add shorter seqs of _2 to _1

   osu_seq_1_update.dt = merge(
      osu_matching_spectra.dt[len_1 >= len_2],
      osu_seq_2_new.dt,
      by.x=c('osu_id_2', 'osu_id_new'),
      by.y=c('osu_id', 'osu_id_new')
   )

   osu_seq_1_update.dt = osu_seq_1_update.dt[
      , .(osu_id=osu_id_1, species, qseqid, copy_number, pctsim, sequence,
          osu_id_new)]

   osu_seq_1_new_all.dt = data.table::rbindlist(list(
      osu_seq_1_new.dt[!(osu_id %in% osu_seq_1_update.dt[, unique(osu_id)])],
      osu_seq_1_update.dt
   ))

   vcat('.')


   # 1 is shorter, so update 2 with seqs from 1
   osu_seq_2_update.dt = merge(
      osu_matching_spectra.dt[len_1 < len_2],
      osu_seq_1_new.dt,
      by.x=c('osu_id_1', 'osu_id_new'),
      by.y=c('osu_id', 'osu_id_new')
   )

   osu_seq_2_update.dt = osu_seq_2_update.dt[
      , .(osu_id=osu_id_2, species, qseqid, copy_number, pctsim, sequence,
          osu_id_new)]

   osu_seq_2_new_all.dt = data.table::rbindlist(list(
      osu_seq_2_new.dt[!(osu_id %in% osu_seq_2_update.dt[, unique(osu_id)])],
      osu_seq_2_update.dt
   ))
   vcat('.')
   osu_seq_1_new_all_f.dt = copy(osu_seq_1_new_all.dt)
   osu_seq_1_new_all_f.dt[, osu_id := osu_id_new]
   osu_seq_1_new_all_f.dt[, osu_id_new := NULL]
   osu_seq_2_new_all_f.dt = copy(osu_seq_2_new_all.dt)
   osu_seq_2_new_all_f.dt[, osu_id := osu_id_new]
   osu_seq_2_new_all_f.dt[, osu_id_new := NULL]

   #------------------ Finalize OSU_SEQ TABLE ----------------------------------
   osu_seq_new.dt = unique(data.table::rbindlist(list(
      osu_seq_1_new_all_f.dt,
      osu_seq_2_new_all_f.dt
   )))
   # Finally, renumerate all qseqid by sequence
   osu_seq_new.dt[, qseqid := NULL]
   osu_seq_new.dt[, qseqid := .GRP, by=sequence]

   osu_seq_1_ref.dt = unique(osu_seq_1_new_all_f.dt[, .(osu_id, pctsim, species)])
   osu_seq_1_ref.dt = osu_seq_1_ref.dt[
      , .(species = print_strains(species),
          pctsim = unique(pctsim[!is.na(pctsim)])[1]), by=osu_id]
   osu_seq_2_ref.dt = unique(osu_seq_2_new_all_f.dt[, .(osu_id, pctsim, species)])
   osu_seq_2_ref.dt = osu_seq_2_ref.dt[
      , .(species = print_strains(species),
          pctsim = unique(pctsim[!is.na(pctsim)])[1]), by=osu_id]
   # Add new OSU annotations to OSU AB table
   osu_ab_1_new.dt = merge(
      osu_ab_1_new.dt,
      osu_seq_1_ref.dt,   # <-------------- change this
      by='osu_id' #, allow.cartesian=T
   )
   setcolorder(osu_ab_1_new.dt, osu_ab_1_colorder)
   osu_ab_2_new.dt = merge(
      osu_ab_2_new.dt,
      osu_seq_2_ref.dt,
      #unique(osu_seq_2_new_all_f.dt[, .(osu_id, pctsim, species)]),
      by='osu_id'# , allow.cartesian=T
   )
   setcolorder(osu_ab_2_new.dt, osu_ab_1_colorder)

   vcat(' OK.\n')

   vcat('* merging output...')
   osu_ab_new.dt = data.table::rbindlist(list(osu_ab_1_new.dt, osu_ab_2_new.dt))
   # DO THE SAME WITH OSU_SEQ TABLES
   vcat(' OK.\n')

   vcat('* finalizing output...')

   osu_seq_new.dt[, nqs := lu(qseqid), by=osu_id]
   seqs_collapse_final = osu_seq_new.dt[
      nqs==1, if (lu(osu_id) > 1) 1, by=sequence][, unique(sequence)]
   for (s in seqs_collapse_final) {
      vcat('\n process seq ', which(seqs_collapse_final==s), '...')
      osu_seq_new_sub.dt = osu_seq_new.dt[sequence==s & nqs==1]
      s_osuids = osu_seq_new_sub.dt[, sort(unique(osu_id))]
      s_final_pctsim = osu_seq_new_sub.dt[, max(pctsim, na.rm=T)]
      s_final_osuid = osu_seq_new_sub.dt[pctsim==s_final_pctsim,
         min(osu_id)
      ]
      s_final_species = osu_seq_new_sub.dt[pctsim==s_final_pctsim,
         species_relabel(species)
      ]
      s_final_cp = osu_seq_new_sub.dt[
         pctsim==s_final_pctsim & osu_id==s_final_osuid,
         max(copy_number)
      ]
      osu_seq_new.dt[
         osu_id %in% s_osuids & nqs==1 & sequence==s,
         c('osu_id', 'species', 'pctsim', 'copy_number') := list(
            s_final_osuid, s_final_species, s_final_pctsim, s_final_cp
         )
      ]
      osu_ab_new.dt[
         osu_id %in% s_osuids,
         c('osu_id', 'osu_count', 'pctsim', 'species') := list(
            s_final_osuid, sum(osu_count), s_final_pctsim, s_final_species
         ),
         by=.(dataset, sample_id)
      ]
   }
   osu_seq_new.dt = unique(osu_seq_new.dt)
   osu_seq_new.dt[, nqs := NULL]
   osu_ab_new.dt = unique(osu_ab_new.dt)

   vcat(' OK.\n')

   return(list(osu_ab_new.dt, osu_seq_new.dt))
}
