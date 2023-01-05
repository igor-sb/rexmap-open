#' Merge paired-end reads
#'
#' Merges reads using overlap version of the Needleman-Wunsch alignment
#' algorithm and calculate merged quality scores based on the posterior
#' probabilities (cite: Edgar, Fyelberg)
#'
#' @param fq_fwd A character vector for input forward FASTQ files (can be
#'   gzipped).
#' @param fq_rev A character vector for input reverse FASTQ files (can be
#'   gzipped).
#' @param fq_mer A character vector for output merged FASTQ files.
#' @param min_sim Minimum similarity to accept an alignment. A floating point
#'   number between 0 and 1. Default: 0.75, corresponding to minimum 75\%
#'   similarity in the aligned overlap region.
#' @param min_aln_len Minimum alignment length. Ignore alignments below this
#'   alignment length.
#' @param match Score for character match in the alignment (default: 5.
#' @param mismatch Score for character mismatch in the alignment (default: -5)
#' @param gap_p Score for gap penalty, either gap opening on gap extension
#'   (default: -7).
#' @param rc_reverse TRUE/FALSE Reverse complement reverse reads? (default:
#'   TRUE)
#' @param ncpu Number of CPU threads to use for multithreading.
#' @param verbose TRUE/FALSE Display of status messages.
#' @param timing TRUE/FALSE Time merging.
#'
#' @export
merge_pairs <- Vectorize(function(
    # fq_fwd,
    # fq_rev,
    # fq_mer,
    # min_sim = 0.75,
    # min_aln_len = 50,
    forward_file,
    reverse_file,
    merged_file,
    transform_reverse_complement = list("forward" = FALSE, "reverse" = TRUE),
    min_acceptance_criteria = list(similarity = 0.75, overlap = 50),
    alignment_scores = list(match = 7L, mismatch = -5L, gap_penalty = -7L),
    number_threads = rexmap_option("ncpu"),
    logger = rexmap_option("logger")
    # rc_reverse = TRUE,
    # ncpu = rexmap_option("ncpu"),
    # verbose = FALSE,
    # timing = FALSE
) {
  check_merge_pairs_args(forward_file, reverse_file, merged_file, n_threads)
  
  log4r::info(
    logger,
    "Loading FASTQs: ", basename(forward_file), ", ", basename(reverse_file)
  )
  
  withr::with_connection(
    list("forward_stream" = ShortRead::FastqStreamer(forward_file),
         "reverse_stream" = ShortRead::FastqStreamer(reverse_file)),
    {
      repeat {
        seqs_quals_ids <- get_sequences_qualities_ids(
          list("forward" = forward_stream, "reverse" = reverse_stream),
          reverse_complement
        )
        if (length(seqs_quals_ids$ids) == 0) break

        merged_sequences <- seqs_quals_ids |>
          merge_sequences(alignment_scores, number_threads) |>
          filter_merged_sequences(min_acceptance_criteria)
        
        write_chunks_to_merged_file(merged_file, merged_sequences)
        calculate_merge_stats(merged_sequences)
      }
    }
  )
  merged_sread <- ShortRead::ShortReadQ(
    sread = Biostrings::DNAStringSet(final_seqs),
    quality = Biostrings::BStringSet(final_qual),
    id = Biostrings::BStringSet(final_ids)
  )

  # Generate statistics for filtered-out reads
  stats <- list(
    "total" = length(merged_aln_filter),
    "low_pct_sim" = length(merged_aln_filter) - sum(as.logical(merged_list[3, ])),
    "low_aln_len" = length(merged_aln_filter) - sum(as.logical(merged_list[4, ]))
  )

  # Free memory and close file connections
  rm(merged_list)

  # If file exists, delete it, then write new one
  if (file.exists(fq_mer)) file_remove_result <- file.remove(fq_mer)

  ShortRead::writeFastq(merged_sread, fq_mer, compress = F)
  return(stats)
}, vectorize.args = c("forward_file", "reverse_file", "merged_file"))


filter_merged_sequences <- function(merged_sequences, min_acceptance_criteria) {
  is_merging_accepted <- lapply(
    merged_chunks,
    function(chunk) {
      chunk$similarity >= min_acceptance_criteria$similarity &&
        chunk$overlap >= min_acceptance_criteria$overlap &&
        !is.null(chunk$sequence) &&
        !is.null(chunk$quality)
    }
  )
  merged_sequences[is_merging_accepted]
}

merge_sequences <- function(seqs_quals_ids, alignment_scores, number_threads) {
  qual_files <- get_merged_qualities_files()
  parallel::mcmapply(
    align_seqs_and_quals,
    seqs_quals_ids$sequences,
    seqs_quals_ids$qualities,
    alignment_scores,
    posterior_match_file = qual_files$matches,
    posterior_mismatch_file = qual_files$mismatches,
    mc.cores = number_threads
  )
}

get_merged_qualities_files <- function() {
  match_file <- system.file(
    "merge_tables",
    "himap_mergepairs_match_qs.txt",
    package = pname_l
  )
  mismatch_file <- system.file(
    "merge_tables",
    "himap_mergepairs_mismatch_qs.txt",
    package = pname_l
  )
  c("matches" = match_file, "mismatches" = mismatch_file)
}

check_merge_pairs_args <- function(
    forward_file,
    reverse_file,
    merged_file,
    number_threads
) {
  stopifnot(
    "forward_file does not exist" = file.exists(forward_file),
    "reverse_file does not exist" = file.exists(reverse_file),
    "merged_file path does not exist" = dir.exists(dirname(merged_file)),
    "invalid number_threads" = number_threads >= 1
  )
}

# Reads chunks from forward and reverse streams in parallel
get_fastq_entries_from_streams <- function(
    stream_list,
    transform_revcomp_list
) {
  mapply(
    function(stream, transform_reverse_complement) {
      sequence_data <- ShortRead::yield(stream)
      if (transform_reverse_complement) {
        sequence_data <- ShortRead::reverseComplement(sequence_data)
      }
      sequence_data
    },
    stream_list,
    transform_revcomp_list,
    SIMPLIFY = FALSE
  )
}

get_sequences_qualities_ids <- function(stream_list, transform_revcomp_list) {
  paired_fastq_entries <- get_fastq_entries_from_streams(
    stream_list, transform_revcomp_list
  )
  if (length(paired_fastq_entries$forward) == 0 ||
      length(paired_fastq_entries$reverse) == 0) {
    return(list(sequences = NULL, qualities = NULL, ids = NULL))
  }
  lapply(
    list("sequences", "qualities", "ids"),
    function(type) get_data_from_fastq_entries(paired_fastq_entries, type)
  )
}

get_data_from_fastq_entries <- function(paired_fastq_entries, data_type) {
  stopifnot(
    "Invalid data type." = data_type %in% c("sequences", "qualities", "ids")
  )
  function_name <- paste0("get_", data_type, "_from_chunk")
  lapply(paired_fastq_entries, get(function_name)) |> data.table::transpose()
}

get_data_from_chunks <- function(chunk_list) {
  lapply(
    function(chunk) {
      list(
        ids = get_ids_from_chunk(chunk),
        sequences = get_sequences_from_chunk(chunk),
        qualities = get_qualities_from_chunk(chunk)
      )
    },
    chunk_list
  )
}

get_ids_from_sequence_data <- function(sequence_data) {
  raw_ids <- ShortRead::id(sequence_data) |> as.character()
  gsub("^([^ ]+) .*", "\\1", raw_ids)
}

get_sequences_from_sequence_data <- function(sequence_data) {
  ShortRead::sread(sequence_data) |> as.character()
}

get_qualities_from_sequence_data <- function(sequence_data) {
  Biostrings::quality(sequence_data) |> Biostrings::quality() |> as.character()
}

#' Convert mergestats table to a normal data.table
#'
#' @export
mergeout_to_table <- function(mergestats) {
  mergestats.dt <- data.table::data.table(
    matrix(
      unlist(t(mergestats)),
      ncol = nrow(mergestats),
      dimnames = list(sampleids_from_filenames(colnames(mergestats), "_"), rownames(mergestats))
    ),
    keep.rownames = "sample_id"
  )
  return(mergestats.dt)
}
