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
    reverse_complement = list("forward" = FALSE, "reverse" = TRUE),
    min_acceptance_criteria = list(similarity = 0.75, overlap = 50),
    alignment_scores = list(match = 5L, mismatch = -5L, gap_penalty = -7L),
    logger = rexmap_option("logger"),
    n_threads = rexmap_option("ncpu")
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
        chunks <- get_chunks_from_streams(
          list("forward" = forward_stream, "reverse" = reverse_stream),
          reverse_complement
        )
        if (length(chunks$forward) == 0 || length(chunks$reverse) == 0) break
        
        merged_sequences <- chunks |>
          merge_sequences(alignment_scores, n_threads) |>
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

merge_chunks <- function(chunks, alignment_scores) {
  qual_files <- get_merged_qualities_files()
  parallel::mcmapply(
    C_mergepairs,
    chunk_data$forward$sequences,
    chunk_data$reverse$sequences,
    chunk_data$forward$qualities,
    chunk_data$reverse$qualities,
    match = alignment_scores$match,
    mismatch = alignment_scores$mismatch,
    gap_p = alignment_scores$gap_penalty,
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
get_chunks_from_streams <- function(stream_list, revcomp_transform_list) {
  mapply(
    function(stream, revcomp_transform) {
      chunk <- ShortRead::yield(stream)
      if (revcomp_transform) chunk <- ShortRead::reverseComplement(chunk)
    },
    stream_list,
    revcomp_transform_list,
    SIMPLIFY = FALSE
  )
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

get_ids_from_chunk <- function(chunk) {
  raw_ids <- ShortRead::id(chunk) |> as.character()
  gsub("^([^ ]+) .*", "\\1", raw_ids)
}

get_sequences_from_chunk <- function(chunk) {
  ShortRead::sread(chunk) |> as.character()
}

get_qualities_from_chunk <- function(chunk) {
  Biostrings::quality(chunk) |> Biostrings::quality() |> as.character()
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
