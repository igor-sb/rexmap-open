test_that("mergepairs c++: create scores", {
  expect_identical(
    create_scoring_matrix(5, -2),
    list(c(5L, -2L, -2L, -2L, -2L),
         c(-2L, 5L, -2L, -2L, -2L),
         c(-2L, -2L, 5L, -2L, -2L),
         c(-2L, -2L, -2L, 5L, -2L),
         c(-2L, -2L, -2L,-2L, -2L))
  )
})

test_that("mergepairs c++: to_int", {
  expect_identical(
    sapply(c("A", "C", "G", "T", "X"), to_int),
    c("A" = 0L, "C" = 1L, "G" = 2L, "T" = 3L, "X" = 4L)
  )
})

test_that("mergepairs c++: get indexes from score and path", {
  sequences <- c(forward = "ATGCGA", reverse = "CGATT")
  ncol <- nchar(sequences[1]) + 1L
  nrow <- nchar(sequences[2]) + 1L
  # C++ index 1, 2: 
  expect_identical(
    get_indexes(1, 2, ncol),
    c("diag" = 7, "left" = 14, "upper" = 8, "current" = 15)
  )
})

test_that("mergepairs c++: dp score = zero & path = left for first row", {
  sequences <- c(forward = "ATGCGA", reverse = "CGATT")
  ncol <- nchar(sequences[1]) + 1L
  nrow <- nchar(sequences[2]) + 1L
  score <- rep(NA_integer_, ncol * nrow)
  path <- rep(NA_integer_, ncol * nrow)
  score_and_path <- test_calc_score_path_first_row(score, path, ncol)
  
  expect_identical(
    dp_path_matrix(score_and_path$path, ncol)[1, ],
    rep("l", ncol)
  )
  expect_identical(
    dp_score_matrix(score_and_path$score, ncol)[1, ],
    rep(0L, ncol)
  )
})

test_that("mergepairs c++: dp score and path = up for first column", {
  sequences <- c(forward = "ATGCGA", reverse = "CGATT")
  scoring_matrix <- create_scoring_matrix(5L, -5L)
  gap_p <- -7L
  ncol <- nchar(sequences[1]) + 1L
  nrow <- nchar(sequences[2]) + 1L
  score <- rep(NA_integer_, ncol * nrow)
  path <- rep(NA_integer_, ncol * nrow)
  score_and_path <- test_calc_score_path_first_column(score, path, nrow, ncol,
                                                      gap_p)
  expect_identical(
    dp_path_matrix(score_and_path$path, ncol)[-1, 1],
    rep("u", nrow - 1)
  )
  expect_identical(
    dp_score_matrix(score_and_path$score, ncol)[-1, 1],
    seq.int(gap_p, gap_p * (nrow - 1L), by = gap_p)
  )
})

test_that("mergepairs c++: calc_score_path_other ", {
  sequences <- c(forward = "ATGCGAT", reverse = "CGGTTAC")
  scoring_matrix <- create_scoring_matrix(7L, -3L)
  gap_p <- -7L
  ncol <- nchar(sequences[1]) + 1L
  nrow <- nchar(sequences[2]) + 1L
  score <- rep(NA_integer_, ncol * nrow)
  path <- rep(NA_integer_, ncol * nrow)
  score[1:ncol] <- 0L
  path[1:ncol] <- chr_to_int("l")
  first_col <- seq.int(1L + ncol, 1L + ncol * (nrow - 1L), by = ncol)
  score[first_col] <- seq.int(gap_p, gap_p * (nrow - 1L), by = gap_p)
  path[first_col] <- chr_to_int("u")
  score_and_path <- test_calc_score_path_other(score, path, sequences,
                                               scoring_matrix, gap_p)
  
  true_scores <- matrix(
    c(
      0L,   0L,   0L,   0L,   0L,   0L,   0L,   0L,
      -7L,  -3L,  -3L,  -3L,   7L,   0L,  -3L,   0L,
      -14L, -10L,  -6L,   4L,   0L,  14L,   7L,   0L,
      -21L, -17L, -13L,   1L,   1L,   7L,  11L,   4L,
      -28L, -24L, -10L,  -6L,  -2L,   0L,   4L,  18L,
      -35L, -31L, -17L, -13L,  -9L,  -5L,  -3L,  18L,
      -42L, -28L, -24L, -20L, -16L, -12L,   2L,  18L,
      -49L, -35L, -31L, -27L, -13L, -19L,  -5L,  18L
    ),
    ncol = 8,
    byrow = TRUE
  ) |> t() |> as.vector()
  true_path <- matrix(
    c(
      "l", "l", "l", "l", "l", "l", "l", "l",
      "u", "d", "d", "d", "d", "l", "d", "u",
      "u", "u", "d", "d", "u", "d", "l", "u",
      "u", "u", "u", "d", "d", "u", "d", "l",
      "u", "u", "d", "u", "d", "u", "u", "d",
      "u", "u", "u", "u", "u", "d", "u", "u",
      "u", "d", "u", "u", "u", "u", "d", "u",
      "u", "u", "u", "u", "d", "u", "u", "u"
    ),
    ncol = 8,
    byrow = TRUE
  ) |> t() |> as.vector() |> chr_to_int()
  expect_identical(score_and_path$score, true_scores)
  expect_identical(score_and_path$path, true_path)
})

test_that("mergepairs c++: find_best_scoring_path w/ mismatch no indel", {
  # ATGCGAT
  #    CGGTTAC
  sequences <- c(forward = "ATGCGAT", reverse = "CGGTTAC")
  scoring_matrix <- create_scoring_matrix(7L, -3L)
  gap_p <- -7L
  score_and_path <- test_find_best_scoring_path(
    sequences, scoring_matrix, gap_p
  )
  
  true_scores <- matrix(
    c(
        0L,   0L,   0L,   0L,   0L,   0L,   0L,   0L,
       -7L,  -3L,  -3L,  -3L,   7L,   0L,  -3L,   0L,
      -14L, -10L,  -6L,   4L,   0L,  14L,   7L,   0L,
      -21L, -17L, -13L,   1L,   1L,   7L,  11L,   4L,
      -28L, -24L, -10L,  -6L,  -2L,   0L,   4L,  18L,
      -35L, -31L, -17L, -13L,  -9L,  -5L,  -3L,  18L,
      -42L, -28L, -24L, -20L, -16L, -12L,   2L,  18L,
      -49L, -35L, -31L, -27L, -13L, -19L,  -5L,  18L
    ),
    ncol = 8,
    byrow = TRUE
  ) |> t() |> as.vector()
  true_path <- matrix(
    c(
      "l", "l", "l", "l", "l", "l", "l", "l",
      "u", "d", "d", "d", "d", "l", "d", "u",
      "u", "u", "d", "d", "u", "d", "l", "u",
      "u", "u", "u", "d", "d", "u", "d", "l",
      "u", "u", "d", "u", "d", "u", "u", "d",
      "u", "u", "u", "u", "u", "d", "u", "u",
      "u", "d", "u", "u", "u", "u", "d", "u",
      "u", "u", "u", "u", "d", "u", "u", "u"
    ),
    ncol = 8,
    byrow = TRUE
  ) |> t() |> as.vector() |> chr_to_int()
  expect_identical(score_and_path$score, true_scores)
  expect_identical(score_and_path$path, true_path)
})

test_that("mergepairs c++: merge_by_path_backtrack w/ mismatch, no indel", {
  #   A  T  G  C  G  A  T
  #   40 39 37 38 35 24 18 
  #            C  G  G  T  T  A  C
  #            39 38 32 34 33 31 27
  #   A  T  G  C  G  G  T  T  A  C
  #   40 39 37 82 78  9 57 33 31 27  
  sequences <- c(forward = "ATGCGAT", reverse = "CGGTTAC")
  qualities <- c(forward = "IHFGD93", reverse = "HGACB@<")
  scoring_matrix <- create_scoring_matrix(7L, -3L)
  gap_p <- -7L
  true_path <- matrix(
    c(
      "l", "l", "l", "l", "l", "l", "l", "l",
      "u", "d", "d", "d", "d", "l", "d", "u",
      "u", "u", "d", "d", "u", "d", "l", "u",
      "u", "u", "u", "d", "d", "u", "d", "l",
      "u", "u", "d", "u", "d", "u", "u", "d",
      "u", "u", "u", "u", "u", "d", "u", "u",
      "u", "d", "u", "u", "u", "u", "d", "u",
      "u", "u", "u", "u", "d", "u", "u", "u"
    ),
    ncol = 8,
    byrow = TRUE
  ) |> t() |> as.vector() |> chr_to_int()
  alignment <- test_merge_by_path_backtrack(
    true_path,
    sequences,
    qualities,
    merged_qualities$match,
    merged_qualities$mismatch
  )
  true_qualities <- c(40L, 39L, 37L, 82L, 78L, 9L, 57L, 33L, 31L, 27L) |>
    quality_integer_to_string()
  true_alignment <- list(
    sequence = "ATGCGGTTAC",
    quality = true_qualities,
    overlap_length = 4,
    overlap_matches = 3
  )
  expect_mapequal(alignment, true_alignment)
})

test_that("mergepairs c++: find_best_scoring_path w/ mismatch, w/ indel", {
  #   A  T  G  C  G  A  T  T  A
  #   40 39 37 38 35 24 18 19 12
  #            C  G  G  T  -  A  C  G  T
  #            14 17 20 34    31 27 0  41 
  #   A  T  G  C  G  G  T  T  A  C  G  T
  #   40 39 37 82 78  9 57 19 31 27 0  41
  sequences <- c(forward = "ATGCGATTA", reverse = "CGGTACGT")
  qualities_fwd_int <- c(40L, 39L, 37L, 38L, 35L, 24L, 18L, 19L, 12L)
  qualities_rev_int <- c(14L, 17L, 20L, 34L, 31L, 27L, 0L, 41L)
  qualities <- c(
    forward = quality_integer_to_string(qualities_fwd_int),
    reverse = quality_integer_to_string(qualities_rev_int)
  )
  score_and_path <- test_find_best_scoring_path(
    sequences,
    create_scoring_matrix(7L, -3L),
    gap_p = -7L
  )
  
  true_scores <- matrix(
    c(
        0L,   0L,   0L,   0L,   0L,   0L,   0L,   0L,   0L,   0L,
       -7L,  -3L,  -3L,  -3L,   7L,   0L,  -3L,  -3L,  -3L,   0L,
      -14L, -10L,  -6L,   4L,   0L,  14L,   7L,   0L,  -6L,   0L,
      -21L, -17L, -13L,   1L,   1L,   7L,  11L,   4L,  -3L,   0L,
      -28L, -24L, -10L,  -6L,  -2L,   0L,   4L,  18L,  11L,   4L,
      -35L, -21L, -17L, -13L,  -9L,  -5L,   7L,  11L,  15L,   18L,
      -42L, -28L, -24L, -20L,  -6L, -12L,   0L,   4L,   8L,   18L,
      -49L, -35L, -31L, -17L, -13L,   1L,  -6L,  -3L,   1L,   18L,
      -56L, -42L, -28L, -24L, -20L,  -6L,  -2L,   1L,   4L,   18L
    ),
    ncol = nchar(sequences["forward"]) + 1,
    byrow = TRUE
  ) |> t() |> as.vector()
  true_path <- matrix(
    c(
      "l",  "l",  "l",  "l",  "l",  "l",  "l",  "l",  "l",  "l",  
      "u",  "d",  "d",  "d",  "d",  "l",  "d",  "d",  "d",  "u",  
      "u",  "u",  "d",  "d",  "u",  "d",  "l",  "l",  "d",  "u",  
      "u",  "u",  "u",  "d",  "d",  "u",  "d",  "l",  "l",  "u",  
      "u",  "u",  "d",  "u",  "d",  "u",  "u",  "d",  "l",  "l",  
      "u",  "d",  "u",  "u",  "u",  "d",  "d",  "u",  "d",  "d",  
      "u",  "u",  "u",  "u",  "d",  "u",  "u",  "u",  "u",  "u",  
      "u",  "u",  "u",  "d",  "u",  "d",  "l",  "u",  "u",  "u",  
      "u",  "u",  "d",  "u",  "u",  "u",  "d",  "d",  "d",  "u" 
    ),
    ncol = nchar(sequences["forward"]) + 1,
    byrow = TRUE
  )  |> t() |> as.vector() |> chr_to_int()
  expect_identical(score_and_path$score, true_scores)
  expect_identical(score_and_path$path, true_path)
})

# Add test cases for alignments:
#  ATGCATGCAG
#       TGCAG

# ATGCATGCAG
# ATGCA
# This should never occur

# Add test case with indels within the overlap

test_that("mergepairs c++: align_seqs_and_quals full alignment", {
  # ATGCGAT
  #    CGGTTAC
  # IHF
  sequences <- c(forward = "ATGCGAT", reverse = "CGGTTAC")
  qualities <- c(forward = "IHFGD93", reverse = "HGACB@<")
  alignment <- align_seqs_and_quals(
    sequences,
    qualities,
    c("match" = 7L, "mismatch" = -4L, "gap_penalty" = 7L),
    merged_qualities$match,
    merged_qualities$mismatch
  )
  true_alignment <- c(sequence = "ATGCGGTTAC", quality = "IHFso*ZB@<")
  expect_mapequal(alignment, true_alignment)
})