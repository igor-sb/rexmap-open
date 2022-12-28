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
  sequences <- c("ATGCGA", "CGATT")
  ncol <- nchar(sequences[1]) + 1L
  nrow <- nchar(sequences[2]) + 1L
  # C++ index 1, 2: 
  expect_identical(
    get_indexes(1, 2, ncol),
    c("diag" = 7, "left" = 14, "upper" = 8, "current" = 15)
  )
})

test_that("mergepairs c++: dp score = zero & path = left for first row", {
  sequences <- c("ATGCGA", "CGATT")
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
  sequences <- c("ATGCGA", "CGATT")
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

test_that("mergepairs c++: dp score and path for other rows and columns", {
  sequences <- c("ATGCGA", "CGATT")
  scoring_matrix <- create_scoring_matrix(5L, -5L)
  gap_p <- -7L
  ncol <- nchar(sequences[1]) + 1L
  nrow <- nchar(sequences[2]) + 1L
  score <- rep(NA_integer_, ncol * nrow)
  path <- rep(NA_integer_, ncol * nrow)
  score[1:7] <- 0L
  path[1:7] <- chr_to_int("l")
  first_col <- seq.int(1L + ncol, 1L + ncol * (nrow - 1L), by = ncol)
  score[first_col] <- seq.int(gap_p, gap_p * (nrow - 1L), by = gap_p)
  path[first_col] <- chr_to_int("u")
  score_and_path <- test_calc_score_path_other(score, path, sequences,
                                               scoring_matrix, gap_p)
  
  expect_identical(
    score_and_path$score,
    c(0L, 0L, 0L, 0L, 0L, 0L, 0L, -7L, -5L, 5L, -2L, 5L, -2L, 0L, 
      -14L, -12L, -2L, 0L, -2L, 10L, 3L, -21L, -9L, -9L, -7L, -5L, 
      3L, 5L, -28L, -16L, -14L, -14L, -12L, -4L, 5L, -35L, -23L, -21L, 
      -19L, -19L, -11L, 5L)
  )
  expect_identical(
    int_to_chr(score_and_path$path),
    c("l", "l", "l", "l", "l", "l", "l", "u", "d", "d", "l", "d", 
      "l", "u", "u", "u", "u", "d", "u", "d", "l", "u", "d", "u", "u", 
      "d", "u", "d", "u", "u", "d", "u", "u", "u", "u", "u", "u", "u", 
      "d", "u", "u", "u")
  )
})

test_that("mergepairs c++: find_best_scoring_overlap", {
  sequences <- c("ATGCGA", "CGATT")
  scoring_matrix <- create_scoring_matrix(5L, -5L)
  gap_p <- -7L
  score_and_path <- find_best_scoring_overlap(sequences, scoring_matrix, gap_p)
  expect_identical(
    score_and_path$score,
    c(0L, 0L, 0L, 0L, 0L, 0L, 0L, -7L, -5L, 5L, -2L, 5L, -2L, 0L, 
      -14L, -12L, -2L, 0L, -2L, 10L, 3L, -21L, -9L, -9L, -7L, -5L, 
      3L, 5L, -28L, -16L, -14L, -14L, -12L, -4L, 5L, -35L, -23L, -21L, 
      -19L, -19L, -11L, 5L)
  )
  expect_identical(
    int_to_chr(score_and_path$path),
    c("l", "l", "l", "l", "l", "l", "l", "u", "d", "d", "l", "d", 
      "l", "u", "u", "u", "u", "d", "u", "d", "l", "u", "d", "u", "u", 
      "d", "u", "d", "u", "u", "d", "u", "u", "u", "u", "u", "u", "u", 
      "d", "u", "u", "u")
  )
})

test_that("mergepairs c++: create alignment", {
  # ATGCGA
  #    CGATT
  sequences <- c("ATGCGA", "CGATT")
  qualities <- c("IHFGDA", "HGACB")
  scoring_matrix <- create_scoring_matrix(5L, -5L)
  gap_p <- -7L
  path <- chr_to_int(
    c("l", "l", "l", "l", "l", "l", "l", "u", "d", "d", "l", "d", 
      "l", "u", "u", "u", "u", "d", "u", "d", "l", "u", "d", "u", "u", 
      "d", "u", "d", "u", "u", "d", "u", "u", "u", "u", "u", "u", "u", 
      "d", "u", "u", "u")
  )
  alignment <- create_alignment(path, sequences, qualities)
  alignment_strings <- Reduce(rbind, alignment, right = TRUE) |>
    apply(2, paste, collapse = "")
  
  expect_identical(
    alignment_strings,
    c(seq1 = "--AGCGTA", qua2 = "BCAGH   ", seq2 = "TTAGC---", qua1 = "  ADGFHI"
    )
  )
})