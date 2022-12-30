chr_to_int <- Vectorize(function(chr) {
  if (is.na(chr)) return(NA_integer_)
  as.integer(charToRaw(chr))
}, vectorize.args = "chr", USE.NAMES = FALSE)

int_to_chr <- Vectorize(function(int) {
  if (is.na(int)) return("")
  rawToChar(as.raw(int))
}, vectorize.args = "int")

dp_path_matrix <- function(path, ncol) {
  matrix(int_to_chr(path), ncol = ncol, byrow = TRUE)
}

dp_score_matrix <- function(score, ncol) {
  matrix(score, ncol = ncol, byrow = TRUE)
}

quality_string_to_integer <- function(quality) {
  strsplit(quality, "") |> unname() |> unlist() |> chr_to_int() - 33L
}

quality_integer_to_string <- function(integers) {
  paste(int_to_chr(integers + 33L), collapse = "")
}