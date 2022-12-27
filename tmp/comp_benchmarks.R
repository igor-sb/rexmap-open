nsize = 200
nitem = 100

set.seed(42)
x <- lapply(
  1:nitem,
  function(i) {
    paste(sample(LETTERS[1:5], size = nsize, replace = TRUE), collapse = "")
  }
)
y <- lapply(
  1:nitem,
  function(i) {
    paste(sample(LETTERS[1:5], size = nsize, replace = TRUE), collapse = "")
  }
)

set.seed(42)
xstr <- lapply(
  1:nitem,
  function(i) {
    sample(LETTERS[1:5], size = nsize, replace = TRUE)
  }
)
ystr <- lapply(
  1:nitem,
  function(i) {
    sample(LETTERS[1:5], size = nsize, replace = TRUE)
  }
)

set.seed(42)
xint <- lapply(
  1:nitem,
  function(i) {
    sample(1:5, size = nsize, replace = TRUE)
  }
)
yint <- lapply(
  1:nitem,
  function(i) {
    sample(1:5, size = nsize, replace = TRUE)
  }
)

string_compare_apply <- function(s1_vec, s2_vec) {
  mapply(
    function(s1, s2) string_compare(s1, s2),
    s1_vec, s2_vec,
    SIMPLIFY = FALSE
  )
}

integer_vec_compare_apply <- function(s1_vec, s2_vec) {
  mapply(
    function(s1, s2) integer_vec_compare(s1, s2),
    s1_vec, s2_vec,
    SIMPLIFY = FALSE
  )
}

string_vec_compare_apply <- function(s1_vec, s2_vec) {
  mapply(
    function(s1, s2) string_vec_compare(s1, s2),
    s1_vec, s2_vec,
    SIMPLIFY = FALSE
  )
}

microbenchmark(
  string_compare_apply(x, y),
  string_vec_compare_apply(xstr, ystr),
  integer_vec_compare_apply(xint, yint)
)